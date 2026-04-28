"""
    End-to-end pipeline driver `N → List[ClassifiedMTC]`.

Given a conductor `N`, `classify_mtcs_at_conductor` executes the full pipeline and returns a list of fully classified
MTCs, each carrying:

- the SL(2, ℤ/N) stratum (m_λ) decomposition (Phase 0/1),
- F_p-validated modular data (Phase 2),
- Galois-coherent, arithmetic S-matrix in ℤ[√d] (Phase 3),
- exact `(S, T)` lifted to `Q(ζ_N)` with `N` fixed as the base conductor.

This file provides the orchestration entry points:

- `classify_from_group(group, all_primes; ...)`: given a Galois-coherent
  group from `group_mtcs_galois_aware`, performs Phase 3 CRT + Phase 4
  exact lift and returns a single `ClassifiedMTC`.

- `classify_mtcs_at_conductor(N; ...)`: the conductor-first top-level
  classification driver.

Design notes:
- The pipeline is conductor-first: `N` is the outer loop; rank emerges
  from stratum enumeration.
- Modular data is lifted exactly to `Q(ζ_N)`.  The effective conductor is
  fixed to the user-supplied `N`.
- `skip_FR = false` runs exact cyclotomic `(F, R)` reconstruction.
"""
# ============================================================
#  classify_from_group: CRT + (F, R) solve for one Galois sector
# ============================================================

function _fixed_exact_group(group::Dict{Int, MTCCandidate})
    return all(c -> c.U_params == :atomic || c.U_params == :block_invariant,
               values(group))
end

function _exact_stratum_data(catalog, stratum::Stratum)
    S_exact, T_exact, _, _ = build_block_diagonal(catalog, stratum)
    return (S = S_exact, T = T_exact)
end

function _verify_exact_candidate_at_prime(S_exact, T_exact, candidate::MTCCandidate,
                                          N::Int, p::Int)
    zeta_N_Fp = find_zeta_in_Fp(N, p)
    r = nrows(S_exact)
    r == size(candidate.S_Fp, 1) || return false
    length(T_exact) == length(candidate.T_Fp) || return false
    for i in 1:r, j in 1:r
        cyclotomic_to_Fp(S_exact[i, j], zeta_N_Fp, p) == candidate.S_Fp[i, j] ||
            return false
    end
    for i in 1:r
        cyclotomic_to_Fp(T_exact[i], zeta_N_Fp, p) == candidate.T_Fp[i] ||
            return false
    end
    return true
end

"""
    classify_from_group(group, N, stratum, all_primes;
                        quadratic_d, scale_factor = 2,
                        sqrtd_fn = compute_sqrt_d_mod_p,
                        reconstruction_bound = 50,
                        galois_sector = 1,
                        test_primes = nothing,
                        skip_FR = false,
                        verbose = false)
        -> ClassifiedMTC

Given a single Galois-coherent group from `group_mtcs_galois_aware`,
perform Phase 3 CRT reconstruction, lift `(S, T)` to `Q(ζ_N)`, and
return a `ClassifiedMTC`.

Arguments:
- `group::Dict{Int, MTCCandidate}`:  a single group (one sector)
- `N::Int`:                          conductor
- `stratum::Stratum`:                the stratum the group came from
- `all_primes::Vector{Int}`:         primes present in `group`
- `quadratic_d::Int`:                    `d` for ℤ[√d]
- `scale_factor::Int = 2`:           scalar multiplying S before recon
- `sqrtd_fn`:                        `(d, p) -> Int` returning √d mod p.
                                     For SU(2)_4 / d = 3 /quadratic_d = 3,
                                     use `compute_sqrt3_cyclotomic_mod_p`
                                     to ensure Galois-consistent choice.
- `reconstruction_bound::Int = 50`:  ℤ[√d] coefficient bound for rational
                                     reconstruction.
- `galois_sector::Int = 1`:          sector index for provenance.
- `test_primes::Vector{Int}=nothing`: which primes to use for CRT.
                                     If `nothing`, uses first
                                     `max(2, length(all_primes) ÷ 2)`
                                     primes; remaining are fresh.
- `skip_FR::Bool = false`:           when true, stop after exact
                                     cyclotomic modular-data lift.
- `verbose::Bool = false`:           print progress.
"""
function classify_from_group(group::Dict{Int, MTCCandidate},
                             N::Int,
                             stratum::Stratum,
                             all_primes::Vector{Int};
                             N_input::Int = N,
                             quadratic_d::Int = _default_quadratic_d_for_conductor(N),
                             scale_factor::Int = 2,
                             sqrtd_fn = compute_sqrt_d_mod_p,
                             reconstruction_bound::Int = 50,
                             galois_sector::Int = 1,
                             test_primes::Union{Vector{Int}, Nothing} = nothing,
                             catalog = nothing,
                             skip_FR::Bool = false,
                             verbose::Bool = false,
                             kwargs...)
    isempty(kwargs) || error("unsupported keyword arguments: $(collect(keys(kwargs)))")
    group_primes = sort(collect(keys(group)))

    if test_primes === nothing
        half = max(2, length(group_primes) ÷ 2)
        used = group_primes[1:half]
        fresh = group_primes[(half + 1):end]
    else
        used = sort([p for p in test_primes if haskey(group, p)])
        fresh = sort([p for p in group_primes if !(p in used)])
    end
    length(used) >= 2 || error(
        "need ≥ 2 used primes for CRT; got $(length(used)). " *
        "Group has $(length(group_primes)) primes; pass test_primes " *
        "explicitly or add more primes.")

    verbose && println("  sector=$galois_sector: used=$used, fresh=$fresh")

    exact_data = catalog === nothing || !_fixed_exact_group(group) ?
        nothing :
        _exact_stratum_data(catalog, stratum)

    # -------- Phase 3: CRT in ℤ[√d] --------
    used_subgroup = Dict(p => group[p] for p in used)
    recon_S = reconstruct_S_matrix(used_subgroup;
                                    quadratic_d = quadratic_d,
                                    bound = reconstruction_bound,
                                    sqrtd_fn = sqrtd_fn)

    # Verify at fresh primes (if any)
    verify_fresh = isempty(fresh)  # vacuously true if none
    if !isempty(fresh)
        all_ok = true
        for p in fresh
            ok = verify_reconstruction(recon_S, group[p], quadratic_d;
                                        scale = scale_factor,
                                        sqrtd_fn = sqrtd_fn)
            all_ok &= ok
            verbose && println("    fresh p=$p: $(ok ? "✓" : "✗")")
        end
        verify_fresh = all_ok
    end
    if exact_data !== nothing
        exact_ok = all(p -> _verify_exact_candidate_at_prime(exact_data.S, exact_data.T,
                                                             group[p], N, p),
                       group_primes)
        verify_fresh = exact_ok
        verbose && !exact_ok && println("    exact fixed-stratum lift failed finite-field verification")
    elseif !verify_fresh
        error("CRT S reconstruction failed fresh-prime verification; " *
              "the candidate is not stable in the selected reconstruction field")
    end

    # -------- Phase 4 input: exact lift (S, T) to Q(ζ_N) --------
    rep = first(values(group))
    if exact_data === nothing
        zeta_Fp = find_zeta_in_Fp(N, rep.p)
        S_cyc, T_cyc, Nijk = lift_mtc_candidate(rep, recon_S;
                                                d = quadratic_d,
                                                N = N,
                                                zeta_Fp = zeta_Fp,
                                                scale = scale_factor)
    else
        S_cyc = exact_data.S
        T_cyc = exact_data.T
        Nijk = rep.N
    end
    recon_S_phase4 = recon_S

    # TensorCategories assumes the unit object is at index 1.
    # MTCCandidate stores `unit_index` explicitly and may keep a different
    # basis ordering; permute all lifted data coherently before Phase 4.
    if rep.unit_index != 1
        perm = vcat(rep.unit_index, [i for i in 1:length(T_cyc) if i != rep.unit_index])
        S_cyc = S_cyc[perm, perm]
        T_cyc = T_cyc[perm]
        Nijk = Nijk[perm, perm, perm]
        recon_S_phase4 = recon_S[perm, perm]
    end
    T_for_phase4 = _normalize_T_by_unit(T_cyc)

    # -------- Exact (F, R) layer --------
    rank = size(Nijk, 1)

    if skip_FR
        return ClassifiedMTC(N, N_input, rank, stratum, Nijk,
                             recon_S_phase4, quadratic_d, scale_factor,
                             used, fresh, verify_fresh,
                             S_cyc, T_for_phase4, nothing, nothing, nothing,
                             galois_sector)
    end
    verbose && println("  exact (F,R) layer requested on rank=$rank...")

    fr_result = compute_FR_from_ST(Nijk;
                                   context = CyclotomicContext(N),
                                   S = S_cyc,
                                   T = T_for_phase4,
                                   return_all = true,
                                   primes = all_primes,
                                   verbose = verbose)
    fr_result.F === nothing && error("Phase 4 could not produce any (F,R) solution")
    isempty(fr_result.candidates) && error("Phase 4 produced no valid (F,R) candidates")

    # OBSOLETE: candidate loop with direct `_modular_data_roundtrip`
    # calls is replaced by Phase 5 API `_select_fr_for_st`.
    selection = _select_fr_for_st(fr_result.candidates, Nijk, S_cyc, T_for_phase4, N)
    selected = selection.selected
    md_roundtrip = _modular_data_roundtrip(selected.F, selected.R, Nijk,
                                           selection.score.S_roundtrip,
                                           selection.score.T_roundtrip, N)
    verbose && println("  modular-data roundtrip: perm=$(md_roundtrip.best_perm), " *
                       "S_err=$(md_roundtrip.S_max), T_err=$(md_roundtrip.T_max), " *
                       "ok=$(md_roundtrip.ok)")

    return ClassifiedMTC(N, N_input, rank, stratum, Nijk, recon_S_phase4,
                         quadratic_d, scale_factor,
                         used, fresh, verify_fresh,
                         md_roundtrip.S_roundtrip, md_roundtrip.T_roundtrip,
                         selected.F, selected.R, md_roundtrip,
                         galois_sector)
end

# ============================================================
#  classify_mtcs_at_conductor: top-level driver
# ============================================================

"""
    classify_mtcs_at_conductor(N; max_rank = 5, primes = nothing, strata = nothing,
                                quadratic_d = 3, scale_factor = 2,
                                conductor_mode = :full_mtc,
                                min_primes = 4, prime_start = 29, prime_window = 2000,
                                sqrtd_fn = nothing,
                                verlinde_threshold = 3,
                                max_block_dim = 3,
                                reconstruction_bound = 50,
                                skip_FR = false,
                                verbose = true)
        -> Vector{ClassifiedMTC}

Fully automated MTC classification for a given conductor input `N`.
For a user-friendly wrapper that auto-selects parameters, see
`classify_mtcs_auto`.

Pipeline:
  Phase 0: build SL(2, ℤ/N) atomic irrep catalog (≤ max_rank)
  Phase 1: enumerate strata (m_λ) with Σ m_λ d_λ = r for each r ≤ max_rank
  Phase 2: for each stratum, sweep block-U at each prime → MTC candidates
  Phase 3: Galois-aware grouping + CRT reconstruction in ℤ[√d]
  Phase 4: lift to Q(ζ_N), then compute exact `(F, R)`

Returns one `ClassifiedMTC` per Galois sector per stratum that yields a
valid modular-data candidate and exact cyclotomic `(F, R)` data.

Note on conductor: `N` is the cyclotomic base conductor. Internally,
`N_effective = N`; the conductor is not enlarged after discovering an
S-field.

Arguments:
- `N::Int`:                        cyclotomic base conductor. Internal
                                   search uses `N_effective = N`.
- `max_rank::Int = 5`:             maximum rank to consider
- `primes::Union{Nothing, Vector{Int}} = nothing`:
                                   good primes (must satisfy `N_effective | p-1`).
                                   If `nothing`, automatically selected via
                                   `select_admissible_primes(N_effective; ...)`.
                                   Split automatically into "used" (first
                                   half) and "fresh" (second half) for
                                   CRT + cross-validation. Minimum 4
                                   recommended (2 used + 2 fresh).
- `strata::Union{Nothing, Vector{Stratum}}`:
                                   if `nothing`, enumerates all strata
                                   of rank ≤ max_rank via
                                   `enumerate_strata`. Otherwise uses
                                   the given strata directly (faster for
                                   targeted re-runs).
- `quadratic_d::Int = 3`:              ℤ[√d] the S-matrix is expected to
                                   live in. Must be chosen to match the
                                   MTC's quadratic field: e.g. 3 for
                                   SU(2)_4 / Ising, 5 for MTCs involving
                                   the golden ratio, 2 for √2-pointed
                                   MTCs. A wrong `quadratic_d` will silently
                                   return 0 classified MTCs.
- `scale_factor::Int = 2`:         scalar multiplying `S` at
                                   reconstruction (matches
                                   `reconstruct_S_matrix` convention).
- `conductor_mode::Symbol = :full_mtc`:
                                   compatibility option. Only
                                   `:full_mtc` is supported, and it fixes
                                   `N_effective = N`.
- `min_primes::Int = 4`:           when `primes = nothing`, number of
                                   admissible primes to auto-select.
- `prime_start::Int = 29`:         when `primes = nothing`, start of
                                   prime search range (exclusive).
- `prime_window::Int = 2000`:      when `primes = nothing`, scan width.
- `sqrtd_fn`:                      custom √d-in-F_p function. If
                                   `nothing` (default), chooses:
                                   cyclotomic variant for `quadratic_d ∈
                                   {2,3,5}` and anchored mode for
                                   other `quadratic_d` (anchor prime +
                                   branch-transform layer over raw
                                   Tonelli roots):
                                   `compute_sqrt3_cyclotomic_mod_p`
                                   (needs 24 | p-1) for `quadratic_d = 3`,
                                   `compute_sqrt2_cyclotomic_mod_p`
                                   (needs 24 | p-1) for `quadratic_d = 2`,
                                   `compute_sqrt5_cyclotomic_mod_p`
                                   (needs 5 | p-1) for `quadratic_d = 5`,
                                   else anchored mode (with per-prime
                                   sign alignment and pre-group
                                   contradiction check/split).
- `verlinde_threshold::Int = 3`:   max absolute fusion coefficient
                                   allowed (beyond which the candidate
                                   is rejected). See
                                   `find_mtcs_at_prime`.
- `max_block_dim::Int = 3`:        cap on the degenerate T-eigenspace
                                   dimension for block-U solving.
- `search_mode::Symbol = :groebner`:
                                   Phase 2 block-U backend mode passed
                                   to `find_mtcs_at_prime` (`:groebner`
                                   or `:exhaustive`).
- `max_units_for_groebner::Int = typemax(Int)`:
                                   cap on fixed-unit Gröbner systems
                                   tried per stratum/prime in Phase 2.
- `groebner_allow_fallback::Bool = false`:
                                   when `search_mode=:groebner`,
                                   whether to fall back to explicit
                                   Cayley/reflection enumeration if solver
                                   extraction is empty.
- `precheck_unit_axiom::Bool = true`:
                                   run fast unit-axiom precheck before
                                   full Verlinde tensor evaluation in
                                   Phase 2 candidate loop.
- `reconstruction_bound::Int = 50`: coefficient bound for ℤ[√d]
                                   rational reconstruction.
- `skip_FR::Bool = false`:         when true, stop after exact
                                   cyclotomic modular-data lift.
- `verbose::Bool = true`:          print per-phase progress.
"""
function classify_mtcs_at_conductor(N::Int;
                                    max_rank::Int = 5,
                                    primes::Union{Nothing, Vector{Int}} = nothing,
                                    strata::Union{Nothing, Vector{Stratum}} = nothing,
                                    quadratic_d::Union{Int, Nothing} = nothing,
                                    scale_factor::Int = 2,
                                    conductor_mode::Symbol = :full_mtc,
                                    min_primes::Int = 4,
                                    prime_start::Int = 29,
                                    prime_window::Int = 2000,
                                    sqrtd_fn = nothing,
                                    verlinde_threshold::Int = 3,
                                    max_block_dim::Int = 3,
                                    search_mode::Symbol = :groebner,
                                    max_units_for_groebner::Int = typemax(Int),
                                    groebner_allow_fallback::Bool = false,
                                    precheck_unit_axiom::Bool = true,
                                    reconstruction_bound::Int = 50,
                                    skip_FR::Bool = false,
                                    verbose::Bool = true,
                                    kwargs...)
    user_sqrtd_fn = sqrtd_fn
    N_effective = compute_effective_conductor(N; conductor_mode = conductor_mode)
    isempty(kwargs) || error("unsupported keyword arguments: $(collect(keys(kwargs)))")
    quadratic_d = quadratic_d === nothing ?
        _default_quadratic_d_for_conductor(N_effective) :
        quadratic_d

    prime_requirement = _prime_requirement_for_reconstruction(N_effective, quadratic_d, sqrtd_fn)
    chosen_primes = primes === nothing ?
        select_admissible_primes(prime_requirement;
                                 min_count = min_primes,
                                 start_from = prime_start,
                                 window = prime_window) :
        copy(primes)

    verbose && println("Conductor mode: $conductor_mode " *
                       "(input N=$N, N_effective=$N_effective)")
    verbose && primes === nothing &&
        println("PrimeSelection: auto-selected primes = $chosen_primes")
    verbose && primes === nothing && prime_requirement != N_effective &&
        println("PrimeSelection: using reconstruction requirement $prime_requirement " *
                "(N_effective=$N_effective, d=$quadratic_d)")

    # Prime validity check
    for p in chosen_primes
        (p - 1) % N_effective == 0 || error(
            "prime $p does not satisfy N_effective | p-1 " *
            "(input N=$N; N_effective=$N_effective)")
    end
    length(chosen_primes) >= 2 || error(
        "need at least 2 primes (got $(length(chosen_primes)))")

    # ------- Phase 0: atomic catalog -------
    # Enumerate atomic SL(2, Z/N)-irreps at the user-specified base conductor N.
    # The field-extension effect is handled later via N_effective in arithmetic
    # checks (prime admissibility, sqrt(d)-dependent reconstruction).
    verbose && println("=== Phase 0: atomic SL(2, ℤ/$N) irrep catalog ===")
    catalog = build_atomic_catalog(N; max_rank = max_rank, verbose = false)
    verbose && println("  $(length(catalog)) atomic irreps (≤ rank $max_rank)")

    # ------- Phase 1: strata -------
    verbose && println("\n=== Phase 1: stratum enumeration ===")
    strata_list = strata === nothing ?
        vcat([enumerate_strata(catalog, r) for r in 1:max_rank]...) :
        strata
    verbose && println("  $(length(strata_list)) strata")

    # ------- Phase 2: find MTCs at each prime, for each stratum -------
    verbose && println("\n=== Phase 2: block-U sweep at primes $chosen_primes ===")

    # Collect (stratum, Dict{p => candidates}) only for strata that yield
    # something at at least one prime.
    stratum_results = Vector{Tuple{Stratum, Dict{Int, Vector{MTCCandidate}}}}()
    skipped_reasons = Dict{String, Int}()
    for (si, st) in enumerate(strata_list)
        results_by_prime = Dict{Int, Vector{MTCCandidate}}()
        any_found = false
        any_prime_errored = false
        last_err = ""
        for p in chosen_primes
            local cands
            try
                cands = find_mtcs_at_prime(catalog, st, p;
                                            verlinde_threshold = verlinde_threshold,
                                            max_block_dim = max_block_dim,
                                            search_mode = search_mode,
                                            max_units_for_groebner = max_units_for_groebner,
                                            groebner_allow_fallback = groebner_allow_fallback,
                                            precheck_unit_axiom = precheck_unit_axiom)
            catch err
                # e.g. max_block_dim exceeded, or p not admissible
                # for this stratum. Record the first such error for
                # a terse verbose summary at end.
                any_prime_errored = true
                last_err = sprint(showerror, err)
                continue
            end
            if !isempty(cands)
                results_by_prime[p] = cands
                any_found = true
            end
        end
        if any_found
            push!(stratum_results, (st, results_by_prime))
            verbose && println("  stratum $si ($(describe_stratum(st, catalog))): " *
                               "MTCs at $(length(results_by_prime)) primes")
        elseif any_prime_errored && verbose
            # Summarise the error kind (first 80 chars) to help diagnose
            # repeated patterns without flooding output.
            key = first(last_err, 80)
            skipped_reasons[key] = get(skipped_reasons, key, 0) + 1
        end
    end
    verbose && println("  $(length(stratum_results)) strata yielded candidates; " *
                       "$(length(strata_list) - length(stratum_results)) skipped")
    if verbose && !isempty(skipped_reasons)
        println("  skip reasons (first 80 chars of error × count):")
        for (reason, count) in sort(collect(skipped_reasons), by = x -> -x[2])
            println("    [$count×] $reason")
        end
    end

    # ------- Phase 3: per-stratum, per-sector modular-data reconstruction -------
    verbose && println("\n=== Phase 3: CRT modular-data reconstruction ===")

    out = ClassifiedMTC[]
    for (st, results_by_prime) in stratum_results
        # Need candidates at ≥ 2 primes to do CRT at all.
        present_primes = sort(collect(keys(results_by_prime)))
        if length(present_primes) < 2
            verbose && println("  stratum $(describe_stratum(st, catalog)): " *
                               "only $(length(present_primes)) prime(s), skip")
            continue
        end

        # Build d-wise sqrt branch selector:
        #  - d = 2,3,5: cyclotomic
        #  - otherwise: anchored transform layer (not raw Tonelli direct)
        anchor = present_primes[1]
        selector_mode = :custom
        active_sqrtd_fn = user_sqrtd_fn
        branch_sign_getter = nothing
        branch_sign_setter = nothing
        if active_sqrtd_fn === nothing
            selector = build_sqrtd_selector(quadratic_d, present_primes, anchor; verbose = verbose)
            active_sqrtd_fn = selector.sqrtd_fn
            selector_mode = selector.mode
            branch_sign_getter = selector.branch_sign_getter
            branch_sign_setter = selector.branch_sign_setter
        end
        verbose && println("  d=$quadratic_d, sqrt-branch mode=$(selector_mode == :custom ? "custom" : String(selector_mode)), grouping bound=$reconstruction_bound")

        contradictory = _branch_consistency_precheck(results_by_prime, anchor, quadratic_d, active_sqrtd_fn;
                                                     reconstruction_bound = reconstruction_bound,
                                                     branch_sign_getter = branch_sign_getter,
                                                     branch_sign_setter = branch_sign_setter,
                                                     verbose = verbose)
        # If contradictory primes remain, split sector inputs (anchor+main, anchor+each contradictory).
        results_chunks = Dict{Int, Vector{MTCCandidate}}[]
        if isempty(contradictory)
            push!(results_chunks, results_by_prime)
        else
            main_dict = Dict{Int, Vector{MTCCandidate}}()
            for (p, cands) in results_by_prime
                if !(p in contradictory)
                    main_dict[p] = cands
                end
            end
            if length(main_dict) >= 2
                push!(results_chunks, main_dict)
            end
            for p in contradictory
                local_dict = Dict{Int, Vector{MTCCandidate}}(
                    anchor => results_by_prime[anchor],
                    p => results_by_prime[p])
                push!(results_chunks, local_dict)
            end
        end

        # Galois-aware grouping
        local groups
        groups = Vector{Dict{Int, MTCCandidate}}()
        grouping_failed = false
        for chunk in results_chunks
            try
                append!(groups, group_mtcs_galois_aware(chunk, anchor;
                                                        quadratic_d = quadratic_d,
                                                        reconstruction_bound = reconstruction_bound,
                                                        sqrtd_fn = active_sqrtd_fn,
                                                        branch_sign_getter = branch_sign_getter,
                                                        branch_sign_setter = branch_sign_setter))
            catch err
                grouping_failed = true
                verbose && println("  stratum: galois grouping failed: $err")
                break
            end
        end
        if grouping_failed
            continue
        end

        verbose && println("  stratum (rank $(st.total_dim)): " *
                           "$(length(groups)) Galois sectors")

        # used/fresh split: first half / second half of primes in this group
        for (gi, group) in enumerate(groups)
            group_primes = sort(collect(keys(group)))
            if length(group_primes) < 2
                verbose && println("    sector $gi: only $(length(group_primes)) " *
                                   "prime(s) — Galois-aligned CRT needs ≥ 2. skip")
                continue
            end
            half = max(2, length(group_primes) ÷ 2)
            used = copy(group_primes[1:half])
            extra_idx = half + 1
            local cmtc
            local sector_ok = false
            local last_err_msg = ""
            while true
                try
                    cmtc = classify_from_group(group, N_effective, st, group_primes;
                                                N_input = N,
                                                quadratic_d = quadratic_d,
                                                scale_factor = scale_factor,
                                                sqrtd_fn = active_sqrtd_fn,
                                                reconstruction_bound = reconstruction_bound,
                                                galois_sector = gi,
                                                test_primes = used,
                                                catalog = catalog,
                                                skip_FR = true,
                                                verbose = verbose)
                    sector_ok = true
                    break
                catch err
                    last_err_msg = sprint(showerror, err)
                    unstable = _is_reconstruction_unstable_message(last_err_msg)
                    if unstable && extra_idx <= length(group_primes)
                        p_new = group_primes[extra_idx]
                        push!(used, p_new)
                        extra_idx += 1
                        verbose && println("    sector $gi: reconstruction unstable; " *
                                           "retry with additional prime $p_new (used=$used)")
                        continue
                    end
                    break
                end
            end
            if !sector_ok
                verbose && println("    sector $gi: classify_from_group failed: $last_err_msg")
                continue
            end
            push!(out, cmtc)
            verbose && println("    → $cmtc")
        end
    end

    if skip_FR
        verbose && println("\n=== Phase 4 skipped: returning modular-data results only ===")
        verbose && println("\n=== Done: $(length(out)) ClassifiedMTC(s) ===")
        return out
    end

    # ------- Exact (F,R) layer -------
    verbose && println("\n=== Exact (F, R) layer requested ===")
    grouped = _classify_modular_data_by_fusion_rule(out)
    key_to_members = Dict{String, Vector{Tuple{Any, Any}}}()
    for (key, idxs) in grouped
        key_to_members[key] = [(out[i].S_cyclotomic, out[i].T_cyclotomic) for i in idxs]
    end
    verbose && println("  modular-data groups: $(length(out)) (S,T) results → " *
                       "$(length(grouped)) fusion-rule keys")

    for (key, idxs) in grouped
        rep_idx = idxs[1]
        rep = out[rep_idx]
        verbose && println("  key=$key: members=$(length(idxs)), rank=$(rep.rank)")
        local fr_result
        try
            fr_result = compute_FR_from_ST(rep.Nijk;
                                           context = CyclotomicContext(rep.N),
                                           S = rep.S_cyclotomic,
                                           T = rep.T_cyclotomic,
                                           return_all = true,
                                           primes = chosen_primes,
                                           verbose = verbose)
        catch err
            msg = sprint(showerror, err)
            if occursin("exact Phase 4 did not find", msg)
                verbose && println("    Phase 4 found no (F,R) solution; leaving members without F/R")
                continue
            end
            rethrow()
        end
        if fr_result.F === nothing || isempty(fr_result.candidates)
            verbose && println("    Phase 4 produced no valid (F,R) candidates; leaving members without F/R")
            continue
        end

        # OBSOLETE: old behavior selected using only the fusion-rule
        # representative (`rep`). New behavior assigns `(F,R)` per `(S,T)`
        # member within the same fusion-rule group.
        for i in idxs
            mtc_i = out[i]
            selection = _select_fr_for_st(fr_result.candidates, mtc_i.Nijk,
                                          mtc_i.S_cyclotomic, mtc_i.T_cyclotomic, mtc_i.N)
            selected = selection.selected
            best_md = _modular_data_roundtrip(selected.F, selected.R, mtc_i.Nijk,
                                              selection.score.S_roundtrip,
                                              selection.score.T_roundtrip, mtc_i.N)
            verbose && println("    member[$i] branch: cand=$(selection.selected_index), " *
                               "galois=$(selection.score.galois_exponent), " *
                               "perm=$(best_md.best_perm), " *
                               "S_err=$(best_md.S_max), T_err=$(best_md.T_max), ok=$(best_md.ok)")
            out[i] = _with_fr_result(out[i], selected.F, selected.R, best_md;
                                     S = best_md.S_roundtrip,
                                     T = best_md.T_roundtrip)
        end
    end

    verbose && println("\n=== Done: $(length(out)) ClassifiedMTC(s) ===")
    return out
end
