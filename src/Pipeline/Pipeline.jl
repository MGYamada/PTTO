"""
End-to-end conductor-first classification pipeline.

Given a conductor `N`, ACMG searches for modular data whose modular
representation factors through `SL(2, ℤ/N)` and whose entries are reconstructed
over `Q(ζ_N)`.  A successful result records:

- the semisimple `SL(2, ℤ/N)` stratum `(m_λ)` used in the search,
- finite-field modular-data candidates at admissible split primes,
- exact cyclotomic `(S, T)` data reconstructed and cross-checked from those
  residues,
- optional exact `(F, R)` data solving pentagon and hexagon constraints.

The pipeline is conductor-first: `N` is the outer arithmetic parameter, while
rank is bounded by `max_rank` and discovered through stratum enumeration.
"""
# ============================================================
#  classify_from_group: CRT + (F, R) solve for one Galois sector
# ============================================================

const _PIPELINE_CYCLOTOMIC_RECONSTRUCTION_BOUND_CAP = 4

function _pipeline_cyclotomic_reconstruction_bounds(reconstruction_bound::Int)
    reconstruction_bound >= 1 || error("reconstruction_bound must be positive")
    capped = min(reconstruction_bound, _PIPELINE_CYCLOTOMIC_RECONSTRUCTION_BOUND_CAP)
    return (coeff_bound = capped, denominator_bound = capped)
end

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

function _reconstruct_cyclotomic_S_matrix(group::Dict{Int, MTCCandidate},
                                          N::Int;
                                          scale::Int = 2,
                                          reconstruction_bound::Int = 4,
                                          denominator_bound::Int = 4)
    primes = sort(collect(keys(group)))
    first_candidate = first(values(group))
    r = size(first_candidate.S_Fp, 1)
    ctx = CyclotomicContext(N)
    K = field(ctx)
    S = zero_matrix(K, r, r)
    for i in 1:r, j in 1:r
        values_ij = Dict(p => mod(scale * group[p].S_Fp[i, j], p) for p in primes)
        x = reconstruct_cyclotomic_element_from_residues(values_ij, ctx;
                                                         coeff_bound = reconstruction_bound,
                                                         denominator_bound = denominator_bound)
        x === nothing && error("failed to cyclotomic-CRT reconstruct S entry ($i, $j)")
        S[i, j] = x / K(scale)
    end
    return S
end

function _verify_cyclotomic_candidate_at_prime(S_exact, T_exact, candidate::MTCCandidate,
                                               N::Int, p::Int)
    return _verify_exact_candidate_at_prime(S_exact, T_exact, candidate, N, p)
end

function _groups_for_stratum(results_by_prime::Dict{Int, Vector{MTCCandidate}},
                             catalog,
                             stratum::Stratum,
                             N::Int)
    exact_data = catalog === nothing ? nothing : _exact_stratum_data(catalog, stratum)
    if exact_data !== nothing
        group = Dict{Int, MTCCandidate}()
        for (p, cands) in results_by_prime
            for c in cands
                if _verify_exact_candidate_at_prime(exact_data.S, exact_data.T, c, N, p)
                    group[p] = c
                    break
                end
            end
        end
        length(group) >= 2 && return [group]
    end
    return group_mtcs_by_fusion(results_by_prime)
end

"""
    classify_from_group(group, N, stratum, all_primes;
                        scale_factor = 2,
                        reconstruction_bound = 50,
                        galois_sector = 1,
                        test_primes = nothing,
                        skip_FR = false,
                        verbose = false)
        -> ClassifiedMTC

Given one prime-indexed Galois sector of finite-field candidates, perform
cyclotomic reconstruction of modular data over `Q(ζ_N)`, verify it at fresh
primes when available, and return a `ClassifiedMTC`.

This is mainly a lower-level orchestration hook for the conductor pipeline.
Most users should call `classify_mtcs_at_conductor` or `classify_mtcs_auto`.

Arguments:
- `group::Dict{Int, MTCCandidate}`:  a single group (one sector)
- `N::Int`:                          conductor
- `stratum::Stratum`:                the stratum the group came from
- `all_primes::Vector{Int}`:         primes present in `group`
- `scale_factor::Int = 2`:           scalar multiplying S before rational CRT
- `reconstruction_bound::Int = 50`:  requested coefficient/denominator bound
                                     for cyclotomic CRT fallback. Pipeline
                                     fallback caps this to a small search-safe
                                     bound before MITM reconstruction.
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
                             scale_factor::Int = 2,
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

    used_subgroup = Dict(p => group[p] for p in used)

    rep = first(values(group))
    verify_exact_lift = nothing
    if exact_data === nothing
        reconstruction_bounds =
            _pipeline_cyclotomic_reconstruction_bounds(reconstruction_bound)
        if reconstruction_bounds.coeff_bound < reconstruction_bound
            verbose && println("    cyclotomic CRT bound capped: " *
                               "$reconstruction_bound -> " *
                               "$(reconstruction_bounds.coeff_bound)")
        end
        S_cyc = _reconstruct_cyclotomic_S_matrix(used_subgroup, N;
                                                 scale = scale_factor,
                                                 reconstruction_bound = reconstruction_bounds.coeff_bound,
                                                 denominator_bound = reconstruction_bounds.denominator_bound)
        zeta_Fp = find_zeta_in_Fp(N, rep.p)
        T_cyc = lift_T_Fp_to_cyclotomic(rep.T_Fp, N, rep.p, zeta_Fp)
        Nijk = rep.N
    else
        exact_ok = all(p -> _verify_exact_candidate_at_prime(exact_data.S, exact_data.T,
                                                             group[p], N, p),
                       group_primes)
        verify_exact_lift = exact_ok
        if !exact_ok
            verbose && println("    exact fixed-stratum lift failed finite-field verification")
            error("exact fixed-stratum lift failed finite-field verification")
        end
        S_cyc = exact_data.S
        T_cyc = exact_data.T
        Nijk = rep.N
    end

    # Verify at fresh primes (if any)
    verify_fresh = isempty(fresh)  # vacuously true if none
    if !isempty(fresh)
        all_ok = true
        for p in fresh
            ok = _verify_cyclotomic_candidate_at_prime(S_cyc, T_cyc, group[p], N, p)
            all_ok &= ok
            verbose && println("    fresh p=$p: $(ok ? "✓" : "✗")")
        end
        verify_fresh = all_ok
    end
    verify_fresh || error("CRT modular-data reconstruction failed fresh-prime verification")

    # TensorCategories assumes the unit object is at index 1.
    # MTCCandidate stores `unit_index` explicitly and may keep a different
    # basis ordering; permute all lifted data coherently before F/R solving.
    if rep.unit_index != 1
        perm = vcat(rep.unit_index, [i for i in 1:length(T_cyc) if i != rep.unit_index])
        S_cyc = S_cyc[perm, perm]
        T_cyc = T_cyc[perm]
        Nijk = Nijk[perm, perm, perm]
    end
    T_for_fr = _normalize_T_by_unit(T_cyc)
    exact_mtc_check = validate_exact_mtc(nothing, nothing, S_cyc, T_for_fr, Nijk)
    exact_mtc_check.valid ||
        error("exact modular data failed fusion/twist validation: $(exact_mtc_check.reason)")

    # -------- Exact (F, R) layer --------
    rank = size(Nijk, 1)

    if skip_FR
        return ClassifiedMTC(N, N_input, rank, stratum, Nijk,
                             scale_factor, used, fresh, verify_fresh,
                             verify_exact_lift,
                             S_cyc, T_for_fr, nothing, nothing, nothing,
                             galois_sector)
    end
    verbose && println("  exact (F,R) layer requested on rank=$rank...")

    fr_result = compute_FR_from_ST(Nijk;
                                   context = CyclotomicContext(N),
                                   S = S_cyc,
                                   T = T_for_fr,
                                   return_all = true,
                                   primes = all_primes,
                                   verbose = verbose)
    fr_result.F === nothing && error("exact F/R reconstruction could not produce any solution")
    isempty(fr_result.candidates) && error("exact F/R reconstruction produced no valid candidates")

    # Select an exact (F,R) branch against the reconstructed modular data itself.
    selection = _select_fr_for_st(fr_result.candidates, Nijk, S_cyc, T_for_fr, N)
    selected = selection.selected
    md_roundtrip = _modular_data_roundtrip(selected.F, selected.R, Nijk,
                                           S_cyc, T_for_fr, N)
    md_roundtrip = _with_fr_roundtrip_metadata(md_roundtrip;
                                               candidate_index = selection.selected_index,
                                               galois_exponent = selection.score.galois_exponent)
    _fr_roundtrip_attachable(md_roundtrip) ||
        error("selected exact (F,R) data does not roundtrip to exact modular data")
    verbose && println("  modular-data roundtrip: perm=$(md_roundtrip.best_perm), " *
                       "S_err=$(md_roundtrip.S_max), T_err=$(md_roundtrip.T_max), " *
                       "ok=$(md_roundtrip.ok)")

    return ClassifiedMTC(N, N_input, rank, stratum, Nijk,
                         scale_factor, used, fresh, verify_fresh,
                         verify_exact_lift,
                         S_cyc, T_for_fr,
                         selected.F, selected.R, md_roundtrip,
                         galois_sector)
end

# ============================================================
#  classify_mtcs_at_conductor: top-level driver
# ============================================================

"""
    classify_mtcs_at_conductor(N; max_rank = 5, primes = nothing, strata = nothing,
                                scale_factor = 2,
                                conductor_mode = :full_mtc,
                                min_primes = 4, prime_start = 29, prime_window = 2000,
                                verlinde_threshold = 3,
                                max_block_dim = 3,
                                reconstruction_bound = 50,
                                skip_FR = false,
                                verbose = true)
        -> Vector{ClassifiedMTC}

Classify modular-data candidates at a fixed conductor `N`.

For a user-friendly wrapper that auto-selects parameters, see
`classify_mtcs_auto`.

The mathematical search space is organized as follows.  Let
`G_N = SL(2, ℤ/N)` and let the atomic catalog be the irreducible
`G_N`-representations realized over `Q(ζ_N)`.  For each rank up to
`max_rank`, ACMG enumerates semisimple strata
`ρ ≅ ⊕_λ ρ_λ^{⊕m_λ}` with total dimension equal to the rank.  Within each
stratum, Block-U search changes basis inside degenerate `T`-eigenspaces and
tests the modular-data axioms, including Verlinde integrality, over admissible
finite fields.

Pipeline:
  1. build the `SL(2, ℤ/N)` atomic irrep catalog (≤ max_rank)
  2. enumerate strata `(m_λ)` with `Σ m_λ d_λ = r` for each `r ≤ max_rank`
  3. sweep Block-U at each prime for each stratum
  4. group candidates and reconstruct modular data in `Q(ζ_N)`
  5. optionally solve exact `(F, R)` pentagon/hexagon data

Returns one `ClassifiedMTC` per Galois sector per stratum that yields a
valid exact cyclotomic modular-data candidate.  If `skip_FR=false`, the result
also attempts exact `(F, R)` reconstruction; the outcome is recorded in
`fr_status(result)`.

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
- `scale_factor::Int = 2`:         scalar multiplying `S` at
                                   cyclotomic CRT reconstruction.
- `conductor_mode::Symbol = :full_mtc`:
                                   compatibility option. Only
                                   `:full_mtc` is supported, and it fixes
                                   `N_effective = N`.
- `min_primes::Int = 4`:           when `primes = nothing`, number of
                                   admissible primes to auto-select.
- `prime_start::Int = 29`:         when `primes = nothing`, start of
                                   prime search range (exclusive).
- `prime_window::Int = 2000`:      when `primes = nothing`, scan width.
- `verlinde_threshold::Int = 3`:   max absolute fusion coefficient
                                   allowed (beyond which the candidate
                                   is rejected). See
                                   `find_mtcs_at_prime`.
- `max_block_dim::Int = 3`:        cap on the degenerate T-eigenspace
                                   dimension for block-U solving.
- `search_mode::Symbol = :groebner`:
                                   Block-U backend mode passed
                                   to `find_mtcs_at_prime` (`:groebner`
                                   or `:exhaustive`).
- `max_units_for_groebner::Int = typemax(Int)`:
                                   cap on fixed-unit Gröbner systems
                                   tried per stratum/prime in Block-U search.
- `groebner_allow_fallback::Bool = false`:
                                   when `search_mode=:groebner`,
                                   whether to fall back to explicit
                                   Cayley/reflection enumeration if solver
                                   extraction is empty.
- `precheck_unit_axiom::Bool = true`:
                                   run fast unit-axiom precheck before
                                   full Verlinde tensor evaluation in
                                   finite-field candidate loop.
- `reconstruction_bound::Int = 50`: requested coefficient and denominator
                                   bound for cyclotomic power-basis CRT.
                                   Pipeline fallback caps this to a small
                                   search-safe bound before reconstruction.
- `skip_FR::Bool = false`:         when true, stop after exact
                                   cyclotomic modular-data lift.
- `verbose::Bool = true`:          print progress.
"""
function classify_mtcs_at_conductor(N::Int;
                                    max_rank::Int = 5,
                                    primes::Union{Nothing, Vector{Int}} = nothing,
                                    strata::Union{Nothing, Vector{Stratum}} = nothing,
                                    scale_factor::Int = 2,
                                    conductor_mode::Symbol = :full_mtc,
                                    min_primes::Int = 4,
                                    prime_start::Int = 29,
                                    prime_window::Int = 2000,
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
    N_effective = compute_effective_conductor(N; conductor_mode = conductor_mode)
    isempty(kwargs) || error("unsupported keyword arguments: $(collect(keys(kwargs)))")

    prime_requirement = N_effective
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

    # Prime validity check
    for p in chosen_primes
        (p - 1) % N_effective == 0 || error(
            "prime $p does not satisfy N_effective | p-1 " *
            "(input N=$N; N_effective=$N_effective)")
    end
    length(chosen_primes) >= 2 || error(
        "need at least 2 primes (got $(length(chosen_primes)))")

    # ------- atomic catalog -------
    # Enumerate atomic SL(2, Z/N)-irreps at the user-specified base conductor N.
    # The cyclotomic field is fixed by the requested conductor.
    verbose && println("=== Atomic SL(2, ℤ/$N) irrep catalog ===")
    catalog = build_atomic_catalog(N; max_rank = max_rank, verbose = false)
    verbose && println("  $(length(catalog)) atomic irreps (≤ rank $max_rank)")

    # ------- strata -------
    verbose && println("\n=== Stratum enumeration ===")
    strata_list = strata === nothing ?
        vcat([enumerate_strata(catalog, r) for r in 1:max_rank]...) :
        strata
    verbose && println("  $(length(strata_list)) strata")

    # ------- find MTCs at each prime, for each stratum -------
    verbose && println("\n=== Block-U sweep at primes $chosen_primes ===")

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

    # ------- per-stratum, per-sector modular-data reconstruction -------
    verbose && println("\n=== CRT modular-data reconstruction ===")

    out = ClassifiedMTC[]
    for (st, results_by_prime) in stratum_results
        # Need candidates at ≥ 2 primes to do CRT at all.
        present_primes = sort(collect(keys(results_by_prime)))
        if length(present_primes) < 2
            verbose && println("  stratum $(describe_stratum(st, catalog)): " *
                               "only $(length(present_primes)) prime(s), skip")
            continue
        end

        groups = _groups_for_stratum(results_by_prime, catalog, st, N_effective)

        verbose && println("  stratum (rank $(st.total_dim)): " *
                           "$(length(groups)) modular-data sector(s)")

        # used/fresh split: first half / second half of primes in this group
        for (gi, group) in enumerate(groups)
            group_primes = sort(collect(keys(group)))
            if length(group_primes) < 2
                verbose && println("    sector $gi: only $(length(group_primes)) " *
                                   "prime(s) — modular CRT needs ≥ 2. skip")
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
                                                scale_factor = scale_factor,
                                                reconstruction_bound = reconstruction_bound,
                                                galois_sector = gi,
                                                test_primes = used,
                                                catalog = catalog,
                                                skip_FR = true,
                                                verbose = verbose)
                    if !cmtc.verify_fresh && extra_idx <= length(group_primes)
                        p_new = group_primes[extra_idx]
                        push!(used, p_new)
                        extra_idx += 1
                        verbose && println("    sector $gi: fresh-prime check failed; " *
                                           "retry with additional prime $p_new (used=$used)")
                        continue
                    end
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
            if !cmtc.verify_fresh
                verbose && println("    sector $gi: fresh-prime verification failed; discard")
                continue
            end
            push!(out, cmtc)
            verbose && println("    → $cmtc")
        end
    end

    if skip_FR
        verbose && println("\n=== F/R reconstruction skipped: returning modular-data results only ===")
        verbose && println("\n=== Done: $(length(out)) ClassifiedMTC(s) ===")
        return out
    end

    # ------- Exact (F,R) reconstruction -------
    verbose && println("\n=== Exact (F, R) reconstruction requested ===")
    out = filter(m -> m.verify_fresh, out)
    grouped = _classify_modular_data_by_fusion_rule(out)
    key_to_members = Dict{String, Vector{Tuple{Any, Any}}}()
    discard_indices = Set{Int}()
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
            if occursin("exact F/R reconstruction did not find", msg)
                verbose && println("    exact F/R reconstruction found no solution; leaving members without F/R")
                for i in idxs
                    out[i] = _with_fr_status(out[i], FRNoSolutionFound)
                end
                continue
            elseif _is_reconstruction_unstable_message(msg)
                verbose && println("    exact F/R reconstruction failed; leaving members without F/R")
                for i in idxs
                    out[i] = _with_fr_status(out[i], FRReconstructionFailed)
                end
                continue
            end
            rethrow()
        end
        if fr_result.F === nothing || isempty(fr_result.candidates)
            verbose && println("    exact F/R reconstruction produced no valid candidates; leaving members without F/R")
            for i in idxs
                out[i] = _with_fr_status(out[i], FRNoSolutionFound)
            end
            continue
        end

        # Assign `(F,R)` per `(S,T)` member within the same fusion-rule group.
        for i in idxs
            mtc_i = out[i]
            selection = _select_fr_for_st(fr_result.candidates, mtc_i.Nijk,
                                          mtc_i.S_cyclotomic, mtc_i.T_cyclotomic, mtc_i.N)
            selected = selection.selected
            best_md = _modular_data_roundtrip(selected.F, selected.R, mtc_i.Nijk,
                                              mtc_i.S_cyclotomic,
                                              mtc_i.T_cyclotomic, mtc_i.N)
            best_md = _with_fr_roundtrip_metadata(best_md;
                                                  candidate_index = selection.selected_index,
                                                  galois_exponent = selection.score.galois_exponent)
            verbose && println("    member[$i] branch: cand=$(selection.selected_index), " *
                               "galois=$(selection.score.galois_exponent), " *
                               "perm=$(best_md.best_perm), " *
                               "S_err=$(best_md.S_max), T_err=$(best_md.T_max), ok=$(best_md.ok)")
            if _fr_roundtrip_attachable(best_md)
                out[i] = _with_fr_result(out[i], selected.F, selected.R, best_md)
            else
                verbose && println("    member[$i] selected (F,R) failed exact (S,T) roundtrip; discarding")
                push!(discard_indices, i)
            end
        end
    end

    out = ClassifiedMTC[m for (i, m) in enumerate(out) if !(i in discard_indices)]
    verbose && println("\n=== Done: $(length(out)) ClassifiedMTC(s) ===")
    return out
end
