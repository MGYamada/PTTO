"""
    Phase 5: End-to-end pipeline driver `N → List[ClassifiedMTC]`.

Given a conductor `N`, `classify_mtcs_at_conductor` executes the full
Phase 0 → 1 → 2 → 3 → 4 pipeline and returns a list of fully classified
MTCs, each carrying:

- the SL(2, ℤ/N) stratum (m_λ) decomposition (Phase 0/1),
- F_p-validated modular data (Phase 2),
- Galois-coherent, arithmetic S-matrix in ℤ[√d] (Phase 3),
- `(F, R)` symbols in ℂ satisfying pentagon + hexagon + ribbon (Phase 4),
- a `VerifyReport` summarising residuals.

This module also provides two mid-level helpers:

- `compute_FR_from_ST(Nijk, T_complex; ...)`: given a fusion tensor
  and complex T-eigenvalues, finds a `(F, R)` pair realising that
  modular data via pentagon HC → hexagon HC → ribbon match.

- `classify_from_group(group, all_primes; ...)`: given a Galois-coherent
  group from `group_mtcs_galois_aware`, performs Phase 3 CRT + Phase 4
  `(F, R)` solve and returns a single `ClassifiedMTC`.

Design notes:
- The pipeline is conductor-first: `N` is the outer loop; rank emerges
  from stratum enumeration.
- `(F, R)` is computed over ℂ (`ComplexF64`). Algebraic lift to ℚ(ζ_N)
  is not performed; users can PSLQ/LLL downstream.
- For each Galois-coherent group, Phase 4 may return no match (e.g. if
  the fusion ring has no braided structure at the given T) — such
  groups are reported but `classified_mtc.FR === nothing` in that case.
"""

using LinearAlgebra

# ============================================================
#  ClassifiedMTC: the final pipeline output
# ============================================================

"""
    ClassifiedMTC

The full output of `classify_mtcs_at_conductor` for a single MTC.

Arithmetic / F_p layer (Phase 0–3 output):
- `N`:              conductor
- `rank`:           rank
- `stratum`:        the SL(2, ℤ/N) irrep decomposition `(m_λ)` that gave
                    rise to this MTC
- `Nijk`:           fusion tensor (Galois-invariant integer array)
- `S_Zsqrtd`:       S-matrix as `(a, b)` = `a + b·√d` in ℤ[√d]
                    (after Phase 3 CRT; entries before dividing by
                    `scale · √d`)
- `scale_d`:        the `d` such that `S_Zsqrtd` lives in ℤ[√d]
- `scale_factor`:   the scalar multiplying `S` before reconstruction
                    (so `S_ℂ = S_Zsqrtd / (scale_factor · √d)`)
- `used_primes`:    primes used for CRT reconstruction
- `fresh_primes`:   primes used for cross-validation (may be empty)
- `verify_fresh`:   `true` iff all fresh primes cross-check

Complex-lifted modular data (Phase 4 input):
- `S_complex`:      `Matrix{ComplexF64}`, lifted from `S_Zsqrtd`
- `T_complex`:      `Vector{ComplexF64}` of T-eigenvalues (θ_i), N-th
                    roots of unity on the unit circle

(F, R) symbols (Phase 4 output):
- `F_values`:       pentagon solution in ℂ, `Vector{ComplexF64}`.
                    `nothing` if no match was found.
- `R_values`:       hexagon solution in ℂ, `Vector{ComplexF64}` of
                    length `2 · r_var_count` (forward ⊕ reverse).
                    `nothing` if no match was found.
- `verify_report`:  `VerifyReport` with pentagon, hexagon, ribbon
                    max-residuals. `nothing` if no match was found.

Provenance:
- `galois_sector`:  integer index of the Galois orbit element
                    (1, 2, ... for groups returned by
                    `group_mtcs_galois_aware`)
"""
struct ClassifiedMTC
    N::Int
    rank::Int
    stratum::Stratum
    Nijk::Array{Int, 3}
    S_Zsqrtd::Matrix{Tuple{Int, Int}}
    scale_d::Int
    scale_factor::Int
    used_primes::Vector{Int}
    fresh_primes::Vector{Int}
    verify_fresh::Bool
    S_complex::Matrix{ComplexF64}
    T_complex::Vector{ComplexF64}
    F_values::Union{Vector{ComplexF64}, Nothing}
    R_values::Union{Vector{ComplexF64}, Nothing}
    verify_report::Union{VerifyReport, Nothing}
    galois_sector::Int
end

function Base.show(io::IO, m::ClassifiedMTC)
    FR_status = if m.F_values === nothing
        "(F,R)=none"
    else
        rep = m.verify_report
        pent = rep === nothing ? "?" : string(round(rep.pentagon_max, sigdigits = 2))
        hex = rep === nothing ? "?" : string(round(rep.hexagon_max, sigdigits = 2))
        rib = rep === nothing || rep.ribbon_max === nothing ? "?" :
            string(round(rep.ribbon_max, sigdigits = 2))
        "(F,R) pent=$pent hex=$hex rib=$rib"
    end
    fresh_str = isempty(m.fresh_primes) ? "no fresh" :
        (m.verify_fresh ? "fresh✓" : "fresh✗")
    print(io, "ClassifiedMTC(N=$(m.N), rank=$(m.rank), ",
          "sector=$(m.galois_sector), $(length(m.used_primes)) primes, ",
          "$fresh_str, $FR_status)")
end

# ============================================================
#  compute_FR_from_ST: Phase 4 wrapper
# ============================================================

"""
    compute_FR_from_ST(Nijk, T_complex; ribbon_atol = 1e-8,
                        pentagon_slice = 1, show_progress = false,
                        verbose = false)
        -> NamedTuple{(:F, :R, :report, :n_pentagon, :n_tried, :n_matches,
                       :f_idx, :r_idx)}

Given a fusion tensor `Nijk` and complex T-eigenvalues `T_complex`,
find a pair `(F, R)` of complex F- and R-symbols satisfying pentagon,
hexagon, and the ribbon relation `(R^{ij}_k)² = θ_i θ_j / θ_k` against
the given T.

Algorithm:
1. Set up pentagon system from `Nijk` (via TensorCategories).
2. Solve pentagon via homotopy continuation, optionally with random
   linear slice to break gauge symmetry.
3. For each pentagon solution `F`:
   - Polish with damped Newton.
   - Build hexagon system with F fixed.
   - Solve hexagon via HC.
   - For each R solution, compute ribbon residuals against `T_complex`;
     keep the `(F, R)` pair with the smallest ribbon residual below
     `ribbon_atol`.
4. Return the best match (or all-`nothing` fields if none found).

Caveats:
- Pentagon HC is feasible only for small fusion rings (~10 F-variables
  before mixed-volume blow-up; Ising-sized rings with ~14 F-vars may
  require hours).
- Ribbon match is a NECESSARY condition only: `(R²) = θ_i θ_j / θ_k`
  is blind to the overall sign of R. Two sign-conjugate braidings
  (Fibonacci vs its complex conjugate) will both pass; we return the
  first one HC produces.

Returns a NamedTuple with:
- `F`:            best `Vector{ComplexF64}` or `nothing`
- `R`:            best `Vector{ComplexF64}` or `nothing`
- `report`:       `VerifyReport` or `nothing`
- `n_pentagon`:   number of pentagon solutions found
- `n_tried`:      total `(F, R)` pairs examined
- `n_matches`:    how many of those passed ribbon
- `f_idx`, `r_idx`: indices of the chosen pair (0 if none)
"""
function compute_FR_from_ST(Nijk::Array{Int, 3},
                            T_complex::Vector{ComplexF64};
                            ribbon_atol::Float64 = 1e-8,
                            pentagon_slice::Int = 1,
                            show_progress::Bool = false,
                            verbose::Bool = false)
    r = size(Nijk, 1)
    length(T_complex) == r || error(
        "T_complex length $(length(T_complex)) does not match rank $r")

    # Rank-1 (trivial) MTC: the only fusion category structure is the
    # unit object with F = R = [1]. Pentagon / hexagon / ribbon hold
    # vacuously. `TensorCategories.pentagon_equations` emits no non-
    # trivial polynomials at rank 1, so we short-circuit here rather
    # than let `get_pentagon_system` error.
    if r == 1
        F_trivial = ComplexF64[1.0]
        R_trivial = ComplexF64[1.0]
        # Fields in order: pentagon_max, hexagon_max, ribbon_max,
        # n_pentagon_eqs, n_hexagon_eqs, rank
        report = VerifyReport(0.0, 0.0, 0.0, 0, 0, 1)
        return (F = F_trivial, R = R_trivial, report = report,
                n_pentagon = 1, n_tried = 1, n_matches = 1,
                f_idx = 1, r_idx = 1)
    end

    # Pentagon system
    local _R, eqs, n
    try
        _R, eqs, n = get_pentagon_system(Nijk, r)
    catch err
        # Pointed / trivially-pentagon fusion rings at r ≥ 2 can also
        # emit no non-trivial equations. Treat this as "no pentagon
        # constraint; fall through with trivial F" rather than error.
        verbose && println("  get_pentagon_system threw ($err); " *
                           "treating as trivially pointed")
        return (F = nothing, R = nothing, report = nothing,
                n_pentagon = 0, n_tried = 0, n_matches = 0,
                f_idx = 0, r_idx = 0)
    end
    verbose && println("  Pentagon: $n variables, $(length(eqs)) equations")

    F_sols = solve_pentagon_homotopy(eqs, n;
                                             slice = pentagon_slice,
                                             include_singular = false,
                                             show_progress = show_progress)
    verbose && println("  Pentagon HC: $(length(F_sols)) solutions")

    best = (F = nothing,
            R = nothing,
            report = nothing,
            n_pentagon = length(F_sols),
            n_tried = 0,
            n_matches = 0,
            f_idx = 0,
            r_idx = 0,
            ribbon_max = Inf)

    for (fi, F_raw) in enumerate(F_sols)
        local F
        try
            F = refine_solution_newton(eqs, F_raw; tol = 1e-14)
        catch err
            verbose && println("    F[$fi] Newton refinement failed: $err")
            continue
        end

        local R_sols
        try
            _R_ring, hex_eqs, n_r = get_hexagon_system(Nijk, r, F)
            R_sols = solve_hexagon_homotopy(hex_eqs, n_r;
                                                    show_progress = show_progress)
        catch err
            verbose && println("    F[$fi] hexagon failed: $err")
            continue
        end

        verbose && println("    F[$fi]: $(length(R_sols)) hexagon sols")

        for (ri, R) in enumerate(R_sols)
            best = (; best..., n_tried = best.n_tried + 1)

            local rib_max
            try
                rib = ribbon_residuals(R, T_complex, Nijk)
                rib_max = maximum(rib)
            catch err
                verbose && println("      R[$ri] ribbon failed: $err")
                continue
            end

            if rib_max < ribbon_atol
                best = (; best..., n_matches = best.n_matches + 1)
                if rib_max < best.ribbon_max
                    rep = verify_mtc(F, R, Nijk; T = T_complex)
                    best = (; best...,
                            F = F, R = R, report = rep,
                            f_idx = fi, r_idx = ri,
                            ribbon_max = rib_max)
                end
            end
        end
    end

    # Drop the internal ribbon_max field from the public return
    return (F = best.F, R = best.R, report = best.report,
            n_pentagon = best.n_pentagon, n_tried = best.n_tried,
            n_matches = best.n_matches,
            f_idx = best.f_idx, r_idx = best.r_idx)
end

# ============================================================
#  classify_from_group: CRT + (F, R) solve for one Galois sector
# ============================================================

"""
    classify_from_group(group, N, stratum, all_primes;
                        scale_d, scale_factor = 2,
                        sqrtd_fn = compute_sqrt_d_mod_p,
                        reconstruction_bound = 5,
                        galois_sector = 1,
                        test_primes = nothing,
                        ribbon_atol = 1e-8,
                        skip_FR = false,
                        verbose = false)
        -> ClassifiedMTC

Given a single Galois-coherent group from `group_mtcs_galois_aware`,
perform Phase 3 CRT reconstruction, Phase 4 (F, R) solve, and return
a `ClassifiedMTC`.

Arguments:
- `group::Dict{Int, MTCCandidate}`:  a single group (one sector)
- `N::Int`:                          conductor
- `stratum::Stratum`:                the stratum the group came from
- `all_primes::Vector{Int}`:         primes present in `group`
- `scale_d::Int`:                    `d` for ℤ[√d]
- `scale_factor::Int = 2`:           scalar multiplying S before recon
- `sqrtd_fn`:                        `(d, p) -> Int` returning √d mod p.
                                     For SU(2)_4 / d = 3 /scale_d = 3,
                                     use `compute_sqrt3_cyclotomic_mod_p`
                                     to ensure Galois-consistent choice.
- `reconstruction_bound::Int = 5`:   ℤ[√d] coefficient bound for rational
                                     reconstruction.
- `galois_sector::Int = 1`:          sector index for provenance.
- `test_primes::Vector{Int}=nothing`: which primes to use for CRT.
                                     If `nothing`, uses first
                                     `max(2, length(all_primes) ÷ 2)`
                                     primes; remaining are fresh.
- `ribbon_atol::Float64 = 1e-8`:     tolerance for the Phase 4 ribbon
                                     match.
- `skip_FR::Bool = false`:           if true, only do Phase 0–3 and leave
                                     (F, R) as `nothing`. Useful for
                                     rank/complexity beyond the pentagon
                                     HC limit (~10 F-vars).
- `verbose::Bool = false`:           print progress.
"""
function classify_from_group(group::Dict{Int, MTCCandidate},
                             N::Int,
                             stratum::Stratum,
                             all_primes::Vector{Int};
                             scale_d::Int,
                             scale_factor::Int = 2,
                             sqrtd_fn = compute_sqrt_d_mod_p,
                             reconstruction_bound::Int = 5,
                             galois_sector::Int = 1,
                             test_primes::Union{Vector{Int}, Nothing} = nothing,
                             ribbon_atol::Float64 = 1e-8,
                             skip_FR::Bool = false,
                             verbose::Bool = false)
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

    # -------- Phase 3: CRT in ℤ[√d] --------
    used_subgroup = Dict(p => group[p] for p in used)
    recon_S = reconstruct_S_matrix(used_subgroup;
                                    scale_d = scale_d,
                                    bound = reconstruction_bound,
                                    sqrtd_fn = sqrtd_fn)

    # Verify at fresh primes (if any)
    verify_fresh = isempty(fresh)  # vacuously true if none
    if !isempty(fresh)
        all_ok = true
        for p in fresh
            ok = verify_reconstruction(recon_S, group[p], scale_d;
                                        scale = scale_factor,
                                        sqrtd_fn = sqrtd_fn)
            all_ok &= ok
            verbose && println("    fresh p=$p: $(ok ? "✓" : "✗")")
        end
        verify_fresh = all_ok
    end

    # -------- Phase 4 input: lift (S, T) to ℂ --------
    rep = first(values(group))
    zeta_Fp = find_zeta_in_Fp(N, rep.p)
    S_ℂ, T_ℂ, Nijk = lift_mtc_candidate(rep, recon_S;
                                                d = scale_d,
                                                N = N,
                                                zeta_Fp = zeta_Fp,
                                                scale = scale_factor)

    # -------- Phase 4: (F, R) solve + verify --------
    rank = size(Nijk, 1)

    if skip_FR
        return ClassifiedMTC(N, rank, stratum, Nijk,
                             recon_S, scale_d, scale_factor,
                             used, fresh, verify_fresh,
                             S_ℂ, T_ℂ, nothing, nothing, nothing,
                             galois_sector)
    end
    verbose && println("  running pentagon/hexagon on rank=$rank...")

    local fr_result
    try
        fr_result = compute_FR_from_ST(Nijk, T_ℂ;
                                        ribbon_atol = ribbon_atol,
                                        verbose = verbose)
    catch err
        verbose && println("  Phase 4 failed: $err")
        return ClassifiedMTC(N, rank, stratum, Nijk, recon_S,
                             scale_d, scale_factor,
                             used, fresh, verify_fresh,
                             S_ℂ, T_ℂ, nothing, nothing, nothing,
                             galois_sector)
    end

    return ClassifiedMTC(N, rank, stratum, Nijk, recon_S,
                         scale_d, scale_factor,
                         used, fresh, verify_fresh,
                         S_ℂ, T_ℂ,
                         fr_result.F, fr_result.R, fr_result.report,
                         galois_sector)
end

# ============================================================
#  classify_mtcs_at_conductor: top-level driver
# ============================================================

"""
    classify_mtcs_at_conductor(N; max_rank, primes, strata = nothing,
                                scale_d = 3, scale_factor = 2,
                                sqrtd_fn = nothing,
                                verlinde_threshold = 3,
                                max_block_dim = 3,
                                reconstruction_bound = 5,
                                ribbon_atol = 1e-8,
                                skip_FR = false,
                                verbose = true)
        -> Vector{ClassifiedMTC}

Fully automated MTC classification for a given conductor `N`.

Pipeline:
  Phase 0: build SL(2, ℤ/N) atomic irrep catalog (≤ max_rank)
  Phase 1: enumerate strata (m_λ) with Σ m_λ d_λ = r for each r ≤ max_rank
  Phase 2: for each stratum, sweep block-U at each prime → MTC candidates
  Phase 3: Galois-aware grouping + CRT reconstruction in ℤ[√d]
  Phase 4: lift to ℂ, solve pentagon/hexagon, verify ribbon against T

Returns one `ClassifiedMTC` per Galois sector per stratum that yields a
valid MTC. If `skip_FR = true`, Phase 4 is skipped and `(F, R)` is left
as `nothing` (useful for large ranks where pentagon HC is infeasible).

Note on conductor: `N` here is the SL(2, ℤ/N)-rep conductor —
equivalently, the order of `T` as a matrix. An MTC's S-matrix may live
in a larger cyclotomic field than ℚ(ζ_N); in NRWW's convention the
full MTC conductor is `max(cond(S), cond(T))`. Fibonacci, for example,
has `cond(T) = 5` but its S involves `D = √(2+φ)` which is NOT in ℚ(ζ_5).
Calling this function with the T-only conductor will miss such MTCs —
they only appear when `N` is large enough to accommodate `S` too.
SU(2)_4 (conductor 24) is reachable; Fibonacci is typically reached
at a larger conductor.

Arguments:
- `N::Int`:                        conductor
- `max_rank::Int`:                 maximum rank to consider
- `primes::Vector{Int}`:           good primes (must satisfy `N | p-1`).
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
- `scale_d::Int = 3`:              ℤ[√d] the S-matrix is expected to
                                   live in. Must be chosen to match the
                                   MTC's quadratic field: e.g. 3 for
                                   SU(2)_4 / Ising, 5 for MTCs involving
                                   the golden ratio, 2 for √2-pointed
                                   MTCs. A wrong `scale_d` will silently
                                   return 0 classified MTCs.
- `scale_factor::Int = 2`:         scalar multiplying `S` at
                                   reconstruction (matches
                                   `reconstruct_S_matrix` convention).
- `sqrtd_fn`:                      custom √d-in-F_p function. If
                                   `nothing` (default), chooses the
                                   cyclotomic variant:
                                   `compute_sqrt3_cyclotomic_mod_p`
                                   (needs 24 | p-1) for `scale_d = 3`,
                                   `compute_sqrt2_cyclotomic_mod_p`
                                   (needs 24 | p-1) for `scale_d = 2`,
                                   `compute_sqrt5_cyclotomic_mod_p`
                                   (needs 5 | p-1) for `scale_d = 5`,
                                   else the non-cyclotomic
                                   `compute_sqrt_d_mod_p` (Tonelli-
                                   Shanks, NOT Galois-consistent).
- `verlinde_threshold::Int = 3`:   max absolute fusion coefficient
                                   allowed (beyond which the candidate
                                   is rejected). See
                                   `find_mtcs_at_prime`.
- `max_block_dim::Int = 3`:        cap on the degenerate T-eigenspace
                                   dimension for naive O(n) Cayley
                                   sweep. Raise with caution (O(4)(F_p)
                                   has ~10¹¹ points).
- `reconstruction_bound::Int = 5`: coefficient bound for ℤ[√d]
                                   rational reconstruction.
- `ribbon_atol::Float64 = 1e-8`:   Phase 4 ribbon-match tolerance.
- `skip_FR::Bool = false`:         skip Phase 4. Useful if
                                   `max_rank ≥ 5` and pentagon HC would
                                   blow up.
- `verbose::Bool = true`:          print per-phase progress.
"""
function classify_mtcs_at_conductor(N::Int;
                                    max_rank::Int,
                                    primes::Vector{Int},
                                    strata::Union{Nothing, Vector{Stratum}} = nothing,
                                    scale_d::Int = 3,
                                    scale_factor::Int = 2,
                                    sqrtd_fn = nothing,
                                    verlinde_threshold::Int = 3,
                                    max_block_dim::Int = 3,
                                    reconstruction_bound::Int = 5,
                                    ribbon_atol::Float64 = 1e-8,
                                    skip_FR::Bool = false,
                                    verbose::Bool = true)
    # Prime validity check
    for p in primes
        (p - 1) % N == 0 || error(
            "prime $p does not satisfy N | p-1 for N=$N")
    end
    length(primes) >= 2 || error(
        "need at least 2 primes (got $(length(primes)))")

    # Default sqrtd_fn for the common cases. Cyclotomic variants are
    # Galois-consistent across primes; Tonelli-Shanks is not.
    if sqrtd_fn === nothing
        sqrtd_fn = if scale_d == 3
            (d, p) -> compute_sqrt3_cyclotomic_mod_p(p)
        elseif scale_d == 2
            (d, p) -> compute_sqrt2_cyclotomic_mod_p(p)
        elseif scale_d == 5
            (d, p) -> compute_sqrt5_cyclotomic_mod_p(p)
        else
            # Generic Tonelli-Shanks fallback. NOT Galois-consistent —
            # CRT reconstruction may pick the wrong sector.
            verbose && @warn "Using generic (Tonelli-Shanks) √d for scale_d=$scale_d; " *
                             "this is not Galois-consistent. Consider supplying an " *
                             "explicit cyclotomic `sqrtd_fn` if the result looks wrong."
            compute_sqrt_d_mod_p
        end
    end

    # ------- Phase 0: atomic catalog -------
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
    verbose && println("\n=== Phase 2: block-U sweep at primes $primes ===")

    # Collect (stratum, Dict{p => candidates}) only for strata that yield
    # something at at least one prime.
    stratum_results = Vector{Tuple{Stratum, Dict{Int, Vector{MTCCandidate}}}}()
    skipped_reasons = Dict{String, Int}()
    for (si, st) in enumerate(strata_list)
        results_by_prime = Dict{Int, Vector{MTCCandidate}}()
        any_found = false
        any_prime_errored = false
        last_err = ""
        for p in primes
            local cands
            try
                cands = find_mtcs_at_prime(catalog, st, p;
                                            verlinde_threshold = verlinde_threshold,
                                            max_block_dim = max_block_dim)
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

    # ------- Phase 3 + 4: per-stratum, per-sector -------
    verbose && println("\n=== Phase 3 + 4: CRT + (F, R) classification ===")

    out = ClassifiedMTC[]
    for (st, results_by_prime) in stratum_results
        # Need candidates at ≥ 2 primes to do CRT at all.
        present_primes = sort(collect(keys(results_by_prime)))
        if length(present_primes) < 2
            verbose && println("  stratum $(describe_stratum(st, catalog)): " *
                               "only $(length(present_primes)) prime(s), skip")
            continue
        end

        # Galois-aware grouping
        anchor = present_primes[1]
        local groups
        try
            groups = group_mtcs_galois_aware(results_by_prime, anchor;
                                              scale_d = scale_d,
                                              sqrtd_fn = sqrtd_fn)
        catch err
            verbose && println("  stratum: galois grouping failed: $err")
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
            used = group_primes[1:half]
            cur_ribbon_atol = ribbon_atol

            local cmtc
            try
                cmtc = classify_from_group(group, N, st, group_primes;
                                            scale_d = scale_d,
                                            scale_factor = scale_factor,
                                            sqrtd_fn = sqrtd_fn,
                                            reconstruction_bound = reconstruction_bound,
                                            galois_sector = gi,
                                            test_primes = used,
                                            ribbon_atol = cur_ribbon_atol,
                                            skip_FR = skip_FR,
                                            verbose = verbose)
            catch err
                verbose && println("    sector $gi: classify_from_group failed: $err")
                continue
            end
            push!(out, cmtc)
            verbose && println("    → $cmtc")
        end
    end

    verbose && println("\n=== Done: $(length(out)) ClassifiedMTC(s) ===")
    return out
end
