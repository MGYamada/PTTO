"""
    End-to-end pipeline driver `N → List[ClassifiedMTC]`.

Given a conductor `N`, `classify_mtcs_at_conductor` executes the full pipeline and returns a list of fully classified
MTCs, each carrying:

- the SL(2, ℤ/N) stratum (m_λ) decomposition (Phase 0/1),
- F_p-validated modular data (Phase 2),
- Galois-coherent, arithmetic S-matrix in ℤ[√d] (Phase 3),
- exact `(S, T)` lifted to `Q(ζ_N)` with `N` fixed as the base conductor.

This module also provides two mid-level helpers:

- `compute_FR_from_ST(Nijk; ...)`: exact cyclotomic Phase 4 helper.

- `classify_from_group(group, all_primes; ...)`: given a Galois-coherent
  group from `group_mtcs_galois_aware`, performs Phase 3 CRT + Phase 4
  exact lift and returns a single `ClassifiedMTC`.

And one high-level convenience API:

- `classify_mtcs_auto(N; ...)`: automatically chooses
  `conductor_mode`, `primes`, and `max_rank` from candidate
  lists, then runs `classify_mtcs_at_conductor`.

Design notes:
- The pipeline is conductor-first: `N` is the outer loop; rank emerges
  from stratum enumeration.
- Modular data is lifted exactly to `Q(ζ_N)`.  The effective conductor is
  fixed to the user-supplied `N`.
- `skip_FR = false` runs exact cyclotomic `(F, R)` reconstruction.
"""

using LinearAlgebra

# ============================================================
#  ClassifiedMTC: the final pipeline output
# ============================================================

"""
    ClassifiedMTC

The full output of `classify_mtcs_at_conductor` for a single MTC.

Arithmetic / F_p layer (Phase 0–3 output):
- `N`:              conductor used internally by the pipeline
- `N_input`:        user-requested conductor (for provenance)
- `rank`:           rank
- `stratum`:        the SL(2, ℤ/N) irrep decomposition `(m_λ)` that gave
                    rise to this MTC
- `Nijk`:           fusion tensor (Galois-invariant integer array)
- `S_Zsqrtd`:       S-matrix as `(a, b)` = `a + b·√d` in ℤ[√d]
                    (after Phase 3 CRT; entries before dividing by
                    `scale · √d`)
- `quadratic_d`:    the inferred `d` such that `S_Zsqrtd` lives in ℤ[√d]
- `scale_factor`:   the scalar multiplying `S` before reconstruction
                    (so `S_cyclotomic = S_Zsqrtd / (scale_factor · √d)`)
- `used_primes`:    primes used for CRT reconstruction
- `fresh_primes`:   primes used for cross-validation (may be empty)
- `verify_fresh`:   `true` iff all fresh primes cross-check

Exact modular data (Phase 4 input):
- `S_cyclotomic`:   S-matrix over `Q(ζ_N)`
- `T_cyclotomic`:   T-eigenvalues over `Q(ζ_N)`

Exact `(F, R)` layer:
- `F_values`:       reserved for exact pentagon data; currently `nothing`.
- `R_values`:       reserved for exact braiding data; currently `nothing`.
- `verify_report`:  reserved for exact verification data; currently `nothing`.

Provenance:
- `galois_sector`:  integer index of the Galois orbit element
                    (1, 2, ... for groups returned by
                    `group_mtcs_galois_aware`)
"""
struct ClassifiedMTC
    N::Int
    N_input::Int
    rank::Int
    stratum::Stratum
    Nijk::Array{Int, 3}
    S_Zsqrtd::Matrix{Tuple{Int, Int}}
    quadratic_d::Int
    scale_factor::Int
    used_primes::Vector{Int}
    fresh_primes::Vector{Int}
    verify_fresh::Bool
    S_cyclotomic::Any
    T_cyclotomic::Any
    F_values::Union{Vector, Nothing}
    R_values::Union{Vector, Nothing}
    verify_report::Any
    galois_sector::Int
end

function _permute_fusion_tensor(Nijk::Array{Int, 3}, perm::Vector{Int})
    return Nijk[perm, perm, perm]
end

function _all_permutations(v::Vector{Int})
    if length(v) <= 1
        return [copy(v)]
    end
    out = Vector{Vector{Int}}()
    for i in eachindex(v)
        head = v[i]
        tail = Vector{Int}(undef, length(v) - 1)
        ti = 1
        for j in eachindex(v)
            if j != i
                tail[ti] = v[j]
                ti += 1
            end
        end
        for rest in _all_permutations(tail)
            push!(out, vcat(head, rest))
        end
    end
    return out
end

"""
    fusion_rule_key(Nijk) -> String

Backward-compatible alias for `canonical_rule(Nijk)`.
"""
function fusion_rule_key(Nijk::AbstractArray{<:Integer, 3})
    return canonical_rule(Nijk)
end

function _classify_modular_data_by_fusion_rule(classified::Vector{ClassifiedMTC})
    grouped = Dict{String, Vector{Int}}()
    for (i, cmtc) in enumerate(classified)
        key = fusion_rule_key(cmtc.Nijk)
        push!(get!(grouped, key, Int[]), i)
    end
    return grouped
end

function _with_fr_result(c::ClassifiedMTC,
                         F::Union{Vector, Nothing},
                         R::Union{Vector, Nothing},
                         report;
                         S = c.S_cyclotomic,
                         T = c.T_cyclotomic)
    return ClassifiedMTC(c.N, c.N_input, c.rank, c.stratum, c.Nijk,
                         c.S_Zsqrtd, c.quadratic_d, c.scale_factor,
                         c.used_primes, c.fresh_primes, c.verify_fresh,
                         S, T,
                         F, R, report, c.galois_sector)
end

# ============================================================
#  classify_mtcs_auto: user-friendly auto-parameter wrapper
# ============================================================

function cyclotomic_requirement(d::Int)
    d in (2, 3) && return 24
    d == 5 && return 5
    return 1
end

function _prime_requirement_for_reconstruction(N_eff::Int, d::Int, sqrtd_fn)
    req = N_eff
    if sqrtd_fn === nothing
        req = lcm(req, cyclotomic_requirement(d))
    end
    return req
end

function compute_effective_conductor(N::Int, args...;
                                     conductor_mode::Symbol = :full_mtc)
    conductor_mode == :full_mtc ||
        error("conductor_mode=:$(conductor_mode) was removed in v0.5.0. Use :full_mtc.")
    return N
end

function _quadratic_candidates_for_conductor(N::Int)
    ds = Int[1]
    N % 8 == 0 && push!(ds, 2)
    N % 12 == 0 && push!(ds, 3)
    N % 5 == 0 && push!(ds, 5)
    return unique(ds)
end

function _default_quadratic_d_for_conductor(N::Int)
    N % 5 == 0 && return 5
    N % 12 == 0 && return 3
    return 1
end

"""
    classify_mtcs_auto(N::Int;
                       max_rank_candidates = [2, 3, 4, 5],
                       d_candidates = [1, 2, 3, 5, 6, 7, 10],
                       conductor_modes = [:full_mtc],
                       min_primes = 4,
                       prime_start = 29,
                       prime_max = 2000,
                       N_eff_max = typemax(Int),
                       stagnation_k = 3,
                       max_attempts = typemax(Int),
                       strata = nothing,
                       scale_factor = 2,
                       sqrtd_fn = nothing,
                       verlinde_threshold = 3,
                       max_block_dim = 3,
                       search_mode = :groebner,
                       max_units_for_groebner = typemax(Int),
                       groebner_allow_fallback = true,
                       precheck_unit_axiom = true,
                       reconstruction_bound = 50,
                       skip_FR = false,
                       verbose = true)
        -> NamedTuple

Auto-select wrapper around `classify_mtcs_at_conductor`.

This function is intended as the recommended public entry point for
users who do not want to manually specify `max_rank`, `primes`, and
quadratic reconstruction branches. The effective conductor is fixed to
the user-supplied `N`.

For each previously unseen effective conductor, the driver tries
`(conductor_mode, max_rank)` combinations and records stage metadata.
Since `N_effective = N`, duplicate stages are skipped. Search stops
when any of:

- no new MTCs for `stagnation_k` consecutive executed stages,
- `N_eff_candidate > N_eff_max`,
- total attempted runs reaches `max_attempts`.

Returns:

1. `classified`: deduplicated union of MTCs found across all stages.
2. reproducibility metadata:
   - `N_input`
   - `N_effective`
   - `d`
   - `conductor_mode`
   - `primes`
   - `max_rank`
   - `attempts`
   - `history`

Prime selection rule for each attempted `(N_effective, d)` stage:
- Let `req = lcm(N_effective, cyclotomic_requirement(d))` when using
  the built-in √d selector, otherwise `req = N_effective`.
- Choose the first `min_primes` primes `p > prime_start` with
  `(p - 1) % req == 0`.
"""
function classify_mtcs_auto(N::Int;
                            max_rank_candidates::Vector{Int} = [2, 3, 4, 5],
                            d_candidates::Vector{Int} = [1, 2, 3, 5, 6, 7, 10],
                            conductor_modes::Vector{Symbol} = [:full_mtc],
                            min_primes::Int = 4,
                            prime_start::Int = 29,
                            prime_max::Int = 2000,
                            N_eff_max::Int = typemax(Int),
                            stagnation_k::Int = 3,
                            max_attempts::Int = typemax(Int),
                            strata::Union{Nothing, Vector{Stratum}} = nothing,
                            scale_factor::Int = 2,
                            sqrtd_fn = nothing,
                            verlinde_threshold::Int = 3,
                            max_block_dim::Int = 3,
                            search_mode::Symbol = :groebner,
                            max_units_for_groebner::Int = typemax(Int),
                            groebner_allow_fallback::Bool = true,
                            precheck_unit_axiom::Bool = true,
                            reconstruction_bound::Int = 50,
                            skip_FR::Bool = false,
                            verbose::Bool = true)
    N >= 1 || error("N must be positive, got $N")
    min_primes >= 2 || error("min_primes must be ≥ 2, got $min_primes")
    !isempty(max_rank_candidates) || error("max_rank_candidates must be non-empty")
    !isempty(d_candidates) || error("d_candidates must be non-empty")
    !isempty(conductor_modes) || error("conductor_modes must be non-empty")
    stagnation_k >= 1 || error("stagnation_k must be ≥ 1, got $stagnation_k")
    max_attempts >= 1 || error("max_attempts must be ≥ 1, got $max_attempts")

    last_result = ClassifiedMTC[]
    last_meta = (N_input = N, N_effective = N, d = 1,
                 conductor_mode = :full_mtc, primes = Int[],
                 max_rank = 0, attempts = 0)
    attempts = 0
    n_stagnant = 0

    seen_stages = Set{Tuple{Int, Int}}()
    seen_signatures = Set{String}()
    history = NamedTuple[]

    mtc_signature(m::ClassifiedMTC) = begin
        nijk_key = join(vec(m.Nijk), ",")
        t_key = join(string.(m.T_cyclotomic), ",")
        string(m.rank, "|", m.galois_sector, "|", nijk_key, "|", t_key)
    end

    for d in d_candidates
        N_eff_candidate = compute_effective_conductor(N, d)
        if N_eff_candidate > N_eff_max
            push!(history, (d = d, N_effective = N_eff_candidate, executed = false,
                            success = false, reason = "N_eff_max_exceeded",
                            attempts = 0, new_mtcs = 0))
            break
        end
        stage_key = (N_eff_candidate, d)
        if stage_key in seen_stages
            push!(history, (d = d, N_effective = N_eff_candidate, executed = false,
                            success = false, reason = "duplicate_stage",
                            attempts = 0, new_mtcs = 0))
            continue
        end
        push!(seen_stages, stage_key)

        stage_attempts = 0
        stage_success = false
        stage_reason = "no_attempts"
        stage_new_mtcs = 0
        stage_result = ClassifiedMTC[]

        for conductor_mode in conductor_modes
            attempts >= max_attempts && break
            conductor_mode == :full_mtc || error(
                "conductor_mode=:$(conductor_mode) was removed in v0.5.0. Use :full_mtc.")

            req = _prime_requirement_for_reconstruction(N_eff_candidate, d, sqrtd_fn)

            local chosen_primes
            try
                chosen_primes = select_admissible_primes(req;
                                                         min_count = min_primes,
                                                         window = prime_max - prime_start,
                                                         start_from = prime_start)
            catch err
                stage_reason = sprint(showerror, err)
                continue
            end

            for max_rank in max_rank_candidates
                    attempts >= max_attempts && break
                    attempts += 1
                    stage_attempts += 1
                    verbose && println("AUTO attempt #$attempts: " *
                                       "d=$d N_eff=$N_eff_candidate req=$req " *
                                       "mode=$conductor_mode " *
                                       "max_rank=$max_rank primes=$chosen_primes")

                    last_meta = (N_input = N, N_effective = N_eff_candidate, d = d,
                                 conductor_mode = conductor_mode,
                                 primes = copy(chosen_primes), max_rank = max_rank,
                                 attempts = attempts)

                    local classified
                    try
                        classified = classify_mtcs_at_conductor(N_eff_candidate;
                                                                max_rank = max_rank,
                                                                primes = chosen_primes,
                                                                strata = strata,
                                                                quadratic_d = d,
                                                                scale_factor = scale_factor,
                                                                conductor_mode = conductor_mode,
                                                                sqrtd_fn = sqrtd_fn,
                                                                verlinde_threshold = verlinde_threshold,
                                                                max_block_dim = max_block_dim,
                                                                search_mode = search_mode,
                                                                max_units_for_groebner = max_units_for_groebner,
                                                                groebner_allow_fallback = groebner_allow_fallback,
                                                                precheck_unit_axiom = precheck_unit_axiom,
                                                                reconstruction_bound = reconstruction_bound,
                                                                skip_FR = skip_FR,
                                                                verbose = verbose)
                    catch err
                        stage_reason = sprint(showerror, err)
                        verbose && println("  AUTO attempt #$attempts failed: $stage_reason")
                        continue
                    end

                    if !isempty(classified)
                        stage_success = true
                        stage_reason = "ok"
                        stage_result = classified
                        break
                    end
                    stage_reason = "no_classified"
            end
            stage_success && break
        end

        if attempts >= max_attempts && !stage_success
            stage_reason = "budget_reached"
        end

        if stage_success
            for m in stage_result
                sig = mtc_signature(m)
                if !(sig in seen_signatures)
                    push!(seen_signatures, sig)
                    push!(last_result, m)
                    stage_new_mtcs += 1
                end
            end
        end

        n_stagnant = stage_new_mtcs == 0 ? n_stagnant + 1 : 0
        push!(history, (d = d, N_effective = N_eff_candidate, executed = true,
                        success = stage_success, reason = stage_reason,
                        attempts = stage_attempts, new_mtcs = stage_new_mtcs))

        if n_stagnant >= stagnation_k || attempts >= max_attempts
            break
        end
    end

    return (classified = last_result,
            N_input = last_meta.N_input,
            N_effective = last_meta.N_effective,
            d = last_meta.d,
            conductor_mode = last_meta.conductor_mode,
            primes = last_meta.primes,
            max_rank = last_meta.max_rank,
            attempts = last_meta.attempts,
            history = history)
end


_is_reconstruction_unstable_message(msg::AbstractString) = begin
    low = lowercase(msg)
    return occursin("reconstruct", low) ||
           occursin("inconsistent", low) ||
           occursin("fresh-prime", low) ||
           occursin("crt", low)
end

function _phase4_removed_error()
    error("exact Phase 4 did not find a cyclotomic F/R solution for this input")
end

function _select_fr_for_st(candidates, Nijk, S_cyc, T_cyc, N)
    isempty(candidates) && _phase4_removed_error()
    all_scores = NamedTuple[]
    best_idx = 0
    best_score = nothing
    best_candidate = nothing
    ctx = CyclotomicContext(N)
    galois_units = [a for a in 1:N if gcd(a, N) == 1]
    for (ci, cand) in enumerate(candidates)
        for a in galois_units
            F = a == 1 ? cand.F : galois_action(ctx, cand.F, a)
            R = a == 1 ? cand.R : galois_action(ctx, cand.R, a)
            score = _score_fr_st_match(F, R, Nijk, S_cyc, T_cyc, N;
                                       candidate_index = ci)
            score = merge(score, (galois_exponent = a,
                                  order_key = (score.ok ? 0 : 1,
                                               string(score.T_max),
                                               string(score.S_max),
                                               ci,
                                               a)))
            push!(all_scores, score)
            if best_score === nothing || score.order_key < best_score.order_key
                best_score = score
                best_idx = ci
                best_candidate = (F = F, R = R, report = cand.report)
            end
        end
    end
    return (selected = best_candidate,
            score = best_score,
            selected_index = best_idx,
            selected_ok = best_score.ok,
            all_scores = all_scores)
end

function _normalize_T_by_unit(T_cyc::AbstractVector)
    isempty(T_cyc) && return collect(T_cyc)
    t0 = T_cyc[1]
    iszero(t0) && return collect(T_cyc)
    return [t / t0 for t in T_cyc]
end


function _matrix_dim(M)
    return M isa MatElem ? (nrows(M), ncols(M)) : size(M)
end

_matrix_entry(M, i::Int, j::Int) = M[i, j]

function _forward_r_var_count(Nijk::Array{Int,3})
    r = size(Nijk, 1)
    return sum(Nijk[i, j, k]^2 for i in 1:r, j in 1:r, k in 1:r)
end

function _forward_R_values(R_values::Vector, Nijk::Array{Int,3})
    n = _forward_r_var_count(Nijk)
    length(R_values) == n && return R_values
    length(R_values) == 2n && return R_values[1:n]
    error("R_values has length $(length(R_values)); expected $n or $(2n)")
end

function _extract_R_block_exact(R_values::Vector, Nijk::Array{Int,3},
                                i::Int, j::Int, k::Int)
    n = Nijk[i, j, k]
    n == 0 && return nothing
    positions, _ = _braiding_block_positions(Nijk)
    pos = positions[(i, j, k)]
    K = parent(R_values[1])
    M = zero_matrix(K, n, n)
    for a in 1:n, b in 1:n
        M[a, b] = R_values[pos[(a - 1) * n + b]]
    end
    return M
end

function _trace_exact(M)
    n = nrows(M)
    t = zero(base_ring(M))
    for i in 1:n
        t += M[i, i]
    end
    return t
end

function _fusion_automorphisms_fixing_unit(Nijk::Array{Int,3})
    r = size(Nijk, 1)
    r == 1 && return [Int[1]]
    autos = Vector{Vector{Int}}()
    for rest in _all_permutations(collect(2:r))
        perm = vcat(1, rest)
        ok = true
        @inbounds for i in 1:r, j in 1:r, k in 1:r
            if Nijk[perm[i], perm[j], perm[k]] != Nijk[i, j, k]
                ok = false
                break
            end
        end
        ok && push!(autos, perm)
    end
    isempty(autos) && push!(autos, collect(1:r))
    return autos
end

function _modular_data_from_FR(R_values::Vector,
                               Nijk::Array{Int,3},
                               S_target)
    r = size(Nijk, 1)
    dims = _matrix_dim(S_target)
    dims == (r, r) || error("S_target has shape $dims; expected ($r, $r)")

    K = parent(_matrix_entry(S_target, 1, 1))
    R_fwd = _forward_R_values(R_values, Nijk)
    d = [_matrix_entry(S_target, a, 1) / _matrix_entry(S_target, 1, 1) for a in 1:r]
    D = inv(_matrix_entry(S_target, 1, 1))

    S_from_R = zero_matrix(K, r, r)
    for a in 1:r, b in 1:r
        acc = zero(K)
        for c in 1:r
            Nijk[a, b, c] == 0 && continue
            Rab = _extract_R_block_exact(R_fwd, Nijk, a, b, c)
            Rba = _extract_R_block_exact(R_fwd, Nijk, b, a, c)
            acc += d[c] * _trace_exact(Rba * Rab)
        end
        S_from_R[a, b] = acc / D
    end

    T_from_R = Vector{typeof(K(1))}(undef, r)
    for a in 1:r
        acc = zero(K)
        for c in 1:r
            Nijk[a, a, c] == 0 && continue
            Raa = _extract_R_block_exact(R_fwd, Nijk, a, a, c)
            acc += d[c] * _trace_exact(Raa)
        end
        T_from_R[a] = acc / d[a]
    end
    if !iszero(T_from_R[1])
        t0 = T_from_R[1]
        T_from_R = [t / t0 for t in T_from_R]
    end
    return (S = S_from_R, T = T_from_R)
end

function _modular_data_roundtrip(F_values::Vector,
                                 R_values::Vector,
                                 Nijk::Array{Int,3},
                                 S_target,
                                 T_target,
                                 N::Int)
    r = size(Nijk, 1)
    length(T_target) == r || error("T_target has length $(length(T_target)); expected $r")
    md_from_fr = _modular_data_from_FR(R_values, Nijk, S_target)
    S_from_R = md_from_fr.S
    T_from_R = md_from_fr.T
    K = parent(_matrix_entry(S_target, 1, 1))

    automorphisms = _fusion_automorphisms_fixing_unit(Nijk)
    best = nothing
    for perm in automorphisms
        S_diffs = [S_from_R[perm[i], perm[j]] - _matrix_entry(S_target, i, j)
                   for i in 1:r, j in 1:r]
        T_diffs = [T_from_R[perm[i]] - T_target[i] for i in 1:r]
        S_ok = all(iszero, S_diffs)
        T_ok = all(iszero, T_diffs)
        S_err = S_ok ? zero(K) : first(x for x in S_diffs if !iszero(x))
        T_err = T_ok ? zero(K) : first(x for x in T_diffs if !iszero(x))
        score = (ok = S_ok && T_ok,
                 S_max = S_err,
                 T_max = T_err,
                 best_perm = perm,
                 S_roundtrip = S_from_R,
                 T_roundtrip = T_from_R)
        if best === nothing || (score.ok && !best.ok)
            best = score
        end
    end
    return best
end

function _score_fr_st_match(F_values::Vector,
                            R_values::Vector,
                            Nijk::Array{Int,3},
                            S_target,
                            T_target,
                            N::Int;
                            candidate_index::Int)
    md = _modular_data_roundtrip(F_values, R_values, Nijk, S_target, T_target, N)
    return merge(md, (candidate_index = candidate_index,
                      order_key = (md.ok ? 0 : 1,
                                   string(md.T_max),
                                   string(md.S_max),
                                   candidate_index)))
end

function _branch_consistency_precheck(results_by_prime::Dict{Int, Vector{MTCCandidate}},
                                      anchor_prime::Int,
                                      quadratic_d::Int,
                                      sqrtd_fn;
                                      reconstruction_bound::Int = 50,
                                      branch_sign_getter = nothing,
                                      branch_sign_setter = nothing,
                                      verbose::Bool = false)
    haskey(results_by_prime, anchor_prime) || return Int[]
    anchor_cands = results_by_prime[anchor_prime]
    isempty(anchor_cands) && return Int[]
    contradictory = Int[]

    function compatibility_status(cands::Vector{MTCCandidate}, p::Int, trial_signs, orig_sign)
        saw_comparable = false
        for anchor_c in anchor_cands
            nrow = size(anchor_c.S_Fp, 1)
            ncol = size(anchor_c.S_Fp, 2)
            for c in cands
                # NOTE:
                # unit_index can differ across primes for the same fusion
                # ring candidate (ordering conventions / finite-field
                # normalization). Requiring exact unit_index equality here
                # over-filters compatible pairs and causes false
                # "branch contradiction" reports.
                c.N == anchor_c.N || continue
                size(c.S_Fp, 1) == nrow || continue
                size(c.S_Fp, 2) == ncol || continue
                saw_comparable = true
                local any_sign_ok = false
                for sgn in trial_signs
                    if branch_sign_setter !== nothing && sgn !== nothing
                        branch_sign_setter(p, sgn)
                    end
                    try
                        s_anchor = sqrtd_fn(quadratic_d, anchor_prime)
                        s_p = sqrtd_fn(quadratic_d, p)
                        two_s_anchor = mod(2 * s_anchor, anchor_prime)
                        two_s_p = mod(2 * s_p, p)
                        matrix_by_prime = Dict(
                            anchor_prime => [mod(two_s_anchor * anchor_c.S_Fp[i, j], anchor_prime)
                                             for i in 1:nrow, j in 1:ncol],
                            p => [mod(two_s_p * c.S_Fp[i, j], p)
                                  for i in 1:nrow, j in 1:ncol])
                        reconstruct_matrix_in_Z_sqrt_d(matrix_by_prime, quadratic_d;
                                                       bound = reconstruction_bound, sqrtd_fn = sqrtd_fn)
                        any_sign_ok = true
                    catch
                        continue
                    finally
                        if branch_sign_setter !== nothing && orig_sign !== nothing
                            branch_sign_setter(p, orig_sign)
                        end
                    end
                end
                any_sign_ok && return :compatible
            end
        end
        return saw_comparable ? :incompatible : :no_comparable
    end

    for (p, cands) in sort(collect(results_by_prime), by = x -> x[1])
        p == anchor_prime && continue
        trial_signs = [nothing]
        orig_sign = nothing
        if branch_sign_getter !== nothing && branch_sign_setter !== nothing
            orig_sign = branch_sign_getter(p)
            trial_signs = orig_sign == 1 ? [1, -1] : [-1, 1]
        end
        status = compatibility_status(cands, p, trial_signs, orig_sign)
        if status == :incompatible
            push!(contradictory, p)
        elseif status == :no_comparable && verbose
            println("  precheck: no comparable candidates at prime $p; skip contradiction split")
        end
    end
    if verbose && !isempty(contradictory)
        println("  precheck: branch contradiction detected at primes $contradictory")
    end
    return contradictory
end

function Base.show(io::IO, m::ClassifiedMTC)
    FR_status = if m.F_values === nothing
        "(F,R)=none"
    else
        rep = m.verify_report
        if rep !== nothing && hasproperty(rep, :ok)
            "(F,R) roundtrip=$(rep.ok ? "✓" : "✗")"
        else
            "(F,R)=attached"
        end
    end
    fresh_str = isempty(m.fresh_primes) ? "no fresh" :
        (m.verify_fresh ? "fresh✓" : "fresh✗")
    n_desc = m.N == m.N_input ? string(m.N) : "$(m.N) [input=$(m.N_input)]"
    print(io, "ClassifiedMTC(N=$(n_desc), rank=$(m.rank), ",
          "sector=$(m.galois_sector), $(length(m.used_primes)) primes, ",
          "$fresh_str, $FR_status)")
end

# ============================================================
#  compute_FR_from_ST: exact Phase 4 over Q(ζ_N)
# ============================================================

function compute_FR_from_ST(Nijk::Array{Int,3};
                            context = nothing,
                            conductor = nothing,
                            N = nothing,
                            S = nothing,
                            T = nothing,
                            return_all::Bool = false,
                            primes::Vector{Int} = [101, 103, 107, 109],
                            verbose::Bool = false,
                            kwargs...)
    isempty(kwargs) || error("unsupported keyword arguments: $(collect(keys(kwargs)))")
    ctx = _default_context_from_kwargs(context = context, conductor = conductor, N = N)
    r = size(Nijk, 1)
    known_F = _known_pentagon_solution(Nijk, ctx)
    F_solutions = if known_F === nothing
        _, pentagon_eqs, nF = get_pentagon_system(Nijk, r)
        solve_pentagon_modular_crt(pentagon_eqs, nF;
                                   Nijk = Nijk,
                                   context = ctx,
                                   primes = primes,
                                   show_progress = verbose)
    else
        [known_F]
    end
    candidates = NamedTuple[]
    for F in F_solutions
        _, hex_eqs, nR = get_hexagon_system(Nijk, r, F; context = ctx)
        R_solutions = solve_hexagon_modular_crt(hex_eqs, nR;
                                                Nijk = Nijk,
                                                context = ctx,
                                                primes = primes,
                                                show_progress = verbose)
        for R in R_solutions
            push!(candidates, (F = F, R = R, report = nothing))
        end
    end
    isempty(candidates) && _phase4_removed_error()
    selected = candidates[1]
    return (F = selected.F,
            R = selected.R,
            report = selected.report,
            candidates = candidates)
end

# ============================================================
#  classify_from_group: CRT + (F, R) solve for one Galois sector
# ============================================================

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

    # -------- Phase 4 input: exact lift (S, T) to Q(ζ_N) --------
    rep = first(values(group))
    zeta_Fp = find_zeta_in_Fp(N, rep.p)
    S_cyc, T_cyc, Nijk = lift_mtc_candidate(rep, recon_S;
                                            d = quadratic_d,
                                            N = N,
                                            zeta_Fp = zeta_Fp,
                                            scale = scale_factor)
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
                                   dimension for naive O(n) Cayley
                                   sweep. Raise with caution (O(4)(F_p)
                                   has ~10¹¹ points).
- `search_mode::Symbol = :groebner`:
                                   Phase 2 block-U backend mode passed
                                   to `find_mtcs_at_prime` (`:groebner`
                                   or `:exhaustive`).
- `max_units_for_groebner::Int = typemax(Int)`:
                                   cap on fixed-unit Gröbner systems
                                   tried per stratum/prime in Phase 2.
- `groebner_allow_fallback::Bool = true`:
                                   when `search_mode=:groebner`,
                                   whether to fall back to enumeration
                                   if solver extraction is empty.
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
                                    groebner_allow_fallback::Bool = true,
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
        fr_result = compute_FR_from_ST(rep.Nijk;
                                       context = CyclotomicContext(rep.N),
                                       S = rep.S_cyclotomic,
                                       T = rep.T_cyclotomic,
                                       return_all = true,
                                       primes = primes,
                                       verbose = verbose)
        fr_result.F === nothing && error("Phase 4 could not produce any (F,R) solution for key=$key")
        isempty(fr_result.candidates) && error("Phase 4 produced no valid (F,R) candidates for key=$key")

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
