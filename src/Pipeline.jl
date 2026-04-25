"""
    End-to-end pipeline driver `N → List[ClassifiedMTC]`.

Given a conductor `N`, `classify_mtcs_at_conductor` executes the full pipeline and returns a list of fully classified
MTCs, each carrying:

- the SL(2, ℤ/N) stratum (m_λ) decomposition (Phase 0/1),
- F_p-validated modular data (Phase 2),
- Galois-coherent, arithmetic S-matrix in ℤ[√d] (Phase 3),
- `(F, R)` symbols in ℂ satisfying pentagon + hexagon (Phase 4),
- a `VerifyReport` summarising residuals.

This module also provides two mid-level helpers:

- `compute_FR_from_ST(Nijk; ...)`: given a fusion tensor,
  solve for `(F, R)` via pentagon HC → hexagon HC.

- `classify_from_group(group, all_primes; ...)`: given a Galois-coherent
  group from `group_mtcs_galois_aware`, performs Phase 3 CRT + Phase 4
  `(F, R)` solve and returns a single `ClassifiedMTC`.

And one high-level convenience API:

- `classify_mtcs_auto(N; ...)`: automatically chooses
  `conductor_mode`, `scale_d`, `primes`, and `max_rank` from candidate
  lists, then runs `classify_mtcs_at_conductor`.

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
- `N`:              effective conductor used internally by the pipeline
- `N_input`:        user-requested conductor (for provenance)
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
- `verify_report`:  `VerifyReport` with pentagon/hexagon residuals
                    max-residuals. `nothing` if no match was found.

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
                         F::Union{Vector{ComplexF64}, Nothing},
                         R::Union{Vector{ComplexF64}, Nothing},
                         report::Union{VerifyReport, Nothing})
    return ClassifiedMTC(c.N, c.N_input, c.rank, c.stratum, c.Nijk,
                         c.S_Zsqrtd, c.scale_d, c.scale_factor,
                         c.used_primes, c.fresh_primes, c.verify_fresh,
                         c.S_complex, c.T_complex,
                         F, R, report, c.galois_sector)
end

# ============================================================
#  classify_mtcs_auto: user-friendly auto-parameter wrapper
# ============================================================

cyclotomic_requirement(scale_d::Int) =
    (scale_d == 2 || scale_d == 3) ? 24 :
    scale_d == 5 ? 5 : 1

function compute_effective_conductor(N::Int, scale_d::Int;
                                     conductor_mode::Symbol = :full_mtc)
    conductor_mode == :full_mtc ||
        error("conductor_mode=:$(conductor_mode) was removed in v0.5.0. Use :full_mtc.")
    return lcm(N, cyclotomic_requirement(scale_d))
end

"""
    classify_mtcs_auto(N::Int;
                       max_rank_candidates = [2, 3, 4, 5],
                       scale_d_candidates = [3, 5, 2],
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
users who do not want to manually specify `max_rank`, `primes`,
`scale_d`, and `conductor_mode`. It now performs a stage-wise search
over `d_candidates`, where each stage sets:

`N_eff_candidate = compute_effective_conductor(N, d)`.

For each previously unseen `N_eff_candidate`, the driver tries
`(conductor_mode, scale_d, max_rank)` combinations at each effective
conductor and records stage metadata. Search stops when any of:

- no new MTCs for `stagnation_k` consecutive executed stages,
- `N_eff_candidate > N_eff_max`,
- total attempted runs reaches `max_attempts`.

Returns:

1. `classified`: deduplicated union of MTCs found across all stages.
2. reproducibility metadata:
   - `N_input`
   - `N_effective`
   - `d`
   - `scale_d`
   - `conductor_mode`
   - `primes`
   - `max_rank`
   - `attempts`
   - `history`

Prime selection rule for each attempted `(conductor_mode, scale_d)`:
- Let `req = lcm(N_effective, cyclotomic_requirement(scale_d))`, where
  `cyclotomic_requirement(2|3)=24`, `cyclotomic_requirement(5)=5`, else
  `1`.
- Choose the first `min_primes` primes `p > prime_start` with
  `(p - 1) % req == 0`.
"""
function classify_mtcs_auto(N::Int;
                            max_rank_candidates::Vector{Int} = [2, 3, 4, 5],
                            scale_d_candidates::Vector{Int} = [3, 5, 2],
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
    !isempty(scale_d_candidates) || error("scale_d_candidates must be non-empty")
    !isempty(d_candidates) || error("d_candidates must be non-empty")
    !isempty(conductor_modes) || error("conductor_modes must be non-empty")
    stagnation_k >= 1 || error("stagnation_k must be ≥ 1, got $stagnation_k")
    max_attempts >= 1 || error("max_attempts must be ≥ 1, got $max_attempts")

    last_result = ClassifiedMTC[]
    last_meta = (N_input = N, N_effective = N, d = 1, scale_d = 0,
                 conductor_mode = :full_mtc, primes = Int[],
                 max_rank = 0, attempts = 0)
    attempts = 0
    n_stagnant = 0

    seen_N_eff = Set{Int}()
    seen_signatures = Set{String}()
    history = NamedTuple[]

    mtc_signature(m::ClassifiedMTC) = begin
        nijk_key = join(vec(m.Nijk), ",")
        t_key = join(["$(round(real(t), digits = 10)):$(round(imag(t), digits = 10))"
                      for t in m.T_complex], ",")
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
        if N_eff_candidate in seen_N_eff
            push!(history, (d = d, N_effective = N_eff_candidate, executed = false,
                            success = false, reason = "duplicate_N_effective",
                            attempts = 0, new_mtcs = 0))
            continue
        end
        push!(seen_N_eff, N_eff_candidate)

        stage_attempts = 0
        stage_success = false
        stage_reason = "no_attempts"
        stage_new_mtcs = 0
        stage_result = ClassifiedMTC[]

        for conductor_mode in conductor_modes
            attempts >= max_attempts && break
            conductor_mode == :full_mtc || error(
                "conductor_mode=:$(conductor_mode) was removed in v0.5.0. Use :full_mtc.")

            for scale_d in scale_d_candidates
                attempts >= max_attempts && break

                req = lcm(N_eff_candidate, cyclotomic_requirement(scale_d))

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
                                       "d=$d N_eff=$N_eff_candidate " *
                                       "mode=$conductor_mode scale_d=$scale_d " *
                                       "max_rank=$max_rank primes=$chosen_primes")

                    classified = classify_mtcs_at_conductor(N_eff_candidate;
                                                            max_rank = max_rank,
                                                            primes = chosen_primes,
                                                            strata = strata,
                                                            scale_d = scale_d,
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

                    last_meta = (N_input = N, N_effective = N_eff_candidate, d = d,
                                 scale_d = scale_d, conductor_mode = conductor_mode,
                                 primes = copy(chosen_primes), max_rank = max_rank,
                                 attempts = attempts)
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
            scale_d = last_meta.scale_d,
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

# OBSOLETE (Phase 5):
# Directly calling this function for candidate comparison is deprecated.
# Use `_score_fr_st_match` / `_select_fr_for_st`, which centralize
# comparison rules and fallback behavior.
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

_maxabs(M::AbstractArray{<:Number}) = isempty(M) ? 0.0 : maximum(abs.(M))

function _normalize_twists(T::Vector{ComplexF64})
    Tn = copy(T)
    if !iszero(Tn[1])
        Tn ./= Tn[1]
    end
    for i in eachindex(Tn)
        if !iszero(Tn[i])
            Tn[i] /= abs(Tn[i])
        end
    end
    return Tn
end

function _normalize_smatrix(S::Matrix{ComplexF64})
    Sn = copy(S)
    if !iszero(Sn[1, 1])
        Sn ./= Sn[1, 1]
    end
    return Sn
end

function _infer_T_candidates_from_monodromy(R_values::Vector{ComplexF64},
                                            Nijk::Array{Int,3},
                                            N::Int;
                                            max_candidates::Int = 64)
    r = size(Nijk, 1)
    vars = max(0, r - 1)
    vars == 0 && return [ComplexF64[1.0 + 0.0im]]

    # equations: t_a + t_b - t_c = m (mod N), where
    # exp(2πim/N) approximates monodromy M^{ab}_c = R^{ba}_c R^{ab}_c.
    eqs = Tuple{Int, Int, Int, Int}[]
    for a in 1:r, b in 1:r, c in 1:r
        Nijk[a, b, c] == 0 && continue
        R_ab = extract_R_block(R_values, Nijk, a, b, c)
        R_ba = extract_R_block(R_values, Nijk, b, a, c)
        size(R_ab) == (1, 1) || continue
        size(R_ba) == (1, 1) || continue
        z = R_ba[1, 1] * R_ab[1, 1]
        z = iszero(z) ? z : z / abs(z)
        m = mod(round(Int, (N * angle(z)) / (2π)), N)
        push!(eqs, (a, b, c, m))
    end
    isempty(eqs) && return [vcat(ComplexF64[1.0 + 0.0im], fill(1.0 + 0.0im, r - 1))]

    sols = Vector{Vector{Int}}()
    ranges = ntuple(_ -> 0:(N - 1), vars)
    for tup in Iterators.product(ranges...)
        length(sols) >= max_candidates && break
        assignments = collect(tup)
        ok = true
        for (a, b, c, m) in eqs
            ta = a == 1 ? 0 : assignments[a - 1]
            tb = b == 1 ? 0 : assignments[b - 1]
            tc = c == 1 ? 0 : assignments[c - 1]
            if mod(ta + tb - tc - m, N) != 0
                ok = false
                break
            end
        end
        ok && push!(sols, assignments)
    end

    T_candidates = Vector{Vector{ComplexF64}}()
    for sol in sols
        Tvals = Vector{ComplexF64}(undef, r)
        Tvals[1] = 1.0 + 0.0im
        for i in 2:r
            Tvals[i] = exp((2π * im * sol[i - 1]) / N)
        end
        push!(T_candidates, Tvals)
    end
    isempty(T_candidates) && push!(T_candidates, vcat(ComplexF64[1.0 + 0.0im], fill(1.0 + 0.0im, r - 1)))
    return T_candidates
end

function _st_signature_under_gauge(S::Matrix{ComplexF64},
                                   T::Vector{ComplexF64},
                                   automorphisms::Vector{Vector{Int}};
                                   digits::Int = 8)
    best = nothing
    for perm in automorphisms
        S_p = S[perm, perm]
        T_p = T[perm]
        for sgn in (1.0, -1.0)
            S_use = sgn > 0 ? S_p : -S_p
            words = String[]
            push!(words, join([string(round(real(z), digits = digits), ",",
                                      round(imag(z), digits = digits)) for z in T_p], ";"))
            push!(words, join([string(round(real(z), digits = digits), ",",
                                      round(imag(z), digits = digits)) for z in vec(S_use)], ";"))
            sig = join(words, "|")
            if best === nothing || sig < best
                best = sig
            end
        end
    end
    return best
end

"""
    _modular_data_roundtrip(F_values, R_values, Nijk, S_target, T_target, N) -> NamedTuple

Reconstruct `(S,T)` from `(F,R)` using the note's balancing/Hopf-link form:
- infer twist candidates from monodromy phases
  `M^{ab}_c = R^{ba}_c R^{ab}_c = θ_a θ_b / θ_c`
  (multiplicity-free channels),
- rebuild `S` from
  `S_ab = (1/D) Σ_c N_ab^c d_c (θ_a θ_b / θ_c)`.

Then compare reconstructed data against targets while accounting for
non-Galois gauge freedom (fusion-rule automorphisms fixing the unit and
the global `S ↦ -S` convention). We also allow complex-conjugate
convention matching `(S,T) ↔ (conj(S), conj(T))`.
"""
function _modular_data_roundtrip(F_values::Vector{ComplexF64},
                                 R_values::Vector{ComplexF64},
                                 Nijk::Array{Int,3},
                                 S_target::Matrix{ComplexF64},
                                 T_target::Vector{ComplexF64},
                                 N::Int)
    _ = F_values
    r = size(Nijk, 1)

    # Quantum dimensions from PF eigenvector of Σ_i N_i.
    A = zeros(Float64, r, r)
    for i in 1:r
        A .+= Float64.(Nijk[i, :, :])
    end
    eig = eigen(A)
    idx = argmax(real(eig.values))
    d = abs.(real(eig.vectors[:, idx]))
    d ./= d[1]
    D = sqrt(sum(d .^ 2))

    T_target_n = _normalize_twists(T_target)
    S_target_n = _normalize_smatrix(S_target)

    function reconstruct_S_from_T(Tvals::Vector{ComplexF64})
        S = Matrix{ComplexF64}(undef, r, r)
        for a in 1:r, b in 1:r
            acc = 0.0 + 0.0im
            for c in 1:r
                Nijk[a, b, c] == 0 && continue
                # Balancing form:
                #   R^{ba}_c R^{ab}_c = θ_a θ_b / θ_c
                # so S_ab = (1/D) Σ_c N_ab^c d_c (θ_a θ_b / θ_c).
                acc += Nijk[a, b, c] * (Tvals[a] * Tvals[b] / Tvals[c]) * d[c]
            end
            S[a, b] = acc / D
        end
        return S
    end

    automorphisms = _fusion_automorphisms_fixing_unit(Nijk)
    best_perm = automorphisms[1]
    best_s = Inf
    best_t = Inf
    best_Sn = _normalize_smatrix(reconstruct_S_from_T(T_target_n))
    best_Tn = copy(T_target_n)

    for T_raw in _infer_T_candidates_from_monodromy(R_values, Nijk, N)
        Tn = _normalize_twists(T_raw)
        Sn = _normalize_smatrix(reconstruct_S_from_T(Tn))
        for perm in automorphisms
            S_p = Sn[perm, perm]
            T_p = Tn[perm]
            for use_conj in (false, true)
                S_cmp = use_conj ? conj.(S_p) : S_p
                T_cmp = use_conj ? conj.(T_p) : T_p
                s_err = min(_maxabs(S_cmp .- S_target_n), _maxabs(-S_cmp .- S_target_n))
                t_err = _maxabs(T_cmp .- T_target_n)
                if (s_err < best_s) || (isapprox(s_err, best_s; atol = 1e-12) && t_err < best_t)
                    best_s = s_err
                    best_t = t_err
                    best_perm = perm
                    best_Sn = use_conj ? conj.(Sn) : Sn
                    best_Tn = use_conj ? conj.(Tn) : Tn
                end
            end
        end
    end

    ok = best_s < 5e-2 && best_t < 5e-2
    return (ok = ok,
            S_max = best_s,
            T_max = best_t,
            best_perm = best_perm,
            S_from = best_Sn,
            T_from = best_Tn,
            st_signature = _st_signature_under_gauge(best_Sn, best_Tn, automorphisms))
end

"""
    _score_fr_st_match(F_values, R_values, Nijk, S_target, T_target, N;
                       candidate_index)
        -> NamedTuple

Phase-5 API to score `(F,R)` against `(S,T)`.

Returns:
- `ok`: whether roundtrip errors are within threshold
- `S_max`: maximum S-side error
- `T_max`: maximum T-side error
- `order_key`: key for deterministic total ordering
- `candidate_index`: input candidate index (final tie-break)
"""
function _score_fr_st_match(F_values::Vector{ComplexF64},
                            R_values::Vector{ComplexF64},
                            Nijk::Array{Int,3},
                            S_target::Matrix{ComplexF64},
                            T_target::Vector{ComplexF64},
                            N::Int;
                            candidate_index::Int)
    md = _modular_data_roundtrip(F_values, R_values, Nijk, S_target, T_target, N)
    # Fixed comparison rule for Phase 5:
    #   1) prioritize ok=true
    #   2) smaller S_max
    #   3) smaller T_max
    #   4) smaller candidate_index
    # Since false < true for Bool, place !ok first to prefer ok=true.
    order_key = (!md.ok, md.S_max, md.T_max, candidate_index)
    return (ok = md.ok,
            S_max = md.S_max,
            T_max = md.T_max,
            best_perm = md.best_perm,
            st_signature = md.st_signature,
            order_key = order_key,
            candidate_index = candidate_index)
end

function _build_fr_st_dictionary(candidates::Vector{<:NamedTuple},
                                 Nijk::Array{Int,3},
                                 N::Int)
    table = Dict{String, Vector{Int}}()
    reconstructed = NamedTuple[]
    for (ci, cand) in enumerate(candidates)
        md = _modular_data_roundtrip(cand.F, cand.R, Nijk,
                                     zeros(ComplexF64, size(Nijk, 1), size(Nijk, 1)),
                                     ones(ComplexF64, size(Nijk, 1)),
                                     N)
        push!(get!(table, md.st_signature, Int[]), ci)
        push!(reconstructed, (candidate_index = ci,
                              signature = md.st_signature,
                              S = md.S_from,
                              T = md.T_from))
    end
    return (dictionary = table, reconstructed = reconstructed)
end

"""
    _select_fr_for_st(candidates, Nijk, S_target, T_target, N)
        -> (selected, score, selected_index, selected_ok, all_scores)

Phase-5 candidate selection API.
Selects one `(F,R)` candidate using total ordering by
`order_key = (!ok, S_max, T_max, candidate_index)`.

Fallback policy:
- If at least one candidate has `ok=true`, choose the minimal `order_key`.
- If none has `ok=true`, do not raise an explicit error; choose the
  minimum-error candidate (`S_max/T_max`, then index) and expose
  `selected_ok=false` to callers.
"""
function _select_fr_for_st(candidates::Vector{<:NamedTuple},
                           Nijk::Array{Int,3},
                           S_target::Matrix{ComplexF64},
                           T_target::Vector{ComplexF64},
                           N::Int)
    isempty(candidates) && error("cannot select (F,R): empty candidate list")
    fr_st = _build_fr_st_dictionary(candidates, Nijk, N)

    all_scores = NamedTuple[]
    best_idx = 0
    best_score = nothing
    for (ci, cand) in enumerate(candidates)
        score = _score_fr_st_match(cand.F, cand.R, Nijk, S_target, T_target, N;
                                   candidate_index = ci)
        push!(all_scores, score)
        if best_score === nothing || score.order_key < best_score.order_key
            best_score = score
            best_idx = ci
        end
    end

    return (selected = candidates[best_idx],
            score = best_score,
            selected_index = best_idx,
            selected_ok = best_score.ok,
            all_scores = all_scores,
            fr_st_dictionary = fr_st.dictionary)
end

function _branch_consistency_precheck(results_by_prime::Dict{Int, Vector{MTCCandidate}},
                                      anchor_prime::Int,
                                      scale_d::Int,
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
                        s_anchor = sqrtd_fn(scale_d, anchor_prime)
                        s_p = sqrtd_fn(scale_d, p)
                        two_s_anchor = mod(2 * s_anchor, anchor_prime)
                        two_s_p = mod(2 * s_p, p)
                        matrix_by_prime = Dict(
                            anchor_prime => [mod(two_s_anchor * anchor_c.S_Fp[i, j], anchor_prime)
                                             for i in 1:nrow, j in 1:ncol],
                            p => [mod(two_s_p * c.S_Fp[i, j], p)
                                  for i in 1:nrow, j in 1:ncol])
                        reconstruct_matrix_in_Z_sqrt_d(matrix_by_prime, scale_d;
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
        pent = rep === nothing ? "?" : string(round(rep.pentagon_max, sigdigits = 2))
        hex = rep === nothing ? "?" : string(round(rep.hexagon_max, sigdigits = 2))
        "(F,R) pent=$pent hex=$hex"
    end
    fresh_str = isempty(m.fresh_primes) ? "no fresh" :
        (m.verify_fresh ? "fresh✓" : "fresh✗")
    n_desc = m.N == m.N_input ? string(m.N) : "$(m.N) [input=$(m.N_input)]"
    print(io, "ClassifiedMTC(N=$(n_desc), rank=$(m.rank), ",
          "sector=$(m.galois_sector), $(length(m.used_primes)) primes, ",
          "$fresh_str, $FR_status)")
end

# ============================================================
#  compute_FR_from_ST: Phase 4 wrapper
# ============================================================

"""
    compute_FR_from_ST(Nijk;
                        return_all = false,
                        pentagon_slice = 1, show_progress = false,
                        verbose = false)
        -> NamedTuple{(:F, :R, :report, :n_pentagon, :n_tried, :n_matches,
                       :f_idx, :r_idx)}

Given a fusion tensor `Nijk`, find a pair `(F, R)` of complex
F- and R-symbols satisfying pentagon and hexagon.

This stage solves pentagon/hexagon from fusion data only. Branch
selection against reconstructed modular data is done afterwards by
`_modular_data_roundtrip` in `classify_from_group`.

Algorithm:
1. Set up pentagon system from `Nijk` (via TensorCategories).
2. Solve pentagon via homotopy continuation, optionally with random
   linear slice to break gauge symmetry.
3. For each pentagon solution `F`:
   - Polish with damped Newton.
   - Build hexagon system with F fixed.
   - Solve hexagon via HC.
   - Keep the `(F,R)` pair with the smallest pentagon/hexagon residual
     score.
4. Return the best numerical pentagon/hexagon solution.

Caveats:
- Pentagon HC is feasible only for small fusion rings (~10 F-variables
  before mixed-volume blow-up; Ising-sized rings with ~14 F-vars may
  require hours).
- Multiple inequivalent braided branches can share the same fusion ring.
  This function does not disambiguate them via `T`; modular-data
  equivalence is handled in the Phase 4 roundtrip check.

Returns a NamedTuple with:
- `F`:            best `Vector{ComplexF64}` or `nothing`
- `R`:            best `Vector{ComplexF64}` or `nothing`
- `report`:       `VerifyReport` or `nothing`
- `n_pentagon`:   number of pentagon solutions found
- `n_tried`:      total `(F, R)` pairs examined
- `n_matches`:    how many candidates were numerically valid
- `f_idx`, `r_idx`: indices of the chosen pair (0 if none)
- `candidates`:   present only when `return_all=true`; contains all
                  numerically valid `(F,R)` candidates with reports.

"""
function compute_FR_from_ST(Nijk::Array{Int, 3};
                            return_all::Bool = false,
                            pentagon_slice::Int = 1,
                            show_progress::Bool = false,
                            verbose::Bool = false)
    r = size(Nijk, 1)

    # Rank-1 (trivial) MTC: the only fusion category structure is the
    # unit object with F = R = [1]. Pentagon / hexagon hold vacuously.
    # `TensorCategories.pentagon_equations` emits no non-
    # trivial polynomials at rank 1, so we short-circuit here rather
    # than let `get_pentagon_system` error.
    if r == 1
        F_trivial = ComplexF64[1.0]
        R_trivial = ComplexF64[1.0]
        report = VerifyReport(0.0, 0.0, 0, 0, 1)
        cand = (F = F_trivial, R = R_trivial, report = report,
                f_idx = 1, r_idx = 1)
        if return_all
            return (F = F_trivial, R = R_trivial, report = report,
                    n_pentagon = 1, n_tried = 1, n_matches = 1,
                    f_idx = 1, r_idx = 1,
                    candidates = [cand])
        end
        return (F = F_trivial, R = R_trivial, report = report,
                n_pentagon = 1, n_tried = 1, n_matches = 1,
                f_idx = 1, r_idx = 1)
    end

    # Pentagon system
    local _R, eqs, n
    try
        _R, eqs, n = get_pentagon_system(Nijk, r)
    catch err
        msg = sprint(showerror, err)
        if occursin("Object not rigid", msg)
            error("fusion ring is not rigid at rank $r; rejecting candidate")
        end
        rethrow(err)
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
            score = Inf)
    all_candidates = NamedTuple{(:F, :R, :report, :f_idx, :r_idx),
                                Tuple{Vector{ComplexF64}, Vector{ComplexF64},
                                      VerifyReport, Int, Int}}[]

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

        for (ri, R_vals) in enumerate(R_sols)
            best = (; best..., n_tried = best.n_tried + 1)

            local rep
            try
                # First pass verifies pentagon/hexagon consistency.
                rep = verify_mtc(F, R_vals, Nijk)
            catch err
                verbose && println("      R[$ri] verify failed: $err")
                continue
            end
            score = max(rep.pentagon_max, rep.hexagon_max)
            if !isfinite(score)
                verbose && println("      R[$ri] non-finite pent/hex score: $score")
                continue
            end

            # Keep counters for API compatibility. n_matches now counts
            # numerically valid pentagon+hexagon candidates.
            best = (; best..., n_matches = best.n_matches + 1)
            push!(all_candidates, (F = F, R = R_vals, report = rep,
                                   f_idx = fi, r_idx = ri))
            if score < best.score
                best = (; best...,
                        F = F, R = R_vals, report = rep,
                        f_idx = fi, r_idx = ri,
                        score = score)
            end
        end
    end

    # Drop the internal score field from the public return
    if return_all
        return (F = best.F, R = best.R, report = best.report,
                n_pentagon = best.n_pentagon, n_tried = best.n_tried,
                n_matches = best.n_matches,
                f_idx = best.f_idx, r_idx = best.r_idx,
                candidates = all_candidates)
    end
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
                        reconstruction_bound = 50,
                        galois_sector = 1,
                        test_primes = nothing,
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
- `reconstruction_bound::Int = 50`:  ℤ[√d] coefficient bound for rational
                                     reconstruction.
- `galois_sector::Int = 1`:          sector index for provenance.
- `test_primes::Vector{Int}=nothing`: which primes to use for CRT.
                                     If `nothing`, uses first
                                     `max(2, length(all_primes) ÷ 2)`
                                     primes; remaining are fresh.
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
                             N_input::Int = N,
                             scale_d::Int,
                             scale_factor::Int = 2,
                             sqrtd_fn = compute_sqrt_d_mod_p,
                             reconstruction_bound::Int = 50,
                             galois_sector::Int = 1,
                             test_primes::Union{Vector{Int}, Nothing} = nothing,
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
    recon_S_phase4 = recon_S

    # TensorCategories assumes the unit object is at index 1.
    # MTCCandidate stores `unit_index` explicitly and may keep a different
    # basis ordering; permute all lifted data coherently before Phase 4.
    if rep.unit_index != 1
        perm = vcat(rep.unit_index, [i for i in 1:length(T_ℂ) if i != rep.unit_index])
        S_ℂ = S_ℂ[perm, perm]
        T_ℂ = T_ℂ[perm]
        Nijk = Nijk[perm, perm, perm]
        recon_S_phase4 = recon_S[perm, perm]
    end
    T_for_phase4 = T_ℂ

    # -------- Phase 4: (F, R) solve + verify --------
    rank = size(Nijk, 1)

    if skip_FR
        return ClassifiedMTC(N, N_input, rank, stratum, Nijk,
                             recon_S_phase4, scale_d, scale_factor,
                             used, fresh, verify_fresh,
                             S_ℂ, T_for_phase4, nothing, nothing, nothing,
                             galois_sector)
    end
    verbose && println("  running pentagon/hexagon on rank=$rank...")

    fr_result = compute_FR_from_ST(Nijk;
                                   return_all = true,
                                   verbose = verbose)
    fr_result.F === nothing && error("Phase 4 could not produce any (F,R) solution")
    isempty(fr_result.candidates) && error("Phase 4 produced no valid (F,R) candidates")

    # OBSOLETE: candidate loop with direct `_modular_data_roundtrip`
    # calls is replaced by Phase 5 API `_select_fr_for_st`.
    selection = _select_fr_for_st(fr_result.candidates, Nijk, S_ℂ, T_for_phase4, N)
    selected = selection.selected
    md_roundtrip = selection.score
    verbose && println("  modular-data roundtrip: perm=$(md_roundtrip.best_perm), " *
                       "S_err=$(md_roundtrip.S_max), T_err=$(md_roundtrip.T_max), " *
                       "ok=$(md_roundtrip.ok)")

    return ClassifiedMTC(N, N_input, rank, stratum, Nijk, recon_S_phase4,
                         scale_d, scale_factor,
                         used, fresh, verify_fresh,
                         S_ℂ, T_for_phase4,
                         selected.F, selected.R, selected.report,
                         galois_sector)
end

# ============================================================
#  classify_mtcs_at_conductor: top-level driver
# ============================================================

"""
    classify_mtcs_at_conductor(N; max_rank = 5, primes = nothing, strata = nothing,
                                scale_d = 3, scale_factor = 2,
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
  Phase 4: lift to ℂ, solve pentagon/hexagon

Returns one `ClassifiedMTC` per Galois sector per stratum that yields a
valid MTC. If `skip_FR = true`, Phase 4 is skipped and `(F, R)` is left
as `nothing` (useful for large ranks where pentagon HC is infeasible).

Note on conductor: `N` is an input indicator (typically the T-order
conductor). Internally, the pipeline searches at `N_effective`
determined by `conductor_mode`; `:full_mtc` uses
`N_effective = lcm(N, cyclotomic_requirement(scale_d))` where
`cyclotomic_requirement(2|3)=24`, `cyclotomic_requirement(5)=5`, else `1`.
An MTC's S-matrix may
live in a larger cyclotomic field than ℚ(ζ_N); in NRWW's convention the
full MTC conductor is `max(cond(S), cond(T))`. Fibonacci, for example,
has `cond(T) = 5` but its S involves `D = √(2+φ)` which is NOT in ℚ(ζ_5).
Using only the T-order conductor can miss such MTCs — they often appear
only when `N_effective` is large enough to accommodate `S` too.

Arguments:
- `N::Int`:                        input conductor indicator. Internal
                                   search runs at `N_effective`.
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
- `conductor_mode::Symbol = :full_mtc`:
                                   interpretation of `N`. In v0.5.0,
                                   only `:full_mtc` is supported, using
                                   `N_effective = lcm(N, cyclotomic_requirement(scale_d))`
                                   (typically equal to `N`) so
                                   S-side field constraints are
                                   conservatively included.
- `min_primes::Int = 4`:           when `primes = nothing`, number of
                                   admissible primes to auto-select.
- `prime_start::Int = 29`:         when `primes = nothing`, start of
                                   prime search range (exclusive).
- `prime_window::Int = 2000`:      when `primes = nothing`, scan width.
- `sqrtd_fn`:                      custom √d-in-F_p function. If
                                   `nothing` (default), chooses:
                                   cyclotomic variant for `scale_d ∈
                                   {2,3,5}` and anchored mode for
                                   other `scale_d` (anchor prime +
                                   branch-transform layer over raw
                                   Tonelli roots):
                                   `compute_sqrt3_cyclotomic_mod_p`
                                   (needs 24 | p-1) for `scale_d = 3`,
                                   `compute_sqrt2_cyclotomic_mod_p`
                                   (needs 24 | p-1) for `scale_d = 2`,
                                   `compute_sqrt5_cyclotomic_mod_p`
                                   (needs 5 | p-1) for `scale_d = 5`,
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
- `skip_FR::Bool = false`:         skip Phase 4. Useful if
                                   `max_rank ≥ 5` and pentagon HC would
                                   blow up.
- `verbose::Bool = true`:          print per-phase progress.
"""
function classify_mtcs_at_conductor(N::Int;
                                    max_rank::Int = 5,
                                    primes::Union{Nothing, Vector{Int}} = nothing,
                                    strata::Union{Nothing, Vector{Stratum}} = nothing,
                                    scale_d::Int = 3,
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
                                    verbose::Bool = true)
    user_sqrtd_fn = sqrtd_fn
    N_effective = compute_effective_conductor(N, scale_d;
                                              conductor_mode = conductor_mode)

    chosen_primes = primes === nothing ?
        select_admissible_primes(N_effective;
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
            selector = build_sqrtd_selector(scale_d, present_primes, anchor; verbose = verbose)
            active_sqrtd_fn = selector.sqrtd_fn
            selector_mode = selector.mode
            branch_sign_getter = selector.branch_sign_getter
            branch_sign_setter = selector.branch_sign_setter
        end
        verbose && println("  d=$scale_d, sqrt-branch mode=$(selector_mode == :custom ? "custom" : String(selector_mode)), grouping bound=$reconstruction_bound")

        contradictory = _branch_consistency_precheck(results_by_prime, anchor, scale_d, active_sqrtd_fn;
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
                                                        scale_d = scale_d,
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
                                                scale_d = scale_d,
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

    # ------- Phase 4: aggregate by fusion rule and run (F,R) once per key -------
    verbose && println("\n=== Phase 4: fusion-rule aggregation + (F, R) classification ===")
    grouped = _classify_modular_data_by_fusion_rule(out)
    key_to_members = Dict{String, Vector{Tuple{Matrix{ComplexF64}, Vector{ComplexF64}}}}()
    for (key, idxs) in grouped
        key_to_members[key] = [(out[i].S_complex, out[i].T_complex) for i in idxs]
    end
    verbose && println("  modular-data groups: $(length(out)) (S,T) results → " *
                       "$(length(grouped)) fusion-rule keys")

    for (key, idxs) in grouped
        rep_idx = idxs[1]
        rep = out[rep_idx]
        verbose && println("  key=$key: members=$(length(idxs)), rank=$(rep.rank)")
        fr_result = compute_FR_from_ST(rep.Nijk;
                                       return_all = true,
                                       verbose = verbose)
        fr_result.F === nothing && error("Phase 4 could not produce any (F,R) solution for key=$key")
        isempty(fr_result.candidates) && error("Phase 4 produced no valid (F,R) candidates for key=$key")

        # OBSOLETE: old behavior selected using only the fusion-rule
        # representative (`rep`). New behavior assigns `(F,R)` per `(S,T)`
        # member within the same fusion-rule group.
        for i in idxs
            mtc_i = out[i]
            selection = _select_fr_for_st(fr_result.candidates, mtc_i.Nijk,
                                          mtc_i.S_complex, mtc_i.T_complex, mtc_i.N)
            selected = selection.selected
            best_md = selection.score
            verbose && println("    member[$i] branch: cand=$(selection.selected_index), " *
                               "perm=$(best_md.best_perm), " *
                               "S_err=$(best_md.S_max), T_err=$(best_md.T_max), ok=$(best_md.ok)")
            out[i] = _with_fr_result(out[i], selected.F, selected.R, selected.report)
        end
    end

    verbose && println("\n=== Done: $(length(out)) ClassifiedMTC(s) ===")
    return out
end
