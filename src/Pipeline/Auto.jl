"""
Auto-parameter wrapper for conductor-first classification.

This file contains conductor helpers and the public `classify_mtcs_auto`
convenience API.  The mathematical input conductor is kept fixed: ACMG treats
`N` as the cyclotomic level through which the modular representation
`SL(2, ℤ) -> GL_r` factors.
"""

# ============================================================
#  classify_mtcs_auto: user-friendly auto-parameter wrapper
# ============================================================

function compute_effective_conductor(N::Int, args...;
                                     conductor_mode::Symbol = :full_mtc)
    conductor_mode == :full_mtc ||
        error("unsupported conductor_mode=:$(conductor_mode); use :full_mtc")
    return N
end

"""
    classify_mtcs_auto(N::Int;
                       max_rank_candidates = [2, 3, 4, 5],
                       conductor_modes = [:full_mtc],
                       min_primes = 4,
                       prime_start = 29,
                       prime_max = 2000,
                       N_eff_max = typemax(Int),
                       stagnation_k = 3,
                       max_attempts = typemax(Int),
                       strata = nothing,
                       scale_factor = 2,
                       verlinde_threshold = 3,
                       max_block_dim = 3,
                       search_mode = :groebner,
                       max_units_for_groebner = typemax(Int),
                       groebner_allow_fallback = false,
                       precheck_unit_axiom = true,
                       reconstruction_bound = 50,
                       skip_FR = false,
                       verbose = true)
        -> NamedTuple

Auto-select wrapper around `classify_mtcs_at_conductor`.

This function is intended as the recommended public entry point for
users who do not want to manually specify a rank cutoff and admissible
finite-field primes.  The conductor is fixed to the user-supplied `N`.

Mathematically, the search is over modular data whose `SL(2, ℤ)` representation
factors through `SL(2, ℤ/N)`.  For each requested rank cutoff, ACMG builds the
atomic irreducible representation catalog, enumerates strata, searches the
finite-field Block-U equations at admissible primes, and reconstructs exact
cyclotomic modular data when the residues are consistent.

The driver records stage metadata for reproducibility.  Search stops when any
of:

- `N_eff_candidate > N_eff_max`,
- total attempted runs reaches `max_attempts`.

`stagnation_k` is accepted as a compatibility keyword.

Returns:

1. `classified`: deduplicated union of MTCs found across all stages.
2. reproducibility metadata:
   - `N_input`
   - `N_effective`
   - `conductor_mode`
   - `primes`
   - `max_rank`
   - `attempts`
   - `history`

Prime selection chooses the first `min_primes` primes `p > prime_start`
with `N | (p - 1)`, so that the `N`-th roots of unity split in `F_p`.

`reconstruction_bound` is forwarded as the requested cyclotomic CRT
fallback bound. The pipeline applies a small internal cap before MITM
reconstruction to avoid explosive searches in degenerate strata.
"""
function classify_mtcs_auto(N::Int;
                            max_rank_candidates::Vector{Int} = [2, 3, 4, 5],
                            conductor_modes::Vector{Symbol} = [:full_mtc],
                            min_primes::Int = 4,
                            prime_start::Int = 29,
                            prime_max::Int = 2000,
                            N_eff_max::Int = typemax(Int),
                            stagnation_k::Int = 3,
                            max_attempts::Int = typemax(Int),
                            strata::Union{Nothing, Vector{Stratum}} = nothing,
                            scale_factor::Int = 2,
                            verlinde_threshold::Int = 3,
                            max_block_dim::Int = 3,
                            search_mode::Symbol = :groebner,
                            max_units_for_groebner::Int = typemax(Int),
                            groebner_allow_fallback::Bool = false,
                            precheck_unit_axiom::Bool = true,
                            reconstruction_bound::Int = 50,
                            skip_FR::Bool = false,
                            verbose::Bool = true)
    N >= 1 || error("N must be positive, got $N")
    min_primes >= 2 || error("min_primes must be ≥ 2, got $min_primes")
    !isempty(max_rank_candidates) || error("max_rank_candidates must be non-empty")
    !isempty(conductor_modes) || error("conductor_modes must be non-empty")
    stagnation_k >= 1 || error("stagnation_k must be ≥ 1, got $stagnation_k")
    max_attempts >= 1 || error("max_attempts must be ≥ 1, got $max_attempts")

    last_result = ClassifiedMTC[]
    last_meta = (N_input = N, N_effective = N,
                 conductor_mode = :full_mtc, primes = Int[],
                 max_rank = 0, attempts = 0)
    attempts = 0

    seen_signatures = Set{String}()
    history = NamedTuple[]

    mtc_signature(m::ClassifiedMTC) = begin
        nijk_key = join(vec(m.Nijk), ",")
        t_key = join(string.(m.T_cyclotomic), ",")
        string(m.rank, "|", m.galois_sector, "|", nijk_key, "|", t_key)
    end

    N_eff_candidate = compute_effective_conductor(N)
    if N_eff_candidate > N_eff_max
        push!(history, (N_effective = N_eff_candidate, executed = false,
                        success = false, reason = "N_eff_max_exceeded",
                        attempts = 0, new_mtcs = 0))
    else
        stage_attempts = 0
        stage_success = false
        stage_reason = "no_attempts"
        stage_new_mtcs = 0
        stage_result = ClassifiedMTC[]

        for conductor_mode in conductor_modes
            attempts >= max_attempts && break
            conductor_mode == :full_mtc || error(
                "unsupported conductor_mode=:$(conductor_mode); use :full_mtc")

            req = N_eff_candidate

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
                                       "N_eff=$N_eff_candidate req=$req " *
                                       "mode=$conductor_mode " *
                                       "max_rank=$max_rank primes=$chosen_primes")

                    last_meta = (N_input = N, N_effective = N_eff_candidate,
                                 conductor_mode = conductor_mode,
                                 primes = copy(chosen_primes), max_rank = max_rank,
                                 attempts = attempts)

                    local classified
                    try
                        classified = classify_mtcs_at_conductor(N_eff_candidate;
                                                                max_rank = max_rank,
                                                                primes = chosen_primes,
                                                                strata = strata,
                                                                scale_factor = scale_factor,
                                                                conductor_mode = conductor_mode,
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

        push!(history, (N_effective = N_eff_candidate, executed = true,
                        success = stage_success, reason = stage_reason,
                        attempts = stage_attempts, new_mtcs = stage_new_mtcs))

    end

    return (classified = last_result,
            N_input = last_meta.N_input,
            N_effective = last_meta.N_effective,
            conductor_mode = last_meta.conductor_mode,
            primes = last_meta.primes,
            max_rank = last_meta.max_rank,
            attempts = last_meta.attempts,
            history = history)
end
