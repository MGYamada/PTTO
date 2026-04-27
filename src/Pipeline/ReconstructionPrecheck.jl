"""
Branch-consistency prechecks for pipeline reconstruction.

These helpers detect incompatible square-root branches before Galois-aware
grouping and reconstruction.
"""

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
