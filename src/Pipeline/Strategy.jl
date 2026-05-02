"""
Lightweight computation strategy estimates.
"""

"""
    estimate_search_complexity(N, rank, stratum)

Return a coarse finite-field search complexity estimate for planning.
"""
function estimate_search_complexity(N::Int, rank::Int, stratum)
    block_count = stratum isa Stratum ? length(stratum.multiplicities) : 1
    score = max(1, N) * max(1, rank)^2 * max(1, block_count)
    class = score <= 200 ? :small : score <= 1500 ? :medium : :large
    return (score = score, class = class, rank = rank, conductor = N,
            block_count = block_count)
end

"""
    estimate_fr_reconstruction_complexity(Nijk)

Return a coarse exact F/R reconstruction estimate from associator and braiding
slot counts.
"""
function estimate_fr_reconstruction_complexity(Nijk::Array{Int,3})
    r = size(Nijk, 1)
    nF = try
        _, _, n = get_pentagon_system(Nijk, r)
        n
    catch
        Base.count(!iszero, Nijk)
    end
    nR = _forward_r_var_count(Nijk)
    score = max(1, nF) * max(1, nR)
    class = score <= 64 ? :small : score <= 512 ? :medium : :large
    return (score = score, class = class, nF = nF, nR = nR, rank = r)
end

"""
    recommend_primes(N, rank; min_count = 4, window = 2000, start_from = 29)

Recommend admissible CRT primes for the given conductor.
"""
function recommend_primes(N::Int, rank::Int; min_count::Int = 4,
                          window::Int = 2000, start_from::Int = 29)
    count = max(min_count, rank <= 3 ? 4 : 6)
    return select_admissible_primes(N; min_count = count,
                                    window = window,
                                    start_from = start_from)
end

"""
    recommend_skip_FR(Nijk; threshold = :medium)

Recommend whether to skip exact F/R reconstruction for a fusion rule.
"""
function recommend_skip_FR(Nijk::Array{Int,3}; threshold = :medium)
    est = estimate_fr_reconstruction_complexity(Nijk)
    skip = threshold == :small ? est.class != :small :
           threshold == :medium ? est.class == :large :
           false
    return (skip_FR = skip, estimate = est,
            reason = skip ? "estimated exact F/R reconstruction complexity is $(est.class)" :
                     "estimated exact F/R reconstruction complexity is acceptable")
end
