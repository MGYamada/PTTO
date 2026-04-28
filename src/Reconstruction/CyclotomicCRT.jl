"""
Shared bounded modular CRT over the power basis of `Q(ζ_N)`.

This is the common lifting layer used by Phase 3 modular-data
reconstruction and by the exact F/R solvers.
"""

const _CRT_COEFF_VECTOR_CACHE = Dict{Tuple{Int, Int}, Vector{Vector{Int}}}()
const _CRT_MITM_LEFT_CACHE = Dict{Tuple{Int, Int, Tuple{Vararg{Int}}, Tuple{Vararg{Int}}}, Dict{Tuple{Vararg{Int}}, Vector{Int}}}()
const _CRT_MITM_RIGHT_CACHE = Dict{Tuple{Int, Int, Int, Tuple{Vararg{Int}}, Tuple{Vararg{Int}}}, Vector{Tuple{Vector{Int}, Vector{Int}}}}()

function _eval_power_basis_mod(coeffs::AbstractVector{<:Integer},
                               zeta_Fp::Int,
                               p::Int)
    acc = 0
    zpow = 1
    for c in coeffs
        acc = mod(acc + mod(c, p) * zpow, p)
        zpow = mod(zpow * zeta_Fp, p)
    end
    return acc
end

function _coeff_score(v::Vector{Int})
    abs_vals = [abs(x) for x in v]
    max_abs = isempty(abs_vals) ? 0 : maximum(abs_vals)
    return (sum(abs_vals), max_abs, v)
end

function _coeff_vectors(len::Int, bound::Int)
    key = (len, bound)
    haskey(_CRT_COEFF_VECTOR_CACHE, key) && return _CRT_COEFF_VECTOR_CACHE[key]
    len == 0 && return [Int[]]
    ranges = ntuple(_ -> -bound:bound, len)
    vectors = vec([collect(t) for t in Iterators.product(ranges...)])
    sort!(vectors; lt = (a, b) -> isless(_coeff_score(a), _coeff_score(b)))
    _CRT_COEFF_VECTOR_CACHE[key] = vectors
    return vectors
end

function _reconstruct_power_basis_coeffs_mitm(targets::Vector{Int},
                                              zetas::Vector{Int},
                                              primes::Vector{Int},
                                              degree_K::Int,
                                              coeff_bound::Int)
    left_len = degree_K ÷ 2
    right_len = degree_K - left_len
    left_vectors = _coeff_vectors(left_len, coeff_bound)
    right_vectors = _coeff_vectors(right_len, coeff_bound)

    zetas_key = Tuple(zetas)
    primes_key = Tuple(primes)
    left_key = (left_len, coeff_bound, zetas_key, primes_key)
    table = get!(_CRT_MITM_LEFT_CACHE, left_key) do
        built = Dict{Tuple{Vararg{Int}}, Vector{Int}}()
        for c_left in left_vectors
            key = Tuple(_eval_power_basis_mod(c_left, zetas[i], primes[i])
                        for i in eachindex(primes))
            if !haskey(built, key) || _coeff_score(c_left) < _coeff_score(built[key])
                built[key] = c_left
            end
        end
        built
    end

    right_key = (right_len, coeff_bound, left_len, zetas_key, primes_key)
    right_data = get!(_CRT_MITM_RIGHT_CACHE, right_key) do
        zeta_offsets = [powermod(zetas[i], left_len, primes[i]) for i in eachindex(primes)]
        data = Tuple{Vector{Int}, Vector{Int}}[]
        for c_right in right_vectors
            vals = Int[]
            for i in eachindex(primes)
                right_val = _eval_power_basis_mod(c_right, zetas[i], primes[i])
                push!(vals, mod(zeta_offsets[i] * right_val, primes[i]))
            end
            push!(data, (vals, c_right))
        end
        data
    end

    best = nothing
    best_score = nothing
    for (vals, c_right) in right_data
        need = Tuple(mod(targets[i] - vals[i], primes[i]) for i in eachindex(primes))
        if haskey(table, need)
            coeffs = vcat(table[need], c_right)
            score = _coeff_score(coeffs)
            if best === nothing || score < best_score
                best = coeffs
                best_score = score
            end
        end
    end
    return best
end

function _cyclotomic_element_from_coeffs(coeffs::Vector{Int}, denom::Int,
                                         ctx::CyclotomicContext)
    K = field(ctx)
    z = zeta(ctx)
    x = zero(K)
    zpow = one(K)
    for c in coeffs
        x += K(c) * zpow
        zpow *= z
    end
    return x // K(denom)
end

function reconstruct_cyclotomic_element_from_residues(values_by_prime::Dict{Int, Int},
                                                      ctx::CyclotomicContext;
                                                      coeff_bound::Int = 4,
                                                      denominator_bound::Int = 4)
    isempty(values_by_prime) && return nothing
    primes = sort(collect(keys(values_by_prime)))
    zetas = [find_zeta_in_Fp(ctx.N, p) for p in primes]
    degree_K = degree(field(ctx))

    best_coeffs = nothing
    best_denom = 0
    best_score = nothing
    for denom in 1:denominator_bound
        any(p -> denom % p == 0, primes) && continue
        targets = [mod(denom * values_by_prime[p], p) for p in primes]
        coeffs = _reconstruct_power_basis_coeffs_mitm(targets, zetas, primes,
                                                      degree_K, coeff_bound)
        coeffs === nothing && continue
        (l1, max_abs, _) = _coeff_score(coeffs)
        score = (l1, max_abs, denom, coeffs)
        if best_coeffs === nothing || score < best_score
            best_coeffs = coeffs
            best_denom = denom
            best_score = score
        end
    end
    best_coeffs === nothing && return nothing
    return _cyclotomic_element_from_coeffs(best_coeffs, best_denom, ctx)
end

function crt_reconstruction_is_well_determined(ctx::CyclotomicContext,
                                               modular_filters;
                                               reconstruction_bound::Int,
                                               denominator_bound::Int)
    modulus = BigInt(1)
    for filter in modular_filters
        modulus *= filter.p
    end
    degree_K = degree(field(ctx))
    search_space = BigInt(denominator_bound) * BigInt(2 * reconstruction_bound + 1)^degree_K
    return modulus >= search_space
end
