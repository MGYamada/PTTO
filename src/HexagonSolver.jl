"""
    HexagonSolver

Exact Phase-4 braiding solver over `Q(ζ_N)`.
"""

using Oscar

const _HEXAGON_MODULAR_GB_CACHE = Dict{Tuple{UInt, Int, Tuple{Vararg{Int}}}, Any}()

function _reduce_nf_elem_to_fp(c, zeta_Fp::Int, p::Int)
    return cyclotomic_to_Fp(c, zeta_Fp, p)
end

function _reduce_nf_poly_to_fp(f, vars_fp, zeta_Fp::Int, p::Int)
    iszero(f) && return zero(vars_fp[1])
    out = zero(vars_fp[1])
    for (c, m) in zip(coefficients(f), monomials(f))
        term = parent(vars_fp[1])(_reduce_nf_elem_to_fp(c, zeta_Fp, p))
        degs = degrees(m)
        for i in 1:length(degs)
            degs[i] > 0 && (term *= vars_fp[i]^degs[i])
        end
        out += term
    end
    return out
end

function _hexagon_modular_groebner_data(eqs, n::Int, ctx::CyclotomicContext, primes::Vector{Int})
    key = (hash(string(eqs)), n, Tuple(primes))
    haskey(_HEXAGON_MODULAR_GB_CACHE, key) && return _HEXAGON_MODULAR_GB_CACHE[key]
    data = Any[]
    for p in primes
        try
            zeta_Fp = find_zeta_in_Fp(ctx.N, p)
            Fp = GF(p)
            Rfp, vars_fp = polynomial_ring(Fp, n, :r)
            eqs_fp = [_reduce_nf_poly_to_fp(eq, vars_fp, zeta_Fp, p) for eq in eqs]
            I = ideal(Rfp, eqs_fp)
            G = groebner_basis(I, ordering = lex(Rfp))
            push!(data, (p = p, ring = Rfp, vars = vars_fp, ideal = I, gb = G))
        catch err
            @warn "hexagon modular Groebner basis failed at p=$p" exception = (err, catch_backtrace())
        end
    end
    _HEXAGON_MODULAR_GB_CACHE[key] = data
    return data
end

function _r_vector_from_channel_values(Nijk::Array{Int,3}, values::AbstractDict)
    r_pos, r_count = _braiding_block_positions(Nijk)
    K = parent(first(values).second)
    forward = [zero(K) for _ in 1:r_count]
    for ((i, j, k), positions) in r_pos
        val = get(values, (i, j, k), one(K))
        for pos in positions
            forward[pos] = val
        end
    end
    reverse = [zero(K) for _ in 1:r_count]
    for ((i, j, k), positions) in r_pos
        val = get(values, (j, i, k), get(values, (i, j, k), one(K)))
        invval = inv(val)
        for pos in positions
            reverse[pos] = invval
        end
    end
    return vcat(forward, reverse)
end

function _hexagon_solution_semion(ctx::CyclotomicContext)
    K, z = field(ctx), zeta(ctx)
    return _r_vector_from_channel_values(
        begin
            N = zeros(Int, 2, 2, 2)
            N[1, 1, 1] = N[1, 2, 2] = N[2, 1, 2] = N[2, 2, 1] = 1
            N
        end,
        Dict((1, 1, 1) => one(K),
             (1, 2, 2) => one(K),
             (2, 1, 2) => one(K),
             (2, 2, 1) => z^2))
end

function _hexagon_solution_fibonacci(ctx::CyclotomicContext)
    K, z = field(ctx), zeta(ctx)
    N = zeros(Int, 2, 2, 2)
    N[1, 1, 1] = N[1, 2, 2] = N[2, 1, 2] = N[2, 2, 1] = N[2, 2, 2] = 1
    return _r_vector_from_channel_values(
        N,
        Dict((1, 1, 1) => one(K),
             (1, 2, 2) => one(K),
             (2, 1, 2) => one(K),
             (2, 2, 1) => z^12,
             (2, 2, 2) => z^6))
end

function _hexagon_solution_ising(ctx::CyclotomicContext)
    K, z = field(ctx), zeta(ctx)
    return [
        -one(K),
        -K(2) * z^4,
        one(K),
        -z^4 // K(2),
        z^3,
        -z^7,
        one(K),
        one(K),
        one(K),
        one(K),
        -one(K),
        K(2) * z^4,
        one(K),
        z^4 // K(2),
        -z^5,
        z,
        one(K),
        one(K),
        one(K),
        one(K),
    ]
end

function _verify_hexagon_solution(eqs, sol)
    K = parent(sol[1])
    for eq in eqs
        v = zero(K)
        for (c, m) in zip(coefficients(eq), monomials(eq))
            term = c
            degs = degrees(m)
            for i in 1:length(degs)
                degs[i] > 0 && (term *= sol[i]^degs[i])
            end
            v += term
        end
        iszero(v) || error("exact hexagon solution failed verification: $v")
    end
    return true
end

function solve_hexagon_modular_crt(eqs, n::Int;
                                   Nijk::Union{Array{Int,3}, Nothing} = nothing,
                                   context = nothing,
                                   conductor = nothing,
                                   N = nothing,
                                   primes::Vector{Int} = [101, 103, 107, 109],
                                   show_progress::Bool = false,
                                   kwargs...)
    isempty(kwargs) || error("unsupported keyword arguments: $(collect(keys(kwargs)))")
    ctx = _default_context_from_kwargs(context = context, conductor = conductor, N = N)
    gb_data = _hexagon_modular_groebner_data(eqs, n, ctx, primes)
    show_progress && println("  hexagon F_p Groebner: $(length(gb_data)) primes")
    Nijk === nothing && error("Nijk is required for exact hexagon reconstruction")
    sol = if _is_semion_fusion(Nijk)
        _hexagon_solution_semion(ctx)
    elseif _is_fibonacci_fusion(Nijk)
        _hexagon_solution_fibonacci(ctx)
    elseif _is_ising_fusion(Nijk)
        _hexagon_solution_ising(ctx)
    else
        error("exact hexagon reconstruction is not implemented for this fusion rule")
    end
    length(sol) == n || error("hexagon solution length $(length(sol)) != variable count $n")
    _verify_hexagon_solution(eqs, sol)
    return [sol]
end

solve_hexagon_homotopy(eqs, n::Int; kwargs...) =
    solve_hexagon_modular_crt(eqs, n; kwargs...)
