"""
Galois actions on cyclotomic elements and exact modular data.

The maps are conductor-local: exponents are interpreted using the
CyclotomicContext carried by the input data.
"""

function _galois_elem(ctx::CyclotomicContext, x, a::Int)
    gcd(a, ctx.N) == 1 || error("Galois exponent $a is not a unit modulo $(ctx.N)")
    coeffs = coordinates(x)
    za = ctx.zeta^mod(a, ctx.N)
    y = zero(ctx.field)
    for k in eachindex(coeffs)
        y += ctx.field(Rational{BigInt}(coeffs[k])) * za^(k - 1)
    end
    return y
end

function galois_action(ctx::CyclotomicContext, x, a::Int)
    return _galois_elem(ctx, x, a)
end

function galois_action(ctx::CyclotomicContext, M::MatElem, a::Int)
    K = ctx.field
    out = zero_matrix(K, nrows(M), ncols(M))
    for i in 1:nrows(M), j in 1:ncols(M)
        out[i, j] = galois_action(ctx, M[i, j], a)
    end
    return out
end

function galois_action(ctx::CyclotomicContext, v::AbstractVector, a::Int)
    return [galois_action(ctx, x, a) for x in v]
end

function galois_action(data::ModularData, a::Int)
    ctx = data.context
    return ModularData(ctx, copy(data.labels),
                       galois_action(ctx, data.S, a),
                       galois_action(ctx, data.T, a),
                       data.cond_S, data.cond_T, data.cond_F)
end

galois_orbit(data::ModularData) =
    [galois_action(data, a) for a in 1:data.context.N if gcd(a, data.context.N) == 1]
