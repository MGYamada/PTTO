"""
Finite-field reduction methods for context-carrying exact modular data.

These methods reduce CyclotomicContext values and ModularData containers
using the cyclotomic reduction primitives.
"""

function reduce_mod_p(ctx::CyclotomicContext, x, p::Int)
    gcd(p, ctx.N) == 1 || error("prime $p must be coprime to conductor $(ctx.N)")
    zeta_Fp = find_zeta_in_Fp(ctx.N, p)
    return cyclotomic_to_Fp(x, zeta_Fp, p)
end

function reduce_mod_p(ctx::CyclotomicContext, M::MatElem, p::Int)
    return [reduce_mod_p(ctx, M[i, j], p) for i in 1:nrows(M), j in 1:ncols(M)]
end

function reduce_mod_p(data::ModularData, p::Int)
    ctx = data.context
    return (S = reduce_mod_p(ctx, data.S, p),
            T = reduce_mod_p(ctx, data.T, p),
            p = p,
            conductor = ctx.N)
end
