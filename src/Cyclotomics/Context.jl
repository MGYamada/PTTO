"""
    CyclotomicContext

Arithmetic context for computations over `Q(ζ_N)`.

The conductor `N` is a choice of ground field, not a value inferred
after constructing modular data.  Objects carrying modular data should
keep this context explicitly so Galois action and finite-field reduction
are tied to the same cyclotomic generator.
"""

using Oscar

struct CyclotomicContext
    N::Int
    field::Any
    zeta::Any
    real_subfield::Any
    conductor::Int
end

function CyclotomicContext(N::Int)
    N >= 1 || error("cyclotomic conductor must be positive, got $N")
    K, z = cyclotomic_field(N)
    return CyclotomicContext(N, K, z, nothing, N)
end
field(ctx::CyclotomicContext) = ctx.field
zeta(ctx::CyclotomicContext) = ctx.zeta
conductor(ctx::CyclotomicContext) = ctx.conductor
