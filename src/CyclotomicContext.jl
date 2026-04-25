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

"""
    ModularData

Exact modular data over a fixed cyclotomic context.
"""
struct ModularData
    context::CyclotomicContext
    labels::Vector{Symbol}
    S::Any
    T::Any
    cond_S::Int
    cond_T::Int
    cond_F::Union{Int, Nothing}
end

field(ctx::CyclotomicContext) = ctx.field
zeta(ctx::CyclotomicContext) = ctx.zeta
conductor(ctx::CyclotomicContext) = ctx.conductor
field(data::ModularData) = field(data.context)
conductor(data::ModularData) = lcm(data.cond_S, data.cond_T)
cond_S(data::ModularData) = data.cond_S
cond_T(data::ModularData) = data.cond_T
cond_F(data::ModularData) = data.cond_F

function _require_same_context(ctx::CyclotomicContext, data_name::Symbol, required::Int)
    ctx.N == required ||
        error("$data_name modular data is defined over Q(ζ_$required); got Q(ζ_$(ctx.N))")
    return nothing
end

function _sqrt2(ctx::CyclotomicContext)
    ctx.N % 8 == 0 || error("√2 requires Q(ζ_N) with 8 | N")
    z = ctx.zeta
    return z^(ctx.N ÷ 8) + z^(7ctx.N ÷ 8)
end

function _sqrt5(ctx::CyclotomicContext)
    ctx.N % 5 == 0 || error("√5 requires Q(ζ_N) with 5 | N")
    z = ctx.zeta
    return one(ctx.field) + 2 * (z^(ctx.N ÷ 5) + z^(4ctx.N ÷ 5))
end

"""
    semion_modular_data(ctx = CyclotomicContext(8))
"""
function semion_modular_data(ctx::CyclotomicContext = CyclotomicContext(8))
    _require_same_context(ctx, :semion, 8)
    K, z = ctx.field, ctx.zeta
    inv_sqrt2 = inv(_sqrt2(ctx))
    S = inv_sqrt2 * matrix(K, 2, 2, [1, 1, 1, -1])
    T = matrix(K, 2, 2, [1, 0, 0, z^2])
    return ModularData(ctx, [:one, :s], S, T, 8, 4, nothing)
end

"""
    fibonacci_modular_data(ctx = CyclotomicContext(20))
"""
function fibonacci_modular_data(ctx::CyclotomicContext = CyclotomicContext(20))
    _require_same_context(ctx, :Fibonacci, 20)
    K, z = ctx.field, ctx.zeta
    sqrt5 = _sqrt5(ctx)
    phi = (one(K) + sqrt5) // K(2)
    invD = z + z^(-1) - z^3 - z^(-3)
    S = invD * matrix(K, 2, 2, [1, phi, phi, -1])
    T = matrix(K, 2, 2, [1, 0, 0, z^8])
    return ModularData(ctx, [:one, :tau], S, T, 20, 5, nothing)
end

"""
    ising_modular_data(ctx = CyclotomicContext(16))
"""
function ising_modular_data(ctx::CyclotomicContext = CyclotomicContext(16))
    _require_same_context(ctx, :Ising, 16)
    K, z = ctx.field, ctx.zeta
    sqrt2 = _sqrt2(ctx)
    S = (K(1) // K(2)) * matrix(K, 3, 3,
                                [1, sqrt2, 1,
                                 sqrt2, 0, -sqrt2,
                                 1, -sqrt2, 1])
    T = matrix(K, 3, 3, [1, 0, 0,
                         0, z, 0,
                         0, 0, -1])
    return ModularData(ctx, [:one, :sigma, :psi], S, T, 16, 16, nothing)
end

function modular_data(name::Symbol; conductor::Union{Int, Nothing} = nothing)
    if name == :Semion || name == :semion
        ctx = CyclotomicContext(conductor === nothing ? 8 : conductor)
        return semion_modular_data(ctx)
    elseif name == :Fibonacci || name == :fibonacci
        ctx = CyclotomicContext(conductor === nothing ? 20 : conductor)
        return fibonacci_modular_data(ctx)
    elseif name == :Ising || name == :ising
        ctx = CyclotomicContext(conductor === nothing ? 16 : conductor)
        return ising_modular_data(ctx)
    end
    error("unknown modular-data family: $name")
end

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

frobenius(data::ModularData, p::Int) = galois_action(data, mod(p, data.context.N))

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
