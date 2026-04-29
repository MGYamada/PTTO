"""
Exact modular-data containers and small built-in examples over Q(zeta_N).

These helpers attach S/T data to an explicit CyclotomicContext without
changing conductor choices implicitly.
"""

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

"""
    toric_code_modular_data(ctx = CyclotomicContext(2))
"""
function toric_code_modular_data(ctx::CyclotomicContext = CyclotomicContext(2))
    _require_same_context(ctx, :toric_code, 2)
    K, z = ctx.field, ctx.zeta
    S = (K(1) // K(2)) * matrix(K, 4, 4,
                                [1, 1, 1, 1,
                                 1, 1, -1, -1,
                                 1, -1, 1, -1,
                                 1, -1, -1, 1])
    T = matrix(K, 4, 4, [1, 0, 0, 0,
                         0, 1, 0, 0,
                         0, 0, 1, 0,
                         0, 0, 0, z])
    return ModularData(ctx, [:one, :e, :m, :epsilon], S, T, 1, 2, nothing)
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
    elseif name == :ToricCode || name == :toric_code || name == :toric
        ctx = CyclotomicContext(conductor === nothing ? 2 : conductor)
        return toric_code_modular_data(ctx)
    end
    error("unknown modular-data family: $name")
end
