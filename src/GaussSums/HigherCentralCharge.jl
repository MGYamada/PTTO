"""
Higher central charges as D²-normalized twist-spectrum moments.
"""

struct HigherCentralChargeResult
    ok::Bool
    n::Int
    normalization::Symbol
    value::Any
    gauss_sum::Any
    denominator::Any
    D_squared::Any
    conductor::Union{Int, Nothing}
    status::Symbol
    message::String
end

struct HCCGeneratingFunction
    weights::Vector{Any}
    twists::Vector{Any}
    period::Union{Int, Nothing}
end

struct HCCLocalFactor
    p::Int
    coefficients::Vector{Int}
    variable::Symbol
end

function Base.show(io::IO, r::HigherCentralChargeResult)
    if r.ok
        print(io, "HigherCentralChargeResult(n=$(r.n), normalization=$(r.normalization), value=$(r.value))")
    else
        print(io, "HigherCentralChargeResult(n=$(r.n), status=$(r.status), message=$(repr(r.message)))")
    end
end

function Base.show(io::IO, gf::HCCGeneratingFunction)
    print(io, "HCCGeneratingFunction(length=$(length(gf.weights)), period=$(gf.period))")
end

function Base.show(io::IO, f::HCCLocalFactor)
    coeffs = f.coefficients
    if length(coeffs) == 2 && coeffs[1] == 1
        a = mod(-coeffs[2], f.p)
        print(io, "1 - $(a)$(f.variable) over F_$(f.p)")
    else
        print(io, "HCCLocalFactor(p=$(f.p), coefficients=$(coeffs), variable=$(f.variable))")
    end
end

Base.length(gf::HCCGeneratingFunction) = length(gf.weights)

function (gf::HCCGeneratingFunction)(n::Integer)
    isempty(gf.weights) && error("empty higher-central-charge generating function")
    K = parent(gf.weights[1])
    total = zero(K)
    nn = gf.period === nothing ? Int(n) : mod(Int(n), gf.period)
    for i in eachindex(gf.weights)
        total += gf.weights[i] * gf.twists[i]^nn
    end
    return total
end

_hcc_data(data::ModularData) = data
_hcc_data(name::Symbol) = modular_data(name)
_hcc_data(name::AbstractString) =
    modular_data(Symbol(lowercase(replace(name, "-" => "_", " " => "_"))))

function _hcc_data(data)
    input = _gauss_input(data)
    labels = [Symbol("x$i") for i in 1:_rank_of(input)]
    return ModularData(input.context, labels, input.S, input.T,
                       conductor(input.context), conductor(input.context), nothing)
end

simple_dimensions(data) = quantum_dimensions(_hcc_data(data))

function simple_twists(data)
    d = _hcc_data(data)
    return [_twist_at(d.T, i) for i in 1:length(d.labels)]
end

global_dimension_squared(data) = total_quantum_dimension_squared(_hcc_data(data))
twist_conductor(data::ModularData) = data.cond_T
twist_conductor(data) = conductor(_hcc_data(data))

_hcc_normalization_error(normalization) =
    error("unsupported higher-central-charge normalization: $(normalization); use :D2")

"""
    higher_central_charge(data, n = 1; normalization = :D2)

Compute the D²-normalized higher central charge
`sum_a d_a^2 * theta_a^n / sum_a d_a^2`.

This moment convention normalizes by the squared global dimension.  In
particular, `higher_central_charge(data, 0) == 1` for valid modular data.
Negative exponents are supported whenever the twists are invertible.
"""
function higher_central_charge(data, n::Integer;
                               normalization::Symbol = :D2,
                               method::Symbol = :modular_data,
                               p::Union{Integer, Nothing} = nothing,
                               embedding = nothing,
                               kwargs...)
    isempty(kwargs) || error("unsupported keyword arguments: $(collect(keys(kwargs)))")
    normalization == :D2 || _hcc_normalization_error(normalization)
    if method == :finite_field
        p === nothing && error("finite-field higher central charge requires keyword p")
        return higher_central_charge_modp(data, n, Int(p); embedding = embedding)
    elseif method != :modular_data
        error("unknown higher central charge method: $method; expected :modular_data or :finite_field")
    end

    md = _hcc_data(data)
    period = higher_central_charge_period(md)
    nn = period === nothing ? Int(n) : mod(Int(n), period)
    ds = simple_dimensions(md)
    twists = simple_twists(md)
    D2 = global_dimension_squared(md)
    K = md.context.field
    numerator = zero(K)
    for i in eachindex(ds)
        numerator += ds[i]^2 * twists[i]^nn
    end
    return numerator / D2
end

higher_central_charge(data; n::Integer = 1, kwargs...) =
    higher_central_charge(data, n; kwargs...)

"""
    higher_central_charges(data, ns; normalization = :D2)

Return the D²-normalized higher central charges for every integer in `ns`.
"""
function higher_central_charges(data, ns; kwargs...)
    return [higher_central_charge(data, n; kwargs...) for n in ns]
end

"""
    higher_central_charge_period(data)

Return the period used by the HCC sequence.  For exact `ModularData` this is
the declared twist conductor.
"""
higher_central_charge_period(data) = twist_conductor(data)

"""
    higher_central_charge_sequence(data)

Return one full period `[ξ_0, ξ_1, ..., ξ_{N-1}]` of the D²-normalized HCC
sequence, where `N = higher_central_charge_period(data)`.
"""
function higher_central_charge_sequence(data; kwargs...)
    period = higher_central_charge_period(data)
    period === nothing && error("higher central charge period is unknown")
    return higher_central_charges(data, 0:(period - 1); kwargs...)
end

"""
    higher_central_charge_generating_function(data)

Return a structured generating function
`sum_a weights[a] / (1 - twists[a] * u)`, where
`weights[a] = d_a^2 / D^2`.
"""
function higher_central_charge_generating_function(data)
    md = _hcc_data(data)
    ds = simple_dimensions(md)
    twists = simple_twists(md)
    D2 = global_dimension_squared(md)
    weights = Any[d^2 / D2 for d in ds]
    return HCCGeneratingFunction(weights, Any[twists...],
                                 higher_central_charge_period(md))
end

function _higher_failure(category_or_modular_data, n::Integer, normalization::Symbol,
                         status::Symbol, message::String)
    input = try
        _gauss_input(category_or_modular_data)
    catch
        nothing
    end
    conductor_value = input === nothing ? nothing : conductor(input.context)
    return HigherCentralChargeResult(false, Int(n), normalization, nothing, nothing, nothing,
                                     nothing, conductor_value, status, message)
end

"""
    higher_central_charge_result(data; n = 1, normalization = :galois)

Legacy structured Gauss-sum result used by compatibility code.  New HCC
moment computations should call `higher_central_charge(data, n)`, which uses
D² normalization and returns the exact value directly.
"""
function higher_central_charge_result(category_or_modular_data;
                                      n::Integer = 1,
                                      normalization::Symbol = :galois)
    input = try
        _gauss_input(category_or_modular_data)
    catch err
        return _higher_failure(category_or_modular_data, n, normalization,
                               :missing_data, sprint(showerror, err))
    end

    try
        _check_gauss_input(input)
        τ = gauss_sum_plus(category_or_modular_data; n = n)
        D = inv(input.S[1, 1])
        D2 = total_quantum_dimension_squared(category_or_modular_data)
        denom = if normalization == :galois
            gcd(Int(n), input.context.N) == 1 ||
                return HigherCentralChargeResult(false, Int(n), normalization, nothing, τ,
                                                 nothing, D2, conductor(input.context),
                                                 :not_galois,
                                                 "n=$(Int(n)) is not a unit modulo conductor $(input.context.N)")
            galois_action(input.context, D, Int(n))
        elseif normalization == :D
            D
        elseif normalization == :D2
            D2
        elseif normalization == :raw
            one(input.context.field)
        else
            return HigherCentralChargeResult(false, Int(n), normalization, nothing, τ,
                                             nothing, D2, conductor(input.context),
                                             :unknown_normalization,
                                             "expected :galois, :D, :D2, or :raw")
        end
        value = τ / denom
        return HigherCentralChargeResult(true, Int(n), normalization, value, τ, denom,
                                         D2, conductor(input.context), :ok, "ok")
    catch err
        return HigherCentralChargeResult(false, Int(n), normalization, nothing, nothing,
                                         nothing, nothing, conductor(input.context),
                                         :error, sprint(showerror, err))
    end
end

"""
    central_charge(data)

Return the ordinary exact central-charge moment using the same D²-normalized
convention as `higher_central_charge(data, 1)`.

This is the moment `p₊ / D²`, not the unit `p₊ / |p₊|`.  Its argument
encodes the central charge modulo 8, while its magnitude is generally `1 / D`.
For the unit form, use `value / abs(value)` numerically or call
`higher_central_charge_result(data; normalization = :galois)`.
"""
function central_charge(category_or_modular_data)
    return higher_central_charge(category_or_modular_data, 1)
end

function _qq_to_rational(c)
    return BigInt(numerator(c)) // BigInt(denominator(c))
end

function _coeff_vector_rational(x, len::Int)
    coeffs = Oscar.coefficients(x)
    v = fill(BigInt(0) // BigInt(1), len)
    for i in 1:min(length(coeffs), len)
        v[i] = _qq_to_rational(coeffs[i])
    end
    if length(coeffs) > len
        for i in (len + 1):length(coeffs)
            iszero(coeffs[i]) || error("coefficient vector is longer than the cyclotomic degree")
        end
    end
    return v
end

function _rational_linear_solve(columns::Vector{Vector{Rational{BigInt}}},
                                b::Vector{Rational{BigInt}})
    m = length(b)
    n = length(columns)
    A = Matrix{Rational{BigInt}}(undef, m, n + 1)
    for i in 1:m
        for j in 1:n
            A[i, j] = columns[j][i]
        end
        A[i, n + 1] = b[i]
    end

    row = 1
    pivots = Int[]
    for col in 1:n
        pivot = 0
        for i in row:m
            if !iszero(A[i, col])
                pivot = i
                break
            end
        end
        pivot == 0 && continue
        if pivot != row
            for j in 1:(n + 1)
                A[row, j], A[pivot, j] = A[pivot, j], A[row, j]
            end
        end
        scale = A[row, col]
        for j in 1:(n + 1)
            A[row, j] /= scale
        end
        for i in 1:m
            i == row && continue
            factor = A[i, col]
            iszero(factor) && continue
            for j in 1:(n + 1)
                A[i, j] -= factor * A[row, j]
            end
        end
        push!(pivots, col)
        row += 1
        row > m && break
    end

    for i in 1:m
        all(iszero(A[i, j]) for j in 1:n) && !iszero(A[i, n + 1]) &&
            error("element is not in the requested cyclotomic subfield")
    end
    length(pivots) == n ||
        error("cyclotomic subfield basis is not independent in the ambient field")

    x = fill(BigInt(0) // BigInt(1), n)
    for (i, col) in enumerate(pivots)
        x[col] = A[i, n + 1]
    end
    return x
end

function _modp_from_rational(q::Rational{BigInt}, p::Int)
    den = Int(mod(denominator(q), p))
    den != 0 || error("bad prime $p: denominator $(denominator(q)) is divisible by p")
    num = mod(numerator(q), p)
    return Int(mod(num * invmod(den, p), p))
end

function _eval_rational_poly_modp(coeffs::AbstractVector{Rational{BigInt}},
                                  root::Int, p::Int)
    value = 0
    power = 1
    for c in coeffs
        value = mod(value + _modp_from_rational(c, p) * power, p)
        power = mod(power * root, p)
    end
    return value
end

function _is_primitive_root_mod(root::Integer, order::Integer, p::Int)
    order == 1 && return mod(root, p) == 1
    powermod(mod(root, p), Int(order), p) == 1 || return false
    for (q, _) in Primes.factor(Int(order))
        powermod(mod(root, p), Int(order) ÷ q, p) == 1 && return false
    end
    return true
end

function _parse_embedding(data, p::Int, embedding)
    period = higher_central_charge_period(data)
    embedding === nothing && return (period, find_zeta_in_Fp(period, p))

    if embedding isa Pair
        order = embedding.first
        root = embedding.second
    elseif embedding isa Integer
        order = period
        root = embedding
    elseif embedding isa NamedTuple && haskey(embedding, :order) && haskey(embedding, :root)
        order = embedding.order
        root = embedding.root
    else
        error("embedding must be `nothing`, an integer root for ζ_$period, `order => root`, or `(order = ..., root = ...)`")
    end
    order isa Integer || error("embedding root order must be an integer")
    root isa Integer || error("embedding root must be an integer residue modulo p")
    order = Int(order)
    root = mod(Int(root), p)
    _is_primitive_root_mod(root, order, p) ||
        error("embedding sends ζ_$order to $root, which is not a primitive $order-th root in F_$p")
    return (order, root)
end

function _reduce_subfield_element_modp(ctx::CyclotomicContext, x, order::Int,
                                       root::Int, p::Int)
    ctx.N % order == 0 ||
        error("cannot use ζ_$order embedding for data over Q(ζ_$(ctx.N)); order must divide $(ctx.N)")
    degN = Int(Primes.totient(ctx.N))
    degm = Int(Primes.totient(order))
    η = ctx.zeta^(ctx.N ÷ order)
    columns = [_coeff_vector_rational(η^j, degN) for j in 0:(degm - 1)]
    b = _coeff_vector_rational(x, degN)
    coeffs = try
        _rational_linear_solve(columns, b)
    catch err
        error("embedding ζ_$order is not enough to reduce this exact value; an extension field or a larger embedding is needed: $(sprint(showerror, err))")
    end
    return _eval_rational_poly_modp(coeffs, root, p)
end

function _hcc_reduce_modp(ctx::CyclotomicContext, x, order::Int, root::Int, p::Int)
    return _reduce_subfield_element_modp(ctx, x, order, root, p)
end

function _hcc_fp_pow_unit(a::Integer, n::Integer, p::Int)
    a = mod(a, p)
    n == 0 && return 1
    n > 0 && return powermod(a, n, p)
    a == 0 && error("cannot raise zero to a negative power in F_$p")
    return powermod(invmod(a, p), -n, p)
end

"""
    higher_central_charge_modp(data, n, p; embedding = nothing)

Reduce the D²-normalized HCC value modulo a good prime `p`.

By default ACMG chooses a primitive root for the twist conductor in `F_p`.
Pass `embedding = order => root` to choose an explicit cyclotomic embedding,
for example `embedding = 5 => 3` for the Fibonacci `p = 11` regression.
Only prime-field reductions are currently implemented; primes needing
extension fields throw a clear error.
"""
function higher_central_charge_modp(data, n::Integer, p::Integer; embedding = nothing)
    md = _hcc_data(data)
    pp = Int(p)
    isprime(pp) || error("higher_central_charge_modp requires a prime p, got $p")
    period = higher_central_charge_period(md)
    gcd(pp, period) == 1 ||
        error("bad prime $pp for HCC reduction: p divides twist conductor $period")
    order, root = _parse_embedding(md, pp, embedding)
    ds = simple_dimensions(md)
    twists = simple_twists(md)

    weights_num = Int[]
    D2 = 0
    for d in ds
        d2 = _hcc_reduce_modp(md.context, d^2, order, root, pp)
        push!(weights_num, d2)
        D2 = mod(D2 + d2, pp)
    end
    D2 != 0 || error("bad prime $pp for HCC reduction: D^2 is zero in F_$pp")

    nn = mod(Int(n), period)
    numerator = 0
    for i in eachindex(weights_num)
        θ = _hcc_reduce_modp(md.context, twists[i], order, root, pp)
        numerator = mod(numerator + weights_num[i] * _hcc_fp_pow_unit(θ, nn, pp), pp)
    end
    return mod(numerator * invmod(D2, pp), pp)
end

"""
    higher_central_charge_sequence_modp(data, p; embedding = nothing)

Return one full HCC period reduced modulo `p`.
"""
function higher_central_charge_sequence_modp(data, p::Integer; embedding = nothing)
    period = higher_central_charge_period(data)
    return [higher_central_charge_modp(data, n, p; embedding = embedding)
            for n in 0:(period - 1)]
end

"""
    hcc_local_factor(data, n, p; embedding = nothing)

Return the local factor `1 - ξ_{n,p} T` over `F_p` for prime-field HCC
reductions.  Extension-field norm factors are future work.
"""
function hcc_local_factor(data, n::Integer, p::Integer; embedding = nothing)
    value = higher_central_charge_modp(data, n, p; embedding = embedding)
    return HCCLocalFactor(Int(p), [1, mod(-value, Int(p))], :T)
end

"""
    hcc_local_factors(data, p; embedding = nothing)

Return one full period of prime-field HCC local factors.
"""
function hcc_local_factors(data, p::Integer; embedding = nothing)
    period = higher_central_charge_period(data)
    return [hcc_local_factor(data, n, p; embedding = embedding)
            for n in 0:(period - 1)]
end
