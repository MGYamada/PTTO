"""
General gauge-group infrastructure for fusion spaces.

For a fusion tensor `N_ab^c`, the general gauge group is the product

    prod_{a,b,c : N_ab^c > 0} GL(N_ab^c).

The existing multiplicity-free toric gauge group is the special case where
each positive factor has dimension one.
"""

struct GeneralGaugeValidationError <: Exception
    message::String
end

Base.showerror(io::IO, err::GeneralGaugeValidationError) = print(io, err.message)

struct GeneralGaugeActionNotImplementedError <: Exception
    message::String
end

Base.showerror(io::IO, err::GeneralGaugeActionNotImplementedError) = print(io, err.message)

"""
    FusionSpaceIndex(a, b, c)

Index for the fusion space `Hom(a ⊗ b, c)`.
"""
struct FusionSpaceIndex
    a::Int
    b::Int
    c::Int
end

FusionSpaceIndex(ch::Tuple{Int,Int,Int}) = FusionSpaceIndex(ch[1], ch[2], ch[3])
Base.Tuple(idx::FusionSpaceIndex) = (idx.a, idx.b, idx.c)
Base.isless(x::FusionSpaceIndex, y::FusionSpaceIndex) = Tuple(x) < Tuple(y)
Base.show(io::IO, idx::FusionSpaceIndex) = print(io, "(", idx.a, ",", idx.b, ",", idx.c, ")")

"""
    GaugeFactor(index, dim)

One factor `GL(dim)` of the general gauge group attached to a fusion space.
"""
struct GaugeFactor
    index::FusionSpaceIndex
    dim::Int
    function GaugeFactor(index::FusionSpaceIndex, dim::Integer)
        d = Int(dim)
        d > 0 || throw(GeneralGaugeValidationError("gauge factor $index has non-positive dimension $d"))
        return new(index, d)
    end
end

"""
    GeneralGaugeData(fusion, factors, field, multiplicity_free)

Structured representation of `prod GL(N_ab^c)`.
"""
struct GeneralGaugeData{FT,F}
    fusion::FT
    factors::Vector{GaugeFactor}
    field::F
    multiplicity_free::Bool
    metadata::Dict{Symbol, Any}
end

function Base.show(io::IO, gauge::GeneralGaugeData)
    kind = is_toric(gauge) ? "toric" : "general"
    print(io, "GeneralGaugeData($(kind), rank=$(gauge_rank(gauge)), factors=$(length(gauge.factors)))")
end

struct GaugeTransformation{G<:GeneralGaugeData,T,D<:AbstractDict}
    gauge::G
    matrices::Dict{FusionSpaceIndex, Matrix{T}}
    metadata::D
end

function GaugeTransformation(gauge::GeneralGaugeData,
                             matrices::Dict{FusionSpaceIndex, Matrix{T}};
                             metadata::Dict{Symbol, Any} = Dict{Symbol, Any}()) where {T}
    return GaugeTransformation{typeof(gauge), T, typeof(metadata)}(gauge, matrices, metadata)
end

function GaugeTransformation(gauge::GeneralGaugeData,
                             matrices::Dict{Tuple{Int,Int,Int}, Matrix{T}};
                             metadata::Dict{Symbol, Any} = Dict{Symbol, Any}()) where {T}
    converted = Dict(FusionSpaceIndex(ch) => M for (ch, M) in matrices)
    return GaugeTransformation(gauge, converted; metadata = metadata)
end

_raw_fusion_tensor(fr::FusionRule) = fr.N
_raw_fusion_tensor(Nijk::Array{Int,3}) = Nijk

function _raw_fusion_tensor(data::FRData)
    return fusion_rule(data).N
end

function _validate_fusion_tensor_shape(Nijk::Array{Int,3})
    r = size(Nijk, 1)
    size(Nijk) == (r, r, r) ||
        throw(GeneralGaugeValidationError("fusion multiplicity tensor must be cubic; got size $(size(Nijk))"))
    all(Nijk .>= 0) ||
        throw(GeneralGaugeValidationError("fusion multiplicities must be non-negative"))
    return r
end

gauge_rank(gauge::GeneralGaugeData) = size(_raw_fusion_tensor(gauge.fusion), 1)

function _gauge_factor_map(gauge::GeneralGaugeData)
    return Dict(f.index => f for f in gauge.factors)
end

"""
    general_gauge_data(fusion; field=nothing, conventions=:tensorcategories)

Construct the gauge-group data `prod_{a,b,c} GL(N_ab^c)`.  Zero-dimensional
fusion spaces are omitted.
"""
function general_gauge_data(fusion; field = nothing, conventions = :tensorcategories)
    Nijk = _raw_fusion_tensor(fusion)
    r = _validate_fusion_tensor_shape(Nijk)
    factors = GaugeFactor[]
    multiplicity_free = true
    for a in 1:r, b in 1:r, c in 1:r
        n = Nijk[a, b, c]
        n == 0 && continue
        n == 1 || (multiplicity_free = false)
        push!(factors, GaugeFactor(FusionSpaceIndex(a, b, c), n))
    end
    metadata = Dict{Symbol, Any}(
        :gauge_group_kind => multiplicity_free ? :toric_rank_one_gl : :general_linear_product,
        :conventions => conventions,
        :formula => "prod_{a,b,c} GL(N_ab^c)",
    )
    return GeneralGaugeData{typeof(fusion), typeof(field)}(fusion, factors, field,
                                                           multiplicity_free, metadata)
end

is_multiplicity_free(gauge::GeneralGaugeData) = gauge.multiplicity_free
is_toric(gauge::GeneralGaugeData) = all(f -> f.dim == 1, gauge.factors)

"""
    gauge_group_dimension(gauge)

Return `sum_{a,b,c} (N_ab^c)^2`, the dimension of the algebraic gauge group.
"""
gauge_group_dimension(gauge::GeneralGaugeData) = sum(f.dim^2 for f in gauge.factors)

function _finite_field_prime(field)
    field isa Integer && return Int(field)
    if field isa Symbol
        m = match(r"^F_(\d+)$", String(field))
        m !== nothing && return parse(Int, m.captures[1])
    end
    return nothing
end

function _gl_order(n::Int, p::Int)
    order = BigInt(1)
    pn = BigInt(p)^n
    for i in 0:(n - 1)
        order *= pn - BigInt(p)^i
    end
    return order
end

function finite_gauge_group_order(gauge::GeneralGaugeData)
    p = _finite_field_prime(gauge.field)
    p === nothing &&
        throw(GeneralGaugeValidationError("finite gauge-group order requires a finite field marker such as :F_101"))
    isprime(p) ||
        throw(GeneralGaugeValidationError("finite gauge-group order requires prime p; got $p"))
    order = BigInt(1)
    for factor in gauge.factors
        order *= _gl_order(factor.dim, p)
    end
    return order
end

function _zero_one_for_field(field)
    p = _finite_field_prime(field)
    p !== nothing && return (zero = FpElem(0, p), one = FpElem(1, p))
    if field !== nothing && !(field isa Symbol)
        try
            return (zero = field(0), one = field(1))
        catch
        end
    end
    return (zero = 0, one = 1)
end

function _identity_matrix_for_factor(factor::GaugeFactor, field)
    elems = _zero_one_for_field(field)
    M = fill(elems.zero, factor.dim, factor.dim)
    for i in 1:factor.dim
        M[i, i] = elems.one
    end
    return M
end

"""
    identity_gauge_transformation(gauge)

Return the identity basis change in every positive-dimensional fusion space.
"""
function identity_gauge_transformation(gauge::GeneralGaugeData)
    matrices = Dict{FusionSpaceIndex, Matrix{typeof(_zero_one_for_field(gauge.field).one)}}()
    for factor in gauge.factors
        matrices[factor.index] = _identity_matrix_for_factor(factor, gauge.field)
    end
    return GaugeTransformation(gauge, matrices; metadata = Dict{Symbol, Any}(:kind => :identity))
end

function _matrix_shape(M)
    try
        return size(M)
    catch
        try
            return (nrows(M), ncols(M))
        catch
            throw(GeneralGaugeValidationError("gauge matrix of type $(typeof(M)) does not expose size or nrows/ncols"))
        end
    end
end

function _matrix_det(M::AbstractMatrix{FpElem})
    m, n = size(M)
    m == n || throw(GeneralGaugeValidationError("determinant requires a square matrix; got size $(size(M))"))
    m == 0 && return FpElem(1, 2)
    p = M[1, 1].p
    all(x -> x.p == p, M) ||
        throw(GeneralGaugeValidationError("finite-field gauge matrix mixes different prime fields"))
    A = copy(M)
    detv = FpElem(1, p)
    for col in 1:n
        pivot = findfirst(row -> !iszero(A[row, col]), col:n)
        pivot === nothing && return FpElem(0, p)
        prow = col + pivot - 1
        if prow != col
            A[col, :], A[prow, :] = copy(A[prow, :]), copy(A[col, :])
            detv = -detv
        end
        pv = A[col, col]
        detv *= pv
        invpv = inv(pv)
        for row in (col + 1):n
            iszero(A[row, col]) && continue
            factor = A[row, col] * invpv
            for j in col:n
                A[row, j] -= factor * A[col, j]
            end
        end
    end
    return detv
end

function _matrix_det(M)
    try
        return det(M)
    catch err
        throw(GeneralGaugeValidationError("could not compute determinant for gauge matrix of type $(typeof(M)): $(sprint(showerror, err))"))
    end
end

function _is_nonzero_scalar(x)
    try
        return !iszero(x)
    catch err
        throw(GeneralGaugeValidationError("could not test determinant for zero: $(sprint(showerror, err))"))
    end
end

function is_invertible_matrix_over_field(M, field = nothing)
    d = _matrix_det(M)
    return _is_nonzero_scalar(d)
end

"""
    validate_gauge_transformation(gauge, transformation)

Check channel coverage, matrix sizes, and invertibility for a gauge
transformation.
"""
function validate_gauge_transformation(gauge::GeneralGaugeData,
                                       transformation::GaugeTransformation)
    transformation.gauge === gauge ||
        throw(GeneralGaugeValidationError("gauge transformation belongs to a different GeneralGaugeData object"))
    factor_map = _gauge_factor_map(gauge)
    required = Set(keys(factor_map))
    supplied = Set(keys(transformation.matrices))
    missing = setdiff(required, supplied)
    !isempty(missing) &&
        throw(GeneralGaugeValidationError("missing gauge matrix for fusion channel $(first(sort(collect(missing))))"))
    extra = setdiff(supplied, required)
    !isempty(extra) &&
        throw(GeneralGaugeValidationError("extra gauge matrix for invalid or zero fusion channel $(first(sort(collect(extra))))"))
    for factor in gauge.factors
        M = transformation.matrices[factor.index]
        _matrix_shape(M) == (factor.dim, factor.dim) ||
            throw(GeneralGaugeValidationError("gauge matrix for fusion channel $(factor.index) has size $(_matrix_shape(M)); expected ($(factor.dim), $(factor.dim))"))
        is_invertible_matrix_over_field(M, gauge.field) ||
            throw(GeneralGaugeValidationError("gauge matrix for fusion channel $(factor.index) is singular"))
    end
    return true
end

validate_gauge_transformation(transformation::GaugeTransformation) =
    validate_gauge_transformation(transformation.gauge, transformation)

function _matrix_is_identity(M)
    m, n = _matrix_shape(M)
    m == n || return false
    for i in 1:m, j in 1:n
        expected = i == j ? one(M[i, j]) : zero(M[i, j])
        M[i, j] == expected || return false
    end
    return true
end

function _is_identity_transformation(transformation::GaugeTransformation)
    return all(_matrix_is_identity(M) for M in values(transformation.matrices))
end

function _toric_gauge_action_from_transformation(transformation::GaugeTransformation)
    is_toric(transformation.gauge) ||
        throw(GeneralGaugeActionNotImplementedError("toric scalar conversion requires all gauge factors to be GL(1)"))
    params = Dict{Tuple{Int,Int,Int}, Any}()
    for (idx, M) in transformation.matrices
        _matrix_shape(M) == (1, 1) ||
            throw(GeneralGaugeValidationError("toric gauge matrix for fusion channel $idx must be 1 x 1"))
        params[Tuple(idx)] = M[1, 1]
    end
    return GaugeAction(params; field = transformation.gauge.field,
                       metadata = Dict{Symbol, Any}(:source => :GeneralGaugeTransformation))
end

function _has_nonabelian_factor(gauge::GeneralGaugeData)
    return any(f -> f.dim > 1, gauge.factors)
end

function _nonabelian_action_error()
    throw(GeneralGaugeActionNotImplementedError(
        "general gauge action is mathematically defined but not yet implemented for nonabelian GL(n) gauge factors with n > 1"))
end

function _fusion_from_gauge(gauge::GeneralGaugeData)
    gauge.fusion isa FusionRule && return gauge.fusion
    gauge.fusion isa Array{Int,3} && return gauge.fusion
    return fusion_rule(gauge.fusion)
end

"""
    apply_gauge_to_F(Fdata, transformation; conventions=nothing)

Apply a general gauge transformation to F-symbol data.  v0.9.3 supports the
identity action in all cases and the toric `GL(1)` scalar action in the
multiplicity-free case.
"""
function apply_gauge_to_F(Fdata, transformation::GaugeTransformation; conventions = nothing)
    validate_gauge_transformation(transformation)
    _is_identity_transformation(transformation) && return copy(Fdata)
    _has_nonabelian_factor(transformation.gauge) && _nonabelian_action_error()
    one_value = isempty(Fdata) ? _zero_one_for_field(transformation.gauge.field).one : one(Fdata[1])
    action = _toric_gauge_action_from_transformation(transformation)
    return _transform_F_values(Fdata, _fusion_from_gauge(transformation.gauge), action, one_value)
end

function apply_gauge_to_R(Rdata, transformation::GaugeTransformation; conventions = nothing)
    validate_gauge_transformation(transformation)
    _is_identity_transformation(transformation) && return copy(Rdata)
    _has_nonabelian_factor(transformation.gauge) && _nonabelian_action_error()
    one_value = isempty(Rdata) ? _zero_one_for_field(transformation.gauge.field).one : one(Rdata[1])
    action = _toric_gauge_action_from_transformation(transformation)
    return _transform_R_values(Rdata, _fusion_from_gauge(transformation.gauge), action, one_value)
end

function apply_gauge_to_FR(data::FRData, transformation::GaugeTransformation; conventions = nothing)
    validate_gauge_transformation(transformation)
    _is_identity_transformation(transformation) && return data
    _has_nonabelian_factor(transformation.gauge) && _nonabelian_action_error()
    return apply_gauge(data, _toric_gauge_action_from_transformation(transformation))
end

function apply_gauge_to_FR(frdata, transformation::GaugeTransformation; conventions = nothing)
    validate_gauge_transformation(transformation)
    _is_identity_transformation(transformation) && return frdata
    _has_nonabelian_factor(transformation.gauge) && _nonabelian_action_error()
    return (F = apply_gauge_to_F(frdata.F, transformation; conventions = conventions),
            R = apply_gauge_to_R(frdata.R, transformation; conventions = conventions))
end

struct GaugeStabilizerNotComputed
    reason::Symbol
    message::String
    gauge::GeneralGaugeData
end

function gauge_stabilizer(frdata, gauge::GeneralGaugeData; kwargs...)
    return GaugeStabilizerNotComputed(:not_computed,
                                      "general gauge stabilizer computation is not implemented in v0.9.3",
                                      gauge)
end

function gauge_orbit_dimension(frdata, gauge::GeneralGaugeData; kwargs...)
    return (upper_bound = gauge_group_dimension(gauge),
            stabilizer_dimension = nothing,
            exact = false)
end
