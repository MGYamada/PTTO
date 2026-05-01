"""
Type definitions for scalar gauge actions on multiplicity-free F/R data.

The current gauge layer acts by one invertible scalar on each basis vector of
`Hom(a ⊗ b, c)`.  Vector-backed `FRData` is multiplicity-free, so the only
supported Hom-basis index is currently `μ = 1`; the explicit four-index key is
kept as the public shape for future multiplicity-aware backends.
"""

"""
    GaugeTransform

Legacy multiplicity-free channel-scalar gauge transform.  Keys are
`(a,b,c)` and represent the scalar on `Hom(a ⊗ b, c)`.
"""
struct GaugeTransform{T}
    scalars::Dict{Tuple{Int,Int,Int}, T}
    fixed_indices::Vector{Int}
    complete::Bool
end

"""
    GaugeParameters

Typed gauge parameters with an explicit Hom-basis index.  Keys are
`(a,b,c,μ)` for a basis vector of `Hom(a ⊗ b, c)`.  The legacy
`GaugeTransform` scalar dictionary is still accepted and represents the
multiplicity-free `μ = 1` case.
"""
struct GaugeParameters{T, D<:AbstractDict{Symbol}}
    values::Dict{Tuple{Int,Int,Int,Int}, T}
    metadata::D
end

GaugeParameters(values::Dict{Tuple{Int,Int,Int,Int}, T};
                metadata::Dict{Symbol, Any} = Dict{Symbol, Any}()) where {T} =
    GaugeParameters{T, typeof(metadata)}(values, metadata)

const GaugeChoice = GaugeParameters

"""
    GaugeAction

FRData-centered gauge action.  The `parameters` dictionary is keyed by
`(a,b,c,μ)` and stores the nonzero scalar acting on the chosen Hom-basis
vector.  `field` records backend metadata such as `:F_101`, a scalar type, or
`nothing`; arithmetic is performed by the scalar values themselves.
"""
struct GaugeAction{T, K, D<:AbstractDict{Symbol}}
    parameters::Dict{Tuple{Int,Int,Int,Int}, T}
    field::K
    metadata::D
end

GaugeAction(parameters::Dict{Tuple{Int,Int,Int,Int}, T};
            field = nothing,
            metadata::Dict{Symbol, Any} = Dict{Symbol, Any}()) where {T} =
    GaugeAction{T, typeof(field), typeof(metadata)}(parameters, field, metadata)

GaugeAction(parameters::Dict{Tuple{Int,Int,Int}, T};
            field = nothing,
            metadata::Dict{Symbol, Any} = Dict{Symbol, Any}()) where {T} =
    GaugeAction(Dict((a, b, c, 1) => value
                     for ((a, b, c), value) in parameters);
                field = field, metadata = metadata)

GaugeAction(parameters::GaugeParameters; field = get(parameters.metadata, :field, nothing),
            metadata::Dict{Symbol, Any} = copy(parameters.metadata)) =
    GaugeAction(copy(parameters.values); field = field, metadata = metadata)

GaugeAction(transform::GaugeTransform; field = nothing,
            metadata::Dict{Symbol, Any} = Dict{Symbol, Any}()) =
    GaugeAction(transform.scalars; field = field,
                metadata = merge(Dict{Symbol, Any}(:legacy_complete => transform.complete,
                                                   :legacy_fixed_indices => transform.fixed_indices),
                                 metadata))

"""
    GaugeFixingResult

Result of a gauge-fixing or normal-form computation.  `F` and `R` are the
transformed vector-backed symbols, `gauge` records the action used, and
`metadata[:frdata]` stores the transformed `FRData` when the input was FRData.
"""
struct GaugeFixingResult{FT,RT,GT,D<:AbstractDict{Symbol}}
    F::FT
    R::RT
    gauge::GT
    fixed_indices::Vector{Int}
    complete::Bool
    metadata::D
end

"""
    ToricGaugeNormalFormResult

Result of a Smith-normal-form guided toric gauge normalization.  `frdata` is
the normalized F/R data, `gauge` is the applied `GaugeAction`, and `split`
records the Smith invariant factors of the relevant character matrix.
"""
struct ToricGaugeNormalFormResult{FT<:FRData,GT,D<:AbstractDict{Symbol}}
    frdata::FT
    gauge::GT
    fixed_indices::Vector{Int}
    split::Any
    stabilizer_size::Any
    stacky_weight::Any
    metadata::D
end

"""
    GaugeDegreeOfFreedom

Descriptor for one scalar gauge coordinate on `Hom(a ⊗ b, c)`.
"""
struct GaugeDegreeOfFreedom
    a::Int
    b::Int
    c::Int
    μ::Int
end

"""
Abstract supertype for explicit gauge-fixing constraints.
"""
abstract type GaugeConstraint end

"""
    FixUnitConstraints()

Constraint marker for preserving the unit-channel convention.  The current
scalar action formulas already treat absent parameters as identity, so this is
recorded as part of the normal-form contract.
"""
struct FixUnitConstraints <: GaugeConstraint end

"""
    FixSelectedFSymbols(indices)

Constraint marker requesting selected F-symbol coordinates to be normalized.
The multiplicity-free normal-form solver currently recomputes the compatible
pivot plan from the input `FRData`.
"""
struct FixSelectedFSymbols <: GaugeConstraint
    indices::Vector{Int}
end

FixSelectedFSymbols(indices::AbstractVector{<:Integer}) =
    FixSelectedFSymbols(Int[i for i in indices])

"""
    FixSelectedRSymbols(indices)

Reserved constraint marker for future R-symbol normalization strategies.
"""
struct FixSelectedRSymbols <: GaugeConstraint
    indices::Vector{Int}
end

FixSelectedRSymbols(indices::AbstractVector{<:Integer}) =
    FixSelectedRSymbols(Int[i for i in indices])

"""
    NormalizationConstraint(strategy)

Records the normal-form strategy used to select a representative.
"""
struct NormalizationConstraint <: GaugeConstraint
    strategy::Symbol
end
