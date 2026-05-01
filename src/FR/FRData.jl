"""
Core F/R data container and accessors.

`FRData` owns TensorCategories-style F/R storage and the accessor boundary for
fusion rules, object labels, Hom-basis indices, F-symbols, and R-symbols.
Consumers such as Gauge and BraidRepresentations should use these accessors
instead of assuming field names, dictionary keys, or raw vector layout.
"""

struct FRData{T, D<:AbstractDict{Symbol}}
    rules::FusionRule
    F_values::Vector{T}
    R_values::Vector{T}
    R_inverse_values::Vector{T}
    metadata::D
end

fusion_rule(data::FRData) = data.rules
F_values(data::FRData) = data.F_values
R_values(data::FRData) = data.R_values
R_inverse_values(data::FRData) = data.R_inverse_values
fr_metadata(data::FRData) = data.metadata
fr_scalar_type(::FRData{T}) where {T} = T

"""
    simples(frdata)

Return the simple-object labels carried by `frdata`.  If no explicit labels
are recorded, ACMG's internal `1:rank` labels are returned.
"""
function simples(data::FRData)
    labels = get(data.metadata, :simple_objects, nothing)
    labels === nothing && return simple_objects(data.rules)
    length(labels) == data.rules.rank ||
        error("FRData object label metadata has length $(length(labels)); expected $(data.rules.rank)")
    return collect(labels)
end

function _frdata_object_index(data::FRData, x)
    if x isa Integer
        idx = Int(x)
    else
        labels = simples(data)
        found = findfirst(==(x), labels)
        if found === nothing && x isa AbstractString
            parsed = tryparse(Int, String(x))
            found = parsed === nothing ? findfirst(==(Symbol(x)), labels) : parsed
        end
        found === nothing &&
            error("unsupported simple-object label $x; pass an integer label in 1:$(data.rules.rank)")
        idx = Int(found)
    end
    _check_object(data.rules, idx)
    return idx
end

_frdata_object_indices(data::FRData, xs...) = map(x -> _frdata_object_index(data, x), xs)

simple_objects(data::FRData) = simples(data)

"""
    fusion_coeff(frdata, a, b, c)

Return `N_ab^c`, accepting either internal integer labels or labels recorded
in `frdata.metadata[:simple_objects]`.
"""
function fusion_coeff(data::FRData, a, b, c)
    ia, ib, ic = _frdata_object_indices(data, a, b, c)
    return data.rules.N[ia, ib, ic]
end

function fusion_coeff(data::FRData, a::Int, b::Int, c::Int)
    _check_object(data.rules, a)
    _check_object(data.rules, b)
    _check_object(data.rules, c)
    return data.rules.N[a, b, c]
end

fusion_coeff(rules::Union{FusionRule, Array{Int, 3}}, a::Int, b::Int, c::Int) =
    _fusion_rule(rules).N[a, b, c]

function fusion_channels(data::FRData, a, b)
    ia, ib = _frdata_object_indices(data, a, b)
    labels = simples(data)
    return [labels[c] for c in 1:data.rules.rank if data.rules.N[ia, ib, c] != 0]
end

function fusion_channels(data::FRData, a::Int, b::Int)
    _check_object(data.rules, a)
    _check_object(data.rules, b)
    labels = simples(data)
    return [labels[c] for c in 1:data.rules.rank if data.rules.N[a, b, c] != 0]
end

"""
    hom_basis(frdata, a, b, c)

Return the explicit Hom-basis indices for `Hom(a ⊗ b, c)`.  ACMG uses
`1:N_ab^c`; multiplicity-free channels therefore have basis `[1]`.
"""
hom_basis(data::FRData, a, b, c) = collect(1:fusion_coeff(data, a, b, c))

gauge_basis_indices(data::FRData, a, b, c) = hom_basis(data, a, b, c)

function _frdata_require_basis_index(data::FRData, idx, a, b, c, role::Symbol)
    idx in hom_basis(data, a, b, c) ||
        error("Hom basis index $idx for $role is outside Hom($a ⊗ $b, $c)")
    return Int(idx)
end

function _frdata_metadata_lookup(data::FRData, table::Symbol, key)
    dict = get(data.metadata, table, nothing)
    dict === nothing && return nothing
    haskey(dict, key) && return dict[key]
    return nothing
end

function fr_value_one(data::FRData{T}) where {T}
    !isempty(data.F_values) && return one(data.F_values[1])
    !isempty(data.R_values) && return one(data.R_values[1])
    return one(T)
end

_fr_value_one(data::FRData) = fr_value_one(data)

"""
    F_symbol(frdata, a, b, c, d; e, f, μ=1, ν=1, ρ=1, σ=1)

Return the F-symbol entry with explicit object labels and Hom-basis indices.
The vector backend uses TensorCategories pentagon variable order.  Higher
multiplicity data can be supplied through `metadata[:F_symbols]` keyed by
`(a,b,c,d,e,f,μ,ν,ρ,σ)`.
"""
function F_symbol(data::FRData, a, b, c, d; e, f,
                  μ::Integer = 1, ν::Integer = 1,
                  ρ::Integer = 1, σ::Integer = 1)
    ia, ib, ic, id, ie, iff = _frdata_object_indices(data, a, b, c, d, e, f)
    key = (ia, ib, ic, id, ie, iff, Int(μ), Int(ν), Int(ρ), Int(σ))
    value = _frdata_metadata_lookup(data, :F_symbols, key)
    value !== nothing && return value

    _frdata_require_basis_index(data, μ, ia, ib, ie, :μ)
    _frdata_require_basis_index(data, ν, ie, ic, id, :ν)
    _frdata_require_basis_index(data, ρ, ib, ic, iff, :ρ)
    _frdata_require_basis_index(data, σ, ia, iff, id, :σ)

    if ia == 1 || ib == 1 || ic == 1
        one_value = fr_value_one(data)
        return (ie == iff && μ == ν == ρ == σ == 1) ? one_value : zero(one_value)
    end
    is_multiplicity_free(data.rules) ||
        error("vector-backed FRData F_symbol access is multiplicity-free only; provide metadata[:F_symbols] for multiplicityful data")
    idx = _f_value_index(data.rules, ia, ib, ic, id, ie, iff)
    idx === nothing &&
        error("missing F-symbol F^($ia,$ib,$ic)_$id[$ie,$iff] in FRData")
    return data.F_values[idx]
end

function has_F_symbol(data::FRData, a, b, c, d; e, f,
                      μ::Integer = 1, ν::Integer = 1,
                      ρ::Integer = 1, σ::Integer = 1)
    try
        F_symbol(data, a, b, c, d; e = e, f = f, μ = μ, ν = ν, ρ = ρ, σ = σ)
        return true
    catch
        return false
    end
end

"""
    R_symbol(frdata, a, b, c; μ=1, ν=1, inverse=false)

Return the R-symbol entry on `Hom(a ⊗ b, c)`.  Multiplicityful values can be
supplied through `metadata[:R_symbols]` keyed by `(a,b,c,μ,ν,inverse)`.
"""
function R_symbol(data::FRData, a, b, c; μ::Integer = 1, ν::Integer = 1,
                  inverse::Bool = false)
    ia, ib, ic = _frdata_object_indices(data, a, b, c)
    key = (ia, ib, ic, Int(μ), Int(ν), inverse)
    value = _frdata_metadata_lookup(data, :R_symbols, key)
    value !== nothing && return value

    _frdata_require_basis_index(data, μ, ia, ib, ic, :μ)
    _frdata_require_basis_index(data, ν, ia, ib, ic, :ν)
    if ia == 1 || ib == 1
        one_value = fr_value_one(data)
        return μ == ν ? one_value : zero(one_value)
    end
    is_multiplicity_free(data.rules) ||
        error("vector-backed FRData R_symbol access is multiplicity-free only; provide metadata[:R_symbols] for multiplicityful data")
    idx = _r_value_index(data.rules, ia, ib, ic)
    idx === nothing && error("missing R-symbol R^($ia,$ib)_$ic in FRData")
    return inverse ? data.R_inverse_values[idx] : data.R_values[idx]
end

function has_R_symbol(data::FRData, a, b, c; μ::Integer = 1, ν::Integer = 1,
                      inverse::Bool = false)
    try
        R_symbol(data, a, b, c; μ = μ, ν = ν, inverse = inverse)
        return true
    catch
        return false
    end
end

function _inverse_R_values_from_forward_values(fr::FusionRule, R_values::Vector{T}) where {T}
    positions, nR = _braiding_block_positions(fr.N)
    length(R_values) == nR || error("R_values has length $(length(R_values)); expected $nR")
    values = Vector{T}(undef, nR)
    for ((a, b, c), pos) in positions
        opp = positions[(b, a, c)]
        values[first(pos)] = inv(R_values[first(opp)])
    end
    return values
end

function _f_value_index(fr::FusionRule, a::Int, b::Int, c::Int, d::Int, e::Int, f::Int)
    _, _, nF = get_pentagon_system(fr.N, fr.rank)
    nF == 0 && return nothing
    for (idx, m) in enumerate(_pentagon_variable_metadata(fr.N, fr.rank, nF))
        (m.i, m.j, m.k, m.o, m.a, m.b) == (a, b, c, d, e, f) && return idx
    end
    return nothing
end

function _r_value_index(fr::FusionRule, a::Int, b::Int, c::Int)
    positions, _ = _braiding_block_positions(fr.N)
    pos = get(positions, (a, b, c), nothing)
    pos === nothing && return nothing
    length(pos) == 1 || error("multiplicity > 1 R-blocks are not supported by vector-backed FRData")
    return first(pos)
end

function _frdata_value_type(Fvals::AbstractVector, Rvals::AbstractVector)
    if !isempty(Fvals)
        return typeof(first(Fvals))
    elseif !isempty(Rvals)
        return typeof(first(Rvals))
    end
    return Int
end

function frdata_from_vectors(rules, Fvals::AbstractVector, Rvals::AbstractVector;
                             metadata = Dict{Symbol, Any}())
    fr = _fusion_rule(rules)
    require_multiplicity_free(fr)
    _, _, nF = get_pentagon_system(fr.N, fr.rank)
    length(Fvals) == nF ||
        error("F_values has length $(length(Fvals)); expected $nF in TensorCategories pentagon order")
    positions, nforward = _braiding_block_positions(fr.N)
    length(Rvals) >= nforward || error("R vector has length $(length(Rvals)); expected at least $nforward")
    length(Rvals) == nforward || length(Rvals) == 2 * nforward ||
        error("R_values has length $(length(Rvals)); expected $nforward or $(2 * nforward)")
    T = _frdata_value_type(Fvals, Rvals)
    Fvec = T[x for x in Fvals]
    Rvec = T[x for x in Rvals[1:nforward]]
    invvals = length(Rvals) >= 2 * nforward ?
              T[x for x in Rvals[(nforward + 1):(2 * nforward)]] :
              _inverse_R_values_from_forward_values(fr, Rvec)
    return FRData{T, typeof(metadata)}(fr, Fvec, Rvec, invvals, metadata)
end

_frdata_from_vectors(rules, Fvals::AbstractVector, Rvals::AbstractVector; kwargs...) =
    frdata_from_vectors(rules, Fvals, Rvals; kwargs...)

function FRData(rules::FusionRule,
                F_values::AbstractVector,
                R_values::AbstractVector,
                R_inverse_values::AbstractVector;
                metadata::Dict{Symbol, Any} = Dict{Symbol, Any}())
    return frdata_from_vectors(rules, F_values, vcat(collect(R_values), collect(R_inverse_values));
                               metadata = metadata)
end

function frdata_from_namedtuple(fr_tuple::NamedTuple;
                                metadata = Dict{Symbol, Any}(:source => :namedtuple))
    if hasproperty(fr_tuple, :rules) && hasproperty(fr_tuple, :F_values) &&
       hasproperty(fr_tuple, :R_values)
        Rinv = hasproperty(fr_tuple, :R_inverse_values) ? fr_tuple.R_inverse_values :
               similar(collect(fr_tuple.R_values), 0)
        return frdata_from_vectors(_fusion_rule(fr_tuple.rules), fr_tuple.F_values,
                                   vcat(collect(fr_tuple.R_values), collect(Rinv));
                                   metadata = metadata)
    end
    error("FR data NamedTuple must contain :rules, :F_values, and :R_values in TensorCategories order")
end

fr_hexagon_values(data::FRData) = vcat(F_values(data), R_values(data), R_inverse_values(data))
fr_pentagon_values(data::FRData) = F_values(data)

_fr_hexagon_values(data::FRData) = fr_hexagon_values(data)
_fr_pentagon_values(data::FRData) = fr_pentagon_values(data)
