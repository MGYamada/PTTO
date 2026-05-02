"""
    ToricGaugeFixingError(message)

Error raised when a requested toric gauge-fixing step is mathematically or
arithmetically unsupported.
"""
struct ToricGaugeFixingError <: Exception
    message::String
end

Base.showerror(io::IO, err::ToricGaugeFixingError) = print(io, err.message)

"""
    ToricGaugeData

Structured toric gauge preconditioning data for multiplicity-free fusion
rules.  `parameters` are the trivalent channels carrying torus coordinates,
and `slice` records the deterministic F-symbol normalization used before
F/R reconstruction.
"""
struct ToricGaugeData
    fusion::FusionRule
    parameters::Vector{GaugeParameter}
    channels::Vector{GaugeParameter}
    slice::NamedTuple
    field::Any
    conventions::Any
    metadata::Dict{Symbol, Any}
end

"""
    Smith normal form summaries for toric gauge actions.
"""

function _smith_invariant_factors(W::AbstractMatrix{<:Integer})
    m, n = size(W)
    (m == 0 || n == 0) && return BigInt[]
    entries = [BigInt(W[i, j]) for i in 1:m for j in 1:n]
    S = AbstractAlgebra.snf(matrix(ZZ, m, n, entries))
    invariants = BigInt[]
    for i in 1:min(m, n)
        d = abs(BigInt(S[i, i]))
        d == 0 && continue
        push!(invariants, d)
    end
    return invariants
end

"""
    smith_gauge_split(W)

Summarize the Smith normal form of the integer gauge weight matrix.  The
returned named tuple records the invariant factors, rank of the effective
torus action, free rank of the ineffective kernel, and finite residual
orders from non-unit invariant factors.
"""
function smith_gauge_split(W::AbstractMatrix{<:Integer})
    invariants = _smith_invariant_factors(W)
    rankW = length(invariants)
    nparams = size(W, 2)
    orders = Int[d for d in invariants if d > 1]
    return (invariant_factors = Int.(invariants),
            effective_rank = rankW,
            ineffective_rank = nparams - rankW,
            residual_orders = orders)
end

"""
    ineffective_kernel_rank(W)

Return the dimension of the connected kernel of the toric gauge action.
"""
ineffective_kernel_rank(W::AbstractMatrix{<:Integer}) = smith_gauge_split(W).ineffective_rank

"""
    residual_gauge_orders(W)

Return the finite diagonal kernel orders visible in the Smith normal form.
"""
residual_gauge_orders(W::AbstractMatrix{<:Integer}) = smith_gauge_split(W).residual_orders

function _fr_f_symbol_coordinates_tensorcategories(fr::FusionRule)
    _, _, nF = get_pentagon_system(fr.N, fr.rank)
    nF == 0 && return NamedTuple[]
    return NamedTuple[(kind = :F, a = m.i, b = m.j, c = m.k, d = m.o, e = m.a, f = m.b)
                      for m in _pentagon_variable_metadata(fr.N, fr.rank, nF)]
end

function _fr_r_symbol_coordinates_tensorcategories(fr::FusionRule)
    positions, nR = _braiding_block_positions(fr.N)
    coords = Vector{NamedTuple}(undef, nR)
    for ((a, b, c), pos) in positions
        length(pos) == 1 ||
            error("toric gauge coordinates currently support multiplicity-free R blocks only")
        coords[first(pos)] = (kind = :R, a = a, b = b, c = c)
    end
    return coords
end

"""
    fr_symbol_coordinates(frdata; include_R=true)

Return symbol-coordinate descriptors in the same vector order used by
`F_values(frdata)` followed by `R_values(frdata)`.
"""
function fr_symbol_coordinates(data::FRData; include_R::Bool = true)
    fr = fusion_rule(data)
    coords = _fr_f_symbol_coordinates_tensorcategories(fr)
    include_R && append!(coords, _fr_r_symbol_coordinates_tensorcategories(fr))
    return coords
end

function fr_symbol_coordinates(rules::Union{FusionRule, Array{Int, 3}}; include_R::Bool = true)
    fr = _fusion_rule(rules)
    coords = _fr_f_symbol_coordinates_tensorcategories(fr)
    include_R && append!(coords, _fr_r_symbol_coordinates_tensorcategories(fr))
    return coords
end

function _toric_symbol_values(data::FRData; include_R::Bool = true)
    return include_R ? vcat(F_values(data), R_values(data)) : copy(F_values(data))
end

"""
    toric_gauge_data(frdata; include_R=true)

Return the ordered toric gauge parameters, symbol coordinates, character
matrix, and Smith-normal-form split for multiplicity-free F/R data.
"""
function toric_gauge_data(data::FRData; include_R::Bool = true)
    coords = fr_symbol_coordinates(data; include_R = include_R)
    params = gauge_parameters(data)
    W = _weight_matrix_for_coordinates(coords, params)
    return (parameters = params,
            coordinates = coords,
            weight_matrix = W,
            split = smith_gauge_split(W))
end

"""
    toric_gauge_data(fusion; field=nothing, conventions=:tensorcategories)

Return deterministic toric gauge preconditioning data for a multiplicity-free
fusion rule.  The object records all channels `(a,b,c)` with `N_ab^c = 1`
and the reproducible F-symbol slice used by the reconstruction layer.
"""
function toric_gauge_data(rules::Union{FusionRule, Array{Int, 3}};
                          field = nothing,
                          conventions = :tensorcategories,
                          coordinate_kind::Symbol = :F)
    is_multiplicity_free(rules) ||
        throw(ToricGaugeFixingError("toric gauge fixing requires multiplicity-free fusion rules; found a fusion coefficient greater than 1"))
    fr = _fusion_rule(rules)
    coordinate_kind == :F ||
        throw(ToricGaugeFixingError("toric gauge fixing currently supports only F-symbol preconditioning; got coordinate_kind=$(repr(coordinate_kind))"))
    slice = toric_gauge_slice(fr; coordinate_kind = coordinate_kind)
    params = gauge_parameters(fr)
    metadata = Dict{Symbol, Any}(
        :gauge_fix_method => :toric_snf_f_slice,
        :fixed_f_indices => slice.fixed_indices,
        :effective_rank => slice.split.effective_rank,
        :complete => slice.complete,
        :nonzero_parameter_domain => field === nothing ? :torus : Symbol("$(field)_units"),
        :gauge_convention => :full_channel_scalar,
        :gauge_group_kind => :full_channel_toric_gauge,
        :includes_unit_channels => true,
        :includes_ineffective_kernel => true,
    )
    return ToricGaugeData(fr, params, copy(params), slice, field, conventions, metadata)
end

function Base.show(io::IO, data::ToricGaugeData)
    print(io, "ToricGaugeData(rank=$(data.fusion.rank), parameters=$(length(data.parameters)), fixed_f=$(length(data.slice.fixed_indices)))")
end

function _rank_int_matrix_rows(rows::Vector{Vector{Int}}, ncols::Int)
    isempty(rows) && return 0
    return rank(matrix(QQ, length(rows), ncols, vcat(rows...)))
end

function _independent_toric_coordinate_rows(W::AbstractMatrix{<:Integer},
                                            target_rank::Int)
    selected = Int[]
    rows = Vector{Int}[]
    current = 0
    for i in 1:size(W, 1)
        row = Int[W[i, j] for j in 1:size(W, 2)]
        any(!iszero, row) || continue
        trial = vcat(rows, [row])
        r = _rank_int_matrix_rows(trial, size(W, 2))
        if r > current
            push!(selected, i)
            push!(rows, row)
            current = r
            current == target_rank && break
        end
    end
    return selected
end

"""
    toric_gauge_slice(rules; coordinate_kind=:F)

Choose a deterministic pre-solver toric gauge slice for a multiplicity-free
fusion rule.  The returned indices are TensorCategories F-coordinate indices
to be fixed to `1`; Smith-normal-form metadata records the effective torus
rank and residual finite kernel for the F-symbol character matrix.
"""
function toric_gauge_slice(rules; coordinate_kind::Symbol = :F)
    coordinate_kind == :F ||
        error("only coordinate_kind=:F is supported for pre-solver toric gauge slices")
    fr = _fusion_rule(rules)
    _require_multiplicity_free_toric(fr)
    coords = fr_symbol_coordinates(fr; include_R = false)
    params = _ordered_gauge_channels(fr, fr.rank)
    nF = length(coords)
    channel_index = Dict(ch => i for (i, ch) in enumerate(params))
    metadata = nF == 0 ? NamedTuple[] : _pentagon_variable_metadata(fr.N, fr.rank, nF)
    W = zeros(Int, nF, length(params))
    for i in 1:nF
        weight = _gauge_weight(metadata[i], channel_index)
        for j in eachindex(weight)
            W[i, j] = weight[j]
        end
    end
    split = smith_gauge_split(W)
    fixed = _independent_toric_coordinate_rows(W, split.effective_rank)
    complete = length(fixed) == split.effective_rank
    return (coordinate_kind = coordinate_kind,
            parameters = params,
            coordinates = coords,
            weight_matrix = W,
            split = split,
            fixed_indices = fixed,
            complete = complete)
end
