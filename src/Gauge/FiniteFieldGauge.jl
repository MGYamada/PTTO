"""
Finite-field toric gauge action on multiplicity-free symbol coordinates.
"""

function _inv_mod_nonzero(a::Integer, p::Integer)
    a = mod(a, p)
    a == 0 && error("cannot invert 0 in F_$p")
    return Int(invmod(a, p))
end

function _pow_mod_signed(a::Integer, exponent::Integer, p::Integer)
    a = mod(a, p)
    exponent == 0 && return 1
    exponent > 0 && return Int(powermod(a, exponent, p))
    return Int(powermod(_inv_mod_nonzero(a, p), -exponent, p))
end

function _modp_scalar_value(x, p::Int)
    if x isa FpElem
        x.p == p || error("cannot use F_$(x.p) scalar over F_$p")
        return x.value
    end
    return mod(Int(x), p)
end

function _gauge_value_mod_p(gauge_element, ch::GaugeParameter, idx::Int, p::Int)
    if gauge_element isa AbstractDict
        return _modp_scalar_value(get(gauge_element, ch, 1), p)
    elseif gauge_element isa GaugeTransform
        return _modp_scalar_value(get(gauge_element.scalars, ch, 1), p)
    elseif gauge_element isa GaugeParameters
        return _modp_scalar_value(get(gauge_element.values, (ch[1], ch[2], ch[3], 1), 1), p)
    elseif gauge_element isa GaugeAction
        return _modp_scalar_value(get(gauge_element.parameters, (ch[1], ch[2], ch[3], 1), 1), p)
    elseif gauge_element isa AbstractVector
        idx <= length(gauge_element) || error("gauge vector is shorter than the parameter list")
        return _modp_scalar_value(gauge_element[idx], p)
    elseif hasproperty(gauge_element, :parameters) && hasproperty(gauge_element, :values)
        pos = findfirst(==(ch), gauge_element.parameters)
        pos === nothing && return 1
        return _modp_scalar_value(gauge_element.values[pos], p)
    end
    error("unsupported gauge element $(typeof(gauge_element))")
end

function _coordinate_factor_mod_p(coord, params, gauge_element, p::Int)
    param_index = Dict(ch => i for (i, ch) in enumerate(params))
    factor = 1
    for (ch, exponent) in _coordinate_weight(coord)
        idx = get(param_index, ch, 0)
        idx == 0 && continue
        value = _gauge_value_mod_p(gauge_element, ch, idx, p)
        factor = mod(factor * _pow_mod_signed(value, exponent, p), p)
    end
    return factor
end

function _infer_parameters_for_symbol_data(symbol_data)
    if hasproperty(symbol_data, :parameters)
        return collect(symbol_data.parameters)
    end
    params = Set{GaugeParameter}()
    coords = if hasproperty(symbol_data, :coordinates)
        symbol_data.coordinates
    elseif symbol_data isa AbstractDict
        keys(symbol_data)
    elseif hasproperty(symbol_data, :F) || hasproperty(symbol_data, :R)
        Iterators.flatten((hasproperty(symbol_data, :F) ? keys(symbol_data.F) : NamedTuple[],
                           hasproperty(symbol_data, :R) ? keys(symbol_data.R) : NamedTuple[]))
    else
        error("symbol_data must carry coordinates or F/R dictionaries")
    end
    for coord in coords
        for ch in keys(_coordinate_weight(coord))
            push!(params, ch)
        end
    end
    return sort(collect(params))
end

function _apply_values_mod_p(values, coords, params, gauge_element, p::Int)
    length(values) == length(coords) || error("values and coordinates length mismatch")
    out = similar(collect(values), Int)
    for (i, coord) in enumerate(coords)
        out[i] = mod(_modp_scalar_value(values[i], p) *
                     _coordinate_factor_mod_p(coord, params, gauge_element, p), p)
    end
    return out
end

"""
    apply_gauge_mod_p(symbol_data, gauge_element, p)

Apply a toric gauge element over `F_p`.  Supported `symbol_data` shapes are:
`(coordinates=..., values=..., parameters=...)`, a dictionary keyed by
symbol coordinates, or `(F=Dict(...), R=Dict(...), parameters=...)`.
"""
function apply_gauge_mod_p(symbol_data, gauge_element, p::Integer)
    p = Int(p)
    params = _infer_parameters_for_symbol_data(symbol_data)
    if hasproperty(symbol_data, :coordinates) && hasproperty(symbol_data, :values)
        values = _apply_values_mod_p(symbol_data.values, symbol_data.coordinates,
                                     params, gauge_element, p)
        return (coordinates = collect(symbol_data.coordinates),
                values = values,
                parameters = params)
    elseif symbol_data isa AbstractDict
        return Dict(coord => mod(_modp_scalar_value(value, p) *
                                 _coordinate_factor_mod_p(coord, params, gauge_element, p), p)
                    for (coord, value) in symbol_data)
    elseif hasproperty(symbol_data, :F) || hasproperty(symbol_data, :R)
        Fout = hasproperty(symbol_data, :F) ?
            Dict(coord => mod(_modp_scalar_value(value, p) *
                              _coordinate_factor_mod_p(coord, params, gauge_element, p), p)
                 for (coord, value) in symbol_data.F) :
            Dict{NamedTuple, Int}()
        Rout = hasproperty(symbol_data, :R) ?
            Dict(coord => mod(_modp_scalar_value(value, p) *
                              _coordinate_factor_mod_p(coord, params, gauge_element, p), p)
                 for (coord, value) in symbol_data.R) :
            Dict{NamedTuple, Int}()
        return (F = Fout, R = Rout, parameters = params)
    end
    error("unsupported symbol_data $(typeof(symbol_data))")
end

function _symbol_data_active_coordinates(symbol_data)
    if hasproperty(symbol_data, :coordinates) && hasproperty(symbol_data, :values)
        return [coord for (coord, value) in zip(symbol_data.coordinates, symbol_data.values)
                if !iszero(value)]
    elseif symbol_data isa AbstractDict
        return [coord for (coord, value) in symbol_data if !iszero(value)]
    elseif hasproperty(symbol_data, :F) || hasproperty(symbol_data, :R)
        coords = NamedTuple[]
        if hasproperty(symbol_data, :F)
            append!(coords, [coord for (coord, value) in symbol_data.F if !iszero(value)])
        end
        if hasproperty(symbol_data, :R)
            append!(coords, [coord for (coord, value) in symbol_data.R if !iszero(value)])
        end
        return coords
    end
    error("unsupported symbol_data $(typeof(symbol_data))")
end

function _weight_matrix_for_coordinates(coords, params)
    param_index = Dict(ch => i for (i, ch) in enumerate(params))
    W = zeros(Int, length(coords), length(params))
    for (row, coord) in enumerate(coords)
        for (ch, exponent) in _coordinate_weight(coord)
            col = get(param_index, ch, 0)
            col != 0 && (W[row, col] += exponent)
        end
    end
    return W
end

"""
    stabilizer_size_mod_p(symbol_data, fusion_rule, p)

Return the size of the stabilizer in `(F_p^*)^n` of the nonzero symbol
coordinates.  Zero coordinates impose no character equation.
"""
function stabilizer_size_mod_p(symbol_data, fusion_rule, p::Integer)
    q = BigInt(p - 1)
    params = gauge_parameters(fusion_rule)
    coords = _symbol_data_active_coordinates(symbol_data)
    W = _weight_matrix_for_coordinates(coords, params)
    invariants = _smith_invariant_factors(W)
    rankW = length(invariants)
    size_kernel = q^(length(params) - rankW)
    for d in invariants
        size_kernel *= gcd(q, BigInt(d))
    end
    return size_kernel
end

"""
    stacky_weight_mod_p(symbol_data, fusion_rule, p)

Return the groupoid-counting weight `1 / |Stab(symbol_data)|`.
"""
function stacky_weight_mod_p(symbol_data, fusion_rule, p::Integer)
    return 1 // stabilizer_size_mod_p(symbol_data, fusion_rule, p)
end

function _toric_prime(data::FRData{FpElem})
    return _frdata_prime(data)
end

function _toric_prime(data::FRData)
    p = get(fr_metadata(data), :p, get(fr_metadata(data), :prime, nothing))
    p === nothing && return nothing
    return Int(p)
end

function _toric_symbol_data(data::FRData; include_R::Bool = true)
    coords = fr_symbol_coordinates(data; include_R = include_R)
    values = _toric_symbol_values(data; include_R = include_R)
    return (coordinates = coords,
            values = values,
            parameters = gauge_parameters(data))
end

"""
    toric_gauge_normal_form(frdata; include_R=true, verify=true)

Inspect the pre-solver Smith-normal-form toric F-slice carried by `FRData`.
This function does not apply a new post-hoc gauge transform and does not fix
R-symbols.  The actual gauge fixing happens before F/R solving through
`gauge_fix(fr_equation_system(...))`; this wrapper reports the slice,
stabilizer, and stacky metadata for downstream braid and Gauss-sum consumers.
"""
function toric_gauge_normal_form(data::FRData; include_R::Bool = true,
                                 verify::Bool = true)
    validate_frdata_for_gauge(data)
    toric = toric_gauge_data(data; include_R = include_R)
    solver_gauge = get(fr_metadata(data), :gauge_fixing, Dict{Symbol, Any}())
    slice = get(solver_gauge, :toric_gauge_slice,
                toric_gauge_slice(fusion_rule(data); coordinate_kind = :F))
    fixed_indices = Int[i for i in get(slice, :fixed_indices, Int[])]
    complete = all(i -> F_values(data)[i] == fr_value_one(data), fixed_indices)
    action = identity_gauge(data)
    p = _toric_prime(data)
    symbol_data = _toric_symbol_data(data; include_R = include_R)
    stabilizer = p === nothing ? nothing :
        stabilizer_size_mod_p(symbol_data, fusion_rule(data), p)
    stacky = stabilizer === nothing ? nothing : 1 // stabilizer
    verification = verify && data isa FRData{FpElem} ? verify_FRData(data) :
                   verify ? validate_frdata_for_gauge(data) : nothing
    meta = Dict{Symbol, Any}(:kind => :toric_snf_normal_form,
                             :include_R => include_R,
                             :parameters => toric.parameters,
                             :coordinates => toric.coordinates,
                             :weight_matrix => toric.weight_matrix,
                             :smith_split => get(slice, :split, toric.split),
                             :fixed_indices => fixed_indices,
                             :complete => complete,
                             :stage => :fr_solve_preprocessing,
                             :gauge_convention => :full_channel_scalar,
                             :gauge_group_kind => :full_channel_toric_gauge,
                             :includes_unit_channels => true,
                             :includes_ineffective_kernel => true,
                             :stabilizer_size => stabilizer,
                             :stacky_weight => stacky,
                             :verification => verification)
    frmeta = copy(fr_metadata(data))
    frmeta[:toric_gauge_normal_form] = meta
    reported = frdata_from_vectors(fusion_rule(data), F_values(data),
                                   vcat(R_values(data), R_inverse_values(data));
                                   metadata = frmeta)
    return ToricGaugeNormalFormResult(reported, action, fixed_indices,
                                      get(slice, :split, toric.split),
                                      stabilizer, stacky, meta)
end

toric_gauge_normal_form_mod_p(data::FRData{FpElem}; kwargs...) =
    toric_gauge_normal_form(data; kwargs...)
