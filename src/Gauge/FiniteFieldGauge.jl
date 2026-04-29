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

function _gauge_value_mod_p(gauge_element, ch::GaugeParameter, idx::Int, p::Int)
    if gauge_element isa AbstractDict
        return mod(Int(get(gauge_element, ch, 1)), p)
    elseif gauge_element isa GaugeTransform
        return mod(Int(get(gauge_element.scalars, ch, 1)), p)
    elseif gauge_element isa AbstractVector
        idx <= length(gauge_element) || error("gauge vector is shorter than the parameter list")
        return mod(Int(gauge_element[idx]), p)
    elseif hasproperty(gauge_element, :parameters) && hasproperty(gauge_element, :values)
        pos = findfirst(==(ch), gauge_element.parameters)
        pos === nothing && return 1
        return mod(Int(gauge_element.values[pos]), p)
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
        out[i] = mod(Int(values[i]) * _coordinate_factor_mod_p(coord, params, gauge_element, p), p)
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
        return Dict(coord => mod(Int(value) *
                                 _coordinate_factor_mod_p(coord, params, gauge_element, p), p)
                    for (coord, value) in symbol_data)
    elseif hasproperty(symbol_data, :F) || hasproperty(symbol_data, :R)
        Fout = hasproperty(symbol_data, :F) ?
            Dict(coord => mod(Int(value) *
                              _coordinate_factor_mod_p(coord, params, gauge_element, p), p)
                 for (coord, value) in symbol_data.F) :
            Dict{NamedTuple, Int}()
        Rout = hasproperty(symbol_data, :R) ?
            Dict(coord => mod(Int(value) *
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
