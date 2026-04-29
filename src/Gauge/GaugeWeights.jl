"""
Toric gauge weights for multiplicity-free fusion rules.

The scalar gauge parameter attached to a nonzero fusion channel
`V_ab^c` is denoted `u[a,b,c]`.  For multiplicity-free symbols,

    F'^{abc}_{d;e,f} = u[a,b,e] u[e,c,d] / (u[b,c,f] u[a,f,d]) F^{abc}_{d;e,f]
    R'^{ab}_c       = u[b,a,c] / u[a,b,c] R^{ab}_c

so every symbol coordinate carries an integral character of the gauge
torus `prod G_m`.
"""

const GaugeParameter = Tuple{Int, Int, Int}
const GaugeWeight = Dict{GaugeParameter, Int}

function _fusion_array(fr::FusionRule)
    return fr.N
end

function _fusion_array(Nijk::Array{Int, 3})
    return Nijk
end

function _require_multiplicity_free_toric(Nijk::Array{Int,3})
    all(x -> x == 0 || x == 1, Nijk) ||
        error("toric gauge infrastructure currently supports multiplicity-free fusion rules only")
    return true
end

"""
    gauge_parameters(fusion_rule)

Return the ordered torus parameters `(a,b,c)` with `N_ab^c = 1`.
"""
function gauge_parameters(fusion_rule)
    Nijk = _fusion_array(fusion_rule)
    _require_multiplicity_free_toric(Nijk)
    r = size(Nijk, 1)
    params = GaugeParameter[]
    for a in 1:r, b in 1:r, c in 1:r
        Nijk[a, b, c] == 0 && continue
        push!(params, (a, b, c))
    end
    return params
end

function _add_weight!(w::GaugeWeight, ch::GaugeParameter, exponent::Int)
    w[ch] = get(w, ch, 0) + exponent
    w[ch] == 0 && delete!(w, ch)
    return w
end

"""
    f_symbol_weight(a,b,c,d,e,f)

Return the gauge character of `F^{abc}_{d;e,f}` as a dictionary from
fusion-channel parameters to integer exponents.
"""
function f_symbol_weight(a::Int, b::Int, c::Int, d::Int, e::Int, f::Int)
    w = GaugeWeight()
    _add_weight!(w, (a, b, e), 1)
    _add_weight!(w, (e, c, d), 1)
    _add_weight!(w, (b, c, f), -1)
    _add_weight!(w, (a, f, d), -1)
    return w
end

"""
    r_symbol_weight(a,b,c)

Return the gauge character of `R^{ab}_c`.
"""
function r_symbol_weight(a::Int, b::Int, c::Int)
    w = GaugeWeight()
    _add_weight!(w, (b, a, c), 1)
    _add_weight!(w, (a, b, c), -1)
    return w
end

function _admissible_f_coordinates(Nijk::Array{Int,3})
    r = size(Nijk, 1)
    coords = NamedTuple[]
    for a in 1:r, b in 1:r, c in 1:r, d in 1:r
        left = Int[e for e in 1:r if Nijk[a, b, e] != 0 && Nijk[e, c, d] != 0]
        right = Int[f for f in 1:r if Nijk[b, c, f] != 0 && Nijk[a, f, d] != 0]
        for e in left, f in right
            push!(coords, (kind = :F, a = a, b = b, c = c, d = d, e = e, f = f))
        end
    end
    return coords
end

function _admissible_r_coordinates(Nijk::Array{Int,3})
    r = size(Nijk, 1)
    coords = NamedTuple[]
    for a in 1:r, b in 1:r, c in 1:r
        Nijk[a, b, c] == 0 && continue
        push!(coords, (kind = :R, a = a, b = b, c = c))
    end
    return coords
end

"""
    symbol_coordinates(fusion_rule; include_R=true)

Return the ordered multiplicity-free F-symbol coordinates and, when
`include_R=true`, R-symbol coordinates determined by the fusion rule.
"""
function symbol_coordinates(fusion_rule; include_R::Bool = true)
    Nijk = _fusion_array(fusion_rule)
    _require_multiplicity_free_toric(Nijk)
    coords = _admissible_f_coordinates(Nijk)
    include_R && append!(coords, _admissible_r_coordinates(Nijk))
    return coords
end

function _coordinate_weight(coord)
    if coord.kind == :F
        return f_symbol_weight(coord.a, coord.b, coord.c, coord.d, coord.e, coord.f)
    elseif coord.kind == :R
        return r_symbol_weight(coord.a, coord.b, coord.c)
    end
    error("unsupported symbol coordinate kind $(coord.kind)")
end

"""
    gauge_weight_matrix(fusion_rule)

Return the integer matrix `W` whose rows are symbol coordinates and whose
columns are `gauge_parameters(fusion_rule)`.
"""
function gauge_weight_matrix(fusion_rule)
    params = gauge_parameters(fusion_rule)
    param_index = Dict(ch => i for (i, ch) in enumerate(params))
    coords = symbol_coordinates(fusion_rule)
    W = zeros(Int, length(coords), length(params))
    for (row, coord) in enumerate(coords)
        for (ch, exponent) in _coordinate_weight(coord)
            col = get(param_index, ch, 0)
            col != 0 && (W[row, col] += exponent)
        end
    end
    return W
end
