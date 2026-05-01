"""
Validation helpers for gauge actions and gauge-fixed FRData.
"""

function _gauge_action_channel_key(key)
    a, b, c, μ = key
    μ == 1 || error("scalar GaugeAction currently supports only μ = 1, got μ = $μ")
    return (a, b, c)
end

function _is_nonzero_gauge_value(value)
    try
        return !iszero(value)
    catch
        return value != zero(value)
    end
end

"""
    validate_gauge_action(frdata, gauge; require_complete=true)

Check that a gauge action is compatible with the scalar multiplicity-free
gauge layer: all keys are admissible Hom-basis coordinates, all values are
nonzero, and, by default, every admissible channel has an explicit parameter.
"""
function validate_gauge_action(data::FRData, gauge; require_complete::Bool = true)
    validate_frdata_for_gauge(data)
    action = _as_gauge_action(gauge; field = _gauge_field(data))
    allowed = Set((dof.a, dof.b, dof.c, dof.μ) for dof in gauge_degrees_of_freedom(data))
    for (key, value) in action.parameters
        key in allowed || error("gauge parameter $key is not an admissible Hom-basis channel")
        _gauge_action_channel_key(key)
        _is_nonzero_gauge_value(value) || error("gauge parameter $key is zero")
    end
    if require_complete
        missing = setdiff(allowed, Set(keys(action.parameters)))
        isempty(missing) || error("gauge action is missing parameters for $(collect(missing))")
    end
    return true
end

function validate_gauge_action(rules::Union{FusionRule, Array{Int,3}}, gauge;
                               require_complete::Bool = true)
    fr = _gauge_rules(rules)
    _require_multiplicity_free_gauge(fr)
    action = _as_gauge_action(gauge)
    allowed = Set((a, b, c, 1) for (a, b, c) in gauge_parameters(fr))
    for (key, value) in action.parameters
        key in allowed || error("gauge parameter $key is not an admissible channel")
        _is_nonzero_gauge_value(value) || error("gauge parameter $key is zero")
    end
    if require_complete
        missing = setdiff(allowed, Set(keys(action.parameters)))
        isempty(missing) || error("gauge action is missing parameters for $(collect(missing))")
    end
    return true
end

"""
    validate_frdata(frdata)

Validate the structural assumptions ACMG's current scalar gauge layer requires.
For finite-field `FRData{FpElem}`, this also re-evaluates pentagon and hexagon
equations.
"""
function validate_frdata(data::FRData)
    validate_frdata_for_gauge(data)
    data isa FRData{FpElem} && return verify_FRData(data)
    return true
end

"""
    validate_gauge_fixed(frdata)

Return whether `frdata` satisfies the current deterministic gauge-fixed
F-symbol conditions and passes available FRData validation.
"""
function validate_gauge_fixed(data::FRData)
    validate_frdata(data) || return false
    return is_gauge_fixed(data)
end
