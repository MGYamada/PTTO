"""
FRData-centered gauge action API.
"""

function _gauge_field(data::FRData)
    return get(fr_metadata(data), :base_field, fr_scalar_type(data))
end

function _as_gauge_action(gauge; field = nothing,
                          metadata::Dict{Symbol, Any} = Dict{Symbol, Any}())
    gauge isa GaugeAction && return gauge
    gauge isa GaugeParameters && return GaugeAction(gauge; field = field,
                                                    metadata = merge(copy(gauge.metadata), metadata))
    gauge isa GaugeTransform && return GaugeAction(gauge; field = field, metadata = metadata)
    if gauge isa AbstractDict
        isempty(gauge) && error("cannot infer scalar type from an empty gauge dictionary")
        first_key = first(keys(gauge))
        if first_key isa Tuple{Int,Int,Int,Int}
            return GaugeAction(Dict(k => v for (k, v) in gauge);
                               field = field, metadata = metadata)
        elseif first_key isa Tuple{Int,Int,Int}
            return GaugeAction(Dict(k => v for (k, v) in gauge);
                               field = field, metadata = metadata)
        end
    end
    error("unsupported gauge object: $(typeof(gauge))")
end

function _gauge_action_scalars(action::GaugeAction)
    return Dict((a, b, c) => value
                for ((a, b, c, μ), value) in action.parameters if μ == 1)
end

"""
    identity_gauge(frdata)
    identity_gauge(fusion_rule; one_value=1)

Return the identity `GaugeAction` on all admissible scalar gauge coordinates.
"""
function identity_gauge(data::FRData)
    one_value = fr_value_one(data)
    params = Dict((a, b, c, μ) => one_value
                  for dof in gauge_degrees_of_freedom(data)
                  for (a, b, c, μ) in ((dof.a, dof.b, dof.c, dof.μ),))
    return GaugeAction(params; field = _gauge_field(data),
                       metadata = Dict{Symbol, Any}(:kind => :identity,
                                                    :fusion_rank => fusion_rule(data).rank))
end

function identity_gauge(rules; one_value = 1, field = nothing)
    params = Dict((a, b, c, 1) => one_value for (a, b, c) in gauge_parameters(rules))
    return GaugeAction(params; field = field,
                       metadata = Dict{Symbol, Any}(:kind => :identity,
                                                    :fusion_rank => _gauge_rules(rules).rank))
end

"""
    inverse_gauge(gauge)

Return the inverse scalar gauge action.
"""
function inverse_gauge(gauge)
    action = _as_gauge_action(gauge)
    params = Dict(key => inv(value) for (key, value) in action.parameters)
    meta = copy(action.metadata)
    meta[:inverse_of] = get(action.metadata, :kind, :gauge_action)
    return GaugeAction(params; field = action.field, metadata = meta)
end

"""
    compose_gauge(first, second)

Compose two gauge actions so that
`apply_gauge(apply_gauge(fr, first), second) == apply_gauge(fr,
compose_gauge(first, second))`.
"""
function compose_gauge(first, second)
    g1 = _as_gauge_action(first)
    g2 = _as_gauge_action(second)
    keys12 = union(keys(g1.parameters), keys(g2.parameters))
    params = Dict{Tuple{Int,Int,Int,Int}, Any}()
    for key in keys12
        if haskey(g1.parameters, key) && haskey(g2.parameters, key)
            params[key] = g1.parameters[key] * g2.parameters[key]
        elseif haskey(g1.parameters, key)
            params[key] = g1.parameters[key]
        else
            params[key] = g2.parameters[key]
        end
    end
    value_types = map(typeof, collect(values(params)))
    T = isempty(value_types) ? Any : promote_type(value_types...)
    typed = Dict{Tuple{Int,Int,Int,Int}, T}(key => value for (key, value) in params)
    field = g2.field === nothing ? g1.field : g2.field
    meta = Dict{Symbol, Any}(:kind => :composition,
                             :left_kind => get(g1.metadata, :kind, :gauge_action),
                             :right_kind => get(g2.metadata, :kind, :gauge_action))
    return GaugeAction(typed; field = field, metadata = meta)
end

function _normalize_gauge_target(target::Symbol)
    target in (:FR, :F, :R) ||
        error("target must be one of :FR, :F, or :R; got $(repr(target))")
    return target
end

"""
    apply_gauge(frdata, gauge; target=:FR)

Apply a gauge action to `FRData`.  `target=:FR` transforms both F- and
R-symbols, while `target=:F` or `target=:R` applies only the corresponding
part of the action.  `gauge` may be `GaugeAction`, `GaugeParameters`,
`GaugeTransform`, a tuple-keyed dictionary, or `nothing`.
"""
function apply_gauge(data::FRData, gauge; target::Symbol = :FR)
    target = _normalize_gauge_target(target)
    gauge === nothing && return target == :FR ? gauge_transform(data, nothing) :
        frdata_from_vectors(fusion_rule(data), F_values(data),
                            vcat(R_values(data), R_inverse_values(data));
                            metadata = copy(fr_metadata(data)))

    validate_frdata_for_gauge(data)
    one_value = fr_value_one(data)
    R_all = vcat(R_values(data), R_inverse_values(data))
    Fout = target in (:FR, :F) ? _transform_F_values(F_values(data), data, gauge, one_value) :
           copy(F_values(data))
    Rout = target in (:FR, :R) ? _transform_R_values(R_all, data, gauge, one_value) :
           copy(R_all)
    metadata = copy(fr_metadata(data))
    metadata[:gauge_action] = gauge
    metadata[:gauge_action_target] = target
    return frdata_from_vectors(fusion_rule(data), Fout, Rout; metadata = metadata)
end

"""
    apply_gauge(F, R, gauge; Nijk, target=:FR)

Legacy vector wrapper for applying a gauge action to raw F/R coordinate
vectors.
"""
function apply_gauge(F, R, gauge; Nijk::Union{Array{Int,3}, Nothing} = nothing,
                     target::Symbol = :FR)
    target = _normalize_gauge_target(target)
    Nijk === nothing && error("Nijk is required for vector-backed gauge actions")
    fr = FusionRule(Nijk)
    _require_multiplicity_free_gauge(fr)
    one_value = _one_like_from_FR(F, R)
    Fout = target in (:FR, :F) ? _transform_F_values(F, fr, gauge, one_value) : copy(F)
    Rout = target in (:FR, :R) ? _transform_R_values(R, fr, gauge, one_value) : copy(R)
    return (F = Fout, R = Rout)
end

gauge_transform(data::FRData, gauge::GaugeAction) = apply_gauge(data, gauge)
