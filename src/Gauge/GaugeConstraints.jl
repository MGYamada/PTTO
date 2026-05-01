"""
Gauge-fixing constraint construction and solving.
"""

"""
    gauge_degrees_of_freedom(frdata)
    gauge_degrees_of_freedom(fusion_rule)

Enumerate scalar gauge coordinates in deterministic channel order.
"""
function gauge_degrees_of_freedom(data::FRData)
    validate_frdata_for_gauge(data)
    out = GaugeDegreeOfFreedom[]
    for a in simples(data), b in simples(data), c in fusion_channels(data, a, b)
        ia, ib, ic = _frdata_object_indices(data, a, b, c)
        for μ in hom_basis(data, ia, ib, ic)
            push!(out, GaugeDegreeOfFreedom(ia, ib, ic, μ))
        end
    end
    sort!(out; by = x -> (x.a, x.b, x.c, x.μ))
    return out
end

function gauge_degrees_of_freedom(rules)
    fr = _gauge_rules(rules)
    _require_multiplicity_free_gauge(fr)
    return [GaugeDegreeOfFreedom(a, b, c, 1) for (a, b, c) in gauge_parameters(fr)]
end

function _constraints_have(::Type{T}, constraints) where {T<:GaugeConstraint}
    return any(c -> c isa T, constraints)
end

function _selected_f_indices_from_plan(data::FRData)
    return Int[entry.var_idx for entry in gauge_fixing_plan(data)]
end

"""
    build_gauge_constraints(frdata; strategy=:default)

Build the explicit constraints used by the normal-form layer.  `:default`
selects the current deterministic F-symbol pivot plan; `:none` records only
the unit-channel convention.
"""
function build_gauge_constraints(data::FRData; strategy::Symbol = :default)
    validate_frdata_for_gauge(data)
    constraints = GaugeConstraint[FixUnitConstraints()]
    if strategy == :default
        push!(constraints, FixSelectedFSymbols(_selected_f_indices_from_plan(data)))
    elseif strategy == :none
        # The none strategy intentionally records no coordinate substitutions.
    else
        error("unknown gauge-fixing strategy $(repr(strategy)); expected :default or :none")
    end
    push!(constraints, NormalizationConstraint(strategy))
    return constraints
end

function build_gauge_constraints(rules; strategy::Symbol = :none)
    strategy == :none ||
        error("fusion-rule gauge constraints without FRData currently support only strategy=:none")
    _require_multiplicity_free_gauge(rules)
    return GaugeConstraint[FixUnitConstraints(), NormalizationConstraint(:none)]
end

"""
    solve_gauge_constraints(frdata, constraints; field=nothing)

Solve explicit gauge constraints and return a `GaugeAction`.  The current
solver delegates deterministic F-symbol normalization to the existing exact
multiplicity-free pivot algorithm; R-symbol constraints are recorded for
future strategies and currently rejected.
"""
function solve_gauge_constraints(data::FRData, constraints::AbstractVector{<:GaugeConstraint};
                                 field = _gauge_field(data))
    validate_frdata_for_gauge(data)
    if _constraints_have(FixSelectedRSymbols, constraints)
        error("FixSelectedRSymbols constraints are not implemented yet")
    end
    if _constraints_have(FixSelectedFSymbols, constraints)
        transform = _canonical_gauge_transform(F_values(data), R_values(data), data)
        meta = Dict{Symbol, Any}(:kind => :normal_form,
                                 :constraints => collect(constraints),
                                 :fixed_indices => transform.fixed_indices,
                                 :complete => transform.complete)
        return GaugeAction(transform; field = field, metadata = meta)
    end
    action = identity_gauge(data)
    meta = copy(action.metadata)
    meta[:constraints] = collect(constraints)
    return GaugeAction(copy(action.parameters); field = field, metadata = meta)
end
