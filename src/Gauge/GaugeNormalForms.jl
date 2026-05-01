"""
Normal-form orchestration for gauge fixing.
"""

function _strategy_from_constraints(constraints)
    idx = findlast(c -> c isa NormalizationConstraint, constraints)
    idx === nothing && return :custom
    return constraints[idx].strategy
end

"""
    gauge_normal_form(frdata; strategy=:default, constraints=nothing, validate=true)

Build an explicit constraint list, solve for a `GaugeAction`, apply it to
`FRData`, and optionally validate the result.  This separates representative
selection from the low-level F/R action so future finite-field and GIT-style
strategies can replace only the constraint solver.
"""
function gauge_normal_form(data::FRData; strategy::Symbol = :default,
                           constraints = nothing, validate::Bool = true)
    validate_frdata_for_gauge(data)
    chosen = constraints === nothing ? build_gauge_constraints(data; strategy = strategy) :
             collect(constraints)
    action = solve_gauge_constraints(data, chosen; field = _gauge_field(data))
    transformed = apply_gauge(data, action)
    fixed_indices = Int[]
    complete = false
    if haskey(action.metadata, :fixed_indices)
        fixed_indices = Int[i for i in action.metadata[:fixed_indices]]
        complete = Bool(action.metadata[:complete])
    elseif _constraints_have(FixSelectedFSymbols, chosen)
        fixed_indices = _selected_f_indices_from_plan(data)
        complete = true
    end
    if validate
        validate_frdata(transformed) ||
            error("gauge normal form failed FRData validation")
        _constraints_have(FixSelectedFSymbols, chosen) && !is_gauge_fixed(transformed) &&
            error("gauge normal form did not satisfy selected F-symbol constraints")
    end
    meta = Dict{Symbol, Any}(:frdata => transformed,
                             :constraints => chosen,
                             :strategy => _strategy_from_constraints(chosen))
    return GaugeFixingResult(F_values(transformed), R_values(transformed), action,
                             fixed_indices, complete, meta)
end

