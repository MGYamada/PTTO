"""
Gauge stabilizer and automorphism infrastructure.

Experimental API introduced in ACMG.jl v0.9.1.  This interface may change
before v1.0.
"""

"""
    StabilizerProblem(solution, gauge_group; metadata = NamedTuple())

Experimental API introduced in ACMG.jl v0.9.1.  This interface may change
before v1.0.

Describe the stabilizer problem `g ⋅ solution == solution` for a gauge group.
The gauge group may be an explicit finite collection of gauge elements, or an
ACMG finite-field gauge group returned by `ACMG.finite_field_gauge_group`.
"""
struct StabilizerProblem{S,G,M}
    solution::S
    gauge_group::G
    metadata::M
end

StabilizerProblem(solution, gauge_group; metadata = NamedTuple()) =
    StabilizerProblem{typeof(solution), typeof(gauge_group), typeof(metadata)}(
        solution, gauge_group, metadata)

"""
    StabilizerEquations(equations; metadata = NamedTuple())

Experimental API introduced in ACMG.jl v0.9.1.  This interface may change
before v1.0.

Container for equations defining `g ⋅ x == x`.  In the current finite-field
toric gauge implementation, equations are stored as character equations on
the nonzero F/R symbol coordinates.
"""
struct StabilizerEquations{E,M}
    equations::E
    metadata::M
end

StabilizerEquations(equations; metadata = NamedTuple()) =
    StabilizerEquations{typeof(equations), typeof(metadata)}(equations, metadata)

"""
    StabilizerResult(automorphisms, order; metadata = NamedTuple())

Experimental API introduced in ACMG.jl v0.9.1.  This interface may change
before v1.0.

Result of a stabilizer computation.  `automorphisms` is either the enumerated
list of gauge transformations fixing the point, or `nothing` when only the
order was requested.
"""
struct StabilizerResult{A,N,M}
    automorphisms::A
    order::N
    metadata::M
end

StabilizerResult(automorphisms, order; metadata = NamedTuple()) =
    StabilizerResult{typeof(automorphisms), typeof(order), typeof(metadata)}(
        automorphisms, order, metadata)

"""
    automorphisms(result)

Return the automorphisms stored in a `StabilizerResult`.
Experimental API introduced in ACMG.jl v0.9.1.  This interface may change
before v1.0.
"""
automorphisms(result::StabilizerResult) = result.automorphisms

"""
    stabilizer_order(result)
    stabilizer_order(solution, gauge_group; kwargs...)

Return `|Aut(x)|` for a stabilizer result or compute it from a stabilizer
problem.  Experimental API introduced in ACMG.jl v0.9.1.
"""
stabilizer_order(result::StabilizerResult) = result.order

"""
    is_trivial_stabilizer(result)

Return `true` when `stabilizer_order(result) == 1`.
Experimental API introduced in ACMG.jl v0.9.1.  This interface may change
before v1.0.
"""
is_trivial_stabilizer(result::StabilizerResult) = stabilizer_order(result) == 1

"""
    stabilizer_metadata(result)

Return metadata recorded by a `StabilizerResult`.
Experimental API introduced in ACMG.jl v0.9.1.  This interface may change
before v1.0.
"""
stabilizer_metadata(result::StabilizerResult) = result.metadata

struct _FiniteFieldGaugeGroup{M}
    p::Int
    parameters::Vector{GaugeParameter}
    metadata::M
end

function _metadata_get(metadata, key::Symbol, default)
    metadata isa AbstractDict && return get(metadata, key, default)
    metadata isa NamedTuple && return haskey(metadata, key) ? getproperty(metadata, key) : default
    return default
end

function _metadata_with(metadata, pairs::Pair{Symbol}...)
    if metadata isa AbstractDict
        out = copy(metadata)
        for (k, v) in pairs
            out[k] = v
        end
        return out
    elseif metadata isa NamedTuple
        return merge(metadata, NamedTuple(pairs))
    end
    out = Dict{Symbol, Any}(:user_metadata => metadata)
    for (k, v) in pairs
        out[k] = v
    end
    return out
end

function _frdata_prime_for_stabilizer(data::FRData{FpElem})
    p = _toric_prime(data)
    p === nothing && error("cannot infer finite field prime for FRData")
    return Int(p)
end

"""
    ACMG.finite_field_gauge_group(frdata::FRData{FpElem})

Internal experimental helper for constructing the finite toric gauge group
`(F_p^*)^n` attached to multiplicity-free finite-field `FRData`.
"""
function finite_field_gauge_group(data::FRData{FpElem}; metadata = NamedTuple())
    validate_frdata_for_gauge(data)
    p = _frdata_prime_for_stabilizer(data)
    params = gauge_parameters(data)
    meta = _metadata_with(metadata, :kind => :finite_field_toric_gauge_group,
                          :p => p, :parameter_count => length(params),
                          :parameters => params)
    return _FiniteFieldGaugeGroup{typeof(meta)}(p, params, meta)
end

_finite_gauge_group_order(group::_FiniteFieldGaugeGroup) =
    BigInt(group.p - 1)^length(group.parameters)

function _gauge_action_from_tuple(group::_FiniteFieldGaugeGroup, values)
    params = Dict{Tuple{Int,Int,Int,Int}, FpElem}()
    for (ch, value) in zip(group.parameters, values)
        params[(ch[1], ch[2], ch[3], 1)] = FpElem(value, group.p)
    end
    return GaugeAction(params; field = Symbol("F_", group.p),
                       metadata = Dict{Symbol, Any}(:kind => :finite_field_gauge_element,
                                                    :p => group.p))
end

function _iter_finite_field_gauge_actions(group::_FiniteFieldGaugeGroup)
    values = collect(1:(group.p - 1))
    isempty(group.parameters) &&
        return (_gauge_action_from_tuple(group, ()) for _ in 1:1)
    products = Iterators.product(ntuple(_ -> values, length(group.parameters))...)
    return (_gauge_action_from_tuple(group, tuple_values) for tuple_values in products)
end

function _default_stabilizer_action(solution, gauge)
    if solution isa FRData && gauge isa GaugeAction
        return apply_gauge(solution, gauge)
    end
    error("no stabilizer action provided for solution $(typeof(solution)) and gauge element $(typeof(gauge)); pass action=(x,g)->...")
end

function _stabilizer_action(problem::StabilizerProblem, action)
    action !== nothing && return action
    stored = _metadata_get(problem.metadata, :action, nothing)
    stored !== nothing && return stored
    return _default_stabilizer_action
end

function _solution_equal(a::FRData, b::FRData)
    return fusion_rule(a).N == fusion_rule(b).N &&
           F_values(a) == F_values(b) &&
           R_values(a) == R_values(b) &&
           R_inverse_values(a) == R_inverse_values(b)
end

_solution_equal(a, b) = a == b

function _stabilizer_metadata(problem, method, order, return_automorphisms, extra = NamedTuple())
    meta = (method = method,
            order = order,
            return_automorphisms = return_automorphisms,
            problem_metadata = problem.metadata)
    return merge(meta, extra)
end

function _enumerated_stabilizer(problem::StabilizerProblem, elements;
                                action = nothing,
                                return_automorphisms::Bool = true)
    act = _stabilizer_action(problem, action)
    autos = Any[]
    count = 0
    for g in elements
        _solution_equal(act(problem.solution, g), problem.solution) || continue
        count += 1
        return_automorphisms && push!(autos, g)
    end
    order = BigInt(count)
    meta = _stabilizer_metadata(problem, :bruteforce, order, return_automorphisms)
    return StabilizerResult(return_automorphisms ? autos : nothing, order; metadata = meta)
end

function _finite_field_gauge_stabilizer(problem::StabilizerProblem{<:FRData,<:_FiniteFieldGaugeGroup};
                                        return_automorphisms::Bool = true,
                                        max_enumeration::Integer = 100_000,
                                        action = nothing)
    group = problem.gauge_group
    group_order = _finite_gauge_group_order(group)
    if group_order > BigInt(max_enumeration)
        error("finite gauge group has order $group_order, exceeding max_enumeration=$max_enumeration; call stabilizer_equations(problem) for algebraic data")
    end
    act = _stabilizer_action(problem, action)
    autos = GaugeAction[]
    count = 0
    for g in _iter_finite_field_gauge_actions(group)
        _solution_equal(act(problem.solution, g), problem.solution) || continue
        count += 1
        return_automorphisms && push!(autos, g)
    end
    order = BigInt(count)
    meta = _stabilizer_metadata(problem, :bruteforce_finite_field_gauge,
                                order, return_automorphisms,
                                (field = Symbol("F_", group.p),
                                 group_order = group_order,
                                 parameter_count = length(group.parameters)))
    return StabilizerResult(return_automorphisms ? autos : nothing, order; metadata = meta)
end

"""
    stabilizer(problem::StabilizerProblem; method = :auto, return_automorphisms = true)
    stabilizer(solution, gauge_group; kwargs...)

Compute the stabilizer `Aut(x) = {g | g ⋅ x = x}` when the gauge group is
small enough to enumerate.  Experimental API introduced in ACMG.jl v0.9.1.

For unsupported or too-large groups, use `stabilizer_equations(problem)` to
obtain algebraic stabilizer equations.
"""
function stabilizer(problem::StabilizerProblem; method::Symbol = :auto,
                    field = nothing,
                    return_automorphisms::Bool = true,
                    max_enumeration::Integer = 100_000,
                    action = nothing)
    method in (:auto, :bruteforce) ||
        error("unsupported stabilizer method $(repr(method)); supported methods are :auto and :bruteforce")
    problem.gauge_group === nothing &&
        error("stabilizer requires a gauge group or finite collection of gauge elements")
    problem.gauge_group isa _FiniteFieldGaugeGroup &&
        return _finite_field_gauge_stabilizer(problem;
                                             return_automorphisms = return_automorphisms,
                                             max_enumeration = max_enumeration,
                                             action = action)
    if applicable(iterate, problem.gauge_group)
        return _enumerated_stabilizer(problem, problem.gauge_group;
                                      action = action,
                                      return_automorphisms = return_automorphisms)
    end
    error("gauge group $(typeof(problem.gauge_group)) is not enumerable; call stabilizer_equations(problem) if equations are available")
end

stabilizer(solution, gauge_group; metadata = NamedTuple(), kwargs...) =
    stabilizer(StabilizerProblem(solution, gauge_group; metadata = metadata); kwargs...)

stabilizer_order(solution, gauge_group; kwargs...) =
    stabilizer_order(stabilizer(solution, gauge_group; return_automorphisms = false,
                                kwargs...))

function _toric_stabilizer_equations(data::FRData, group::_FiniteFieldGaugeGroup)
    symbol_data = _toric_symbol_data(data; include_R = true)
    coords = _symbol_data_active_coordinates(symbol_data)
    equations = NamedTuple[]
    for coord in coords
        weight = _coordinate_weight(coord)
        push!(equations, (coordinate = coord, weight = weight, rhs = FpElem(1, group.p)))
    end
    meta = (kind = :toric_character_equations,
            field = Symbol("F_", group.p),
            parameters = group.parameters,
            coordinate_count = length(coords))
    return StabilizerEquations(equations; metadata = meta)
end

"""
    stabilizer_equations(problem::StabilizerProblem)
    stabilizer_equations(solution, gauge_group)

Return equations defining `g ⋅ x == x` when ACMG knows how to express them.
Experimental API introduced in ACMG.jl v0.9.1.
"""
function stabilizer_equations(problem::StabilizerProblem)
    if problem.solution isa FRData && problem.gauge_group isa _FiniteFieldGaugeGroup
        return _toric_stabilizer_equations(problem.solution, problem.gauge_group)
    end
    stored = _metadata_get(problem.metadata, :equations, nothing)
    stored !== nothing && return StabilizerEquations(stored; metadata = problem.metadata)
    error("stabilizer equations are not implemented for solution $(typeof(problem.solution)) and gauge group $(typeof(problem.gauge_group))")
end

stabilizer_equations(solution, gauge_group; metadata = NamedTuple()) =
    stabilizer_equations(StabilizerProblem(solution, gauge_group; metadata = metadata))
