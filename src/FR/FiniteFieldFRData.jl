"""
Finite-field F/R solving boundary.

This file adapts the existing Phase-4 equation generators to a direct
prime-field fiber: symbolic pentagon and hexagon equations are combined over
`QQ`, reduced modulo a prime, solved by the existing finite-field Groebner
point enumerator, and wrapped as `FRData{FpElem}`.
"""

const FiniteFieldFRData = FRData{FpElem}
const _FR_MODP_POINT_CACHE = Dict{Tuple{UInt, Int, Int}, Vector{Vector{Int}}}()

struct FRModPSolveFailure
    rules::FusionRule
    p::Int
    status::Symbol
    message::String
    metadata::Dict{Symbol, Any}
end

function Base.show(io::IO, r::FRModPSolveFailure)
    print(io, "FRModPSolveFailure(p=$(r.p), status=$(r.status), message=$(repr(r.message)))")
end

_fp_values(xs::AbstractVector{<:Integer}, p::Int) = FpElem[FpElem(x, p) for x in xs]

function _known_fr_modp_conductor_and_labels(fr::FusionRule)
    fusion_isomorphic(fr, semion_fusion_rules()) &&
        return (name = :semion, conductor = 8)
    fusion_isomorphic(fr, fibonacci_fusion_rules()) &&
        return (name = :fibonacci, conductor = 20)
    fusion_isomorphic(fr, ising_fusion_rules()) &&
        return (name = :ising, conductor = 16)
    return (name = :custom, conductor = nothing)
end

function _lift_qq_poly_to_vars(f, vars)
    R = parent(vars[1])
    iszero(f) && return zero(R)
    out = zero(R)
    for (c, m) in zip(coefficients(f), monomials(f))
        term = R(c)
        for (i, d) in enumerate(degrees(m))
            d > 0 && (term *= vars[i]^d)
        end
        out += term
    end
    return out
end

function _combined_fr_equations_QQ(fr::FusionRule)
    _, pentagon_eqs, nF = get_pentagon_system(fr.N, fr.rank)
    _, hexagon_eqs, nH = get_hexagon_fr_system(fr.N, fr.rank)
    H, vars = polynomial_ring(QQ, nH, :h)
    pentagon = [_lift_qq_poly_to_vars(eq, vars[1:nF]) for eq in pentagon_eqs]
    hexagon = [_lift_qq_poly_to_vars(eq, vars) for eq in hexagon_eqs]
    return H, filter(e -> !iszero(e), unique(vcat(pentagon, hexagon))), nF, nH
end

function _fr_modp_groebner_points(fr::FusionRule, p::Int; max_points::Int = 4096)
    _, eqs, _, nH = _combined_fr_equations_QQ(fr)
    key = (hash(string(eqs)), nH, p)
    haskey(_FR_MODP_POINT_CACHE, key) && return _FR_MODP_POINT_CACHE[key]
    data = _compute_modular_groebner_data(eqs, nH, [p])
    isempty(data) && error("finite-field Groebner preprocessing failed at p=$p")
    points = _enumerate_modular_triangular_solutions(first(data); max_points = max_points)
    points.complete || error("finite-field solution enumeration at p=$p exceeded max_points=$max_points")
    _FR_MODP_POINT_CACHE[key] = points.points
    return points.points
end

function _modp_phase4_reference_solution(fr::FusionRule, p::Int, info;
                                         primes::Vector{Int} = [p],
                                         kwargs...)
    info.conductor === nothing &&
        error("direct finite-field solve did not finish and no conductor metadata is known")
    data = info.name == :semion ? semion_modular_data(CyclotomicContext(info.conductor)) :
           info.name == :fibonacci ? fibonacci_modular_data(CyclotomicContext(info.conductor)) :
           info.name == :ising ? ising_modular_data(CyclotomicContext(info.conductor)) :
           error("no Phase-4 branch metadata is available for $(info.name)")
    twists = [data.T[i, i] for i in 1:length(data.labels)]
    result = compute_FR_from_ST(fr.N;
                                conductor = info.conductor,
                                primes = primes,
                                S = data.S,
                                T = twists,
                                kwargs...)
    Fp = [reduce_mod_p(data.context, x, p) for x in result.F]
    Rp = [reduce_mod_p(data.context, x, p) for x in result.R]
    _, nR = _braiding_block_positions(fr.N)
    length(Rp) == 2 * nR || error("Phase-4 result returned $(length(Rp)) R coordinates; expected $(2 * nR)")
    return vcat(Fp, Rp)
end

function _select_modp_branch(points, fr::FusionRule, p::Int, branch, info)
    branch isa Integer && return Int(branch)
    branch == :auto || return 1
    nF = get(fr_equation_system(fr).metadata, :f_variables, 0)
    if info.name == :semion
        idx = findfirst(pt -> nF >= 1 && mod(pt[1], p) != 1, points)
        idx !== nothing && return idx
    end
    return 1
end

function _eval_qq_poly_mod_p(f, vals::AbstractVector{<:Integer}, p::Int)
    total = 0
    for (c, m) in zip(coefficients(f), monomials(f))
        term = _rational_mod(c, p)
        for (i, d) in enumerate(degrees(m))
            d > 0 && (term = mod(term * powermod(mod(vals[i], p), d, p), p))
        end
        total = mod(total + term, p)
    end
    return total
end

function _verify_pentagon_modp(fr::FusionRule, F::AbstractVector{<:Integer}, p::Int)
    _, eqs, nF = get_pentagon_system(fr.N, fr.rank)
    length(F) == nF || return false
    return all(eq -> _eval_qq_poly_mod_p(eq, F, p) == 0, eqs)
end

function _verify_hexagon_modp(fr::FusionRule, values::AbstractVector{<:Integer}, p::Int)
    _, eqs, nH = get_hexagon_fr_system(fr.N, fr.rank)
    length(values) == nH || return false
    return all(eq -> _eval_qq_poly_mod_p(eq, values, p) == 0, eqs)
end

function frdata_from_modp_solution(rules, solution::AbstractVector{<:Integer}, p::Integer;
                                   metadata = Dict{Symbol, Any}())
    fr = _fusion_rule(rules)
    pp = Int(p)
    _, _, nF = get_pentagon_system(fr.N, fr.rank)
    _, nR = _braiding_block_positions(fr.N)
    length(solution) == nF + 2 * nR ||
        error("mod-p F/R solution has length $(length(solution)); expected $(nF + 2 * nR)")
    F = _fp_values(solution[1:nF], pp)
    R = _fp_values(solution[(nF + 1):(nF + nR)], pp)
    Rinv = _fp_values(solution[(nF + nR + 1):(nF + 2nR)], pp)
    info = _known_fr_modp_conductor_and_labels(fr)
    meta = copy(metadata)
    meta[:base_field] = Symbol("F_$pp")
    meta[:p] = pp
    meta[:prime] = pp
    meta[:solver_status] = get(meta, :solver_status, :solved)
    meta[:verification] = get(meta, :verification,
                              Dict(:pentagon => _verify_pentagon_modp(fr, solution[1:nF], pp),
                                   :hexagon => _verify_hexagon_modp(fr, solution, pp)))
    meta[:simple_objects] = get(meta, :simple_objects, simple_objects(fr))
    meta[:name] = get(meta, :name, info.name)
    info.conductor === nothing || (meta[:conductor] = get(meta, :conductor, info.conductor))
    meta[:source] = get(meta, :source, :phase4_modp_solution)
    meta[:format] = :tensorcategories_variable_order
    return FRData(fr, F, R, Rinv; metadata = meta)
end

"""
    solve_fr_mod_p(rules, p; strategy=:safe, branch=:auto, max_points=4096, failure=:throw)

Solve the multiplicity-free pentagon/hexagon equations over `F_p` and return
`FRData{FpElem}`.  The implementation uses the existing Phase-4 equation
generators and finite-field Groebner point enumeration; no hard-coded F/R
coordinates are used.
"""
function solve_fr_mod_p(rules, p::Integer; strategy::Symbol = :safe,
                        branch = :auto, max_points::Int = 4096,
                        solver::Symbol = :auto,
                        failure::Symbol = :throw, kwargs...)
    isempty(kwargs) || error("unsupported keyword arguments: $(collect(keys(kwargs)))")
    pp = Int(p)
    pp > 1 && isprime(pp) || error("p must be prime, got $p")
    fr = _fusion_rule(rules)
    require_multiplicity_free(fr)
    try
        info = _known_fr_modp_conductor_and_labels(fr)
        fixed = gauge_fix(fr_equation_system(fr); strategy = strategy)
        reduced = reduce_mod_p(fixed, pp)
        points = if solver == :phase4
            [_modp_phase4_reference_solution(fr, pp, info; kwargs...)]
        elseif solver == :direct || solver == :auto
            isempty(kwargs) || error("unsupported keyword arguments for direct solver: $(collect(keys(kwargs)))")
            _fr_modp_groebner_points(fr, pp; max_points = max_points)
        else
            error("unknown solver=$(repr(solver)); expected :auto, :direct, or :phase4")
        end
        isempty(points) && error("no F/R solution found over F_$pp")
        selected_index = _select_modp_branch(points, fr, pp, branch, info)
        1 <= selected_index <= length(points) ||
            error("branch index $selected_index outside 1:$(length(points))")
        verification = Dict(:pentagon => _verify_pentagon_modp(fr, points[selected_index][1:fixed.metadata[:f_variables]], pp),
                            :hexagon => _verify_hexagon_modp(fr, points[selected_index], pp))
        metadata = Dict{Symbol, Any}(:solver_status => :solved,
                                     :solver => solver,
                                     :gauge_fixing => fixed.metadata,
                                     :finite_field_system => reduced.metadata,
                                     :branch => branch,
                                     :branch_index => selected_index,
                                     :solution_count => length(points),
                                     :verification => verification)
        return frdata_from_modp_solution(fr, points[selected_index], pp; metadata = metadata)
    catch err
        failure == :return || rethrow(err)
        return FRModPSolveFailure(fr, pp, :failed, sprint(showerror, err),
                                  Dict{Symbol, Any}(:strategy => strategy, :branch => branch))
    end
end

semion_fr_data_mod_p(p::Integer; kwargs...) =
    solve_fr_mod_p(semion_fusion_rules(), p; kwargs...)
fibonacci_fr_data_mod_p(p::Integer; kwargs...) =
    solve_fr_mod_p(fibonacci_fusion_rules(), p; kwargs...)
ising_fr_data_mod_p(p::Integer; kwargs...) =
    solve_fr_mod_p(ising_fusion_rules(), p; kwargs...)

function verify_pentagon(data::FRData{FpElem})
    return _verify_pentagon_modp(fusion_rule(data),
                                 [x.value for x in F_values(data)],
                                 first(F_values(data)).p)
end

function verify_hexagon(data::FRData{FpElem})
    p = !isempty(F_values(data)) ? first(F_values(data)).p : first(R_values(data)).p
    vals = [x.value for x in vcat(F_values(data), R_values(data), R_inverse_values(data))]
    return _verify_hexagon_modp(fusion_rule(data), vals, p)
end

verify_FRData(data::FRData{FpElem}) = verify_pentagon(data) && verify_hexagon(data)
