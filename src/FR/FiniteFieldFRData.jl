"""
Finite-field F/R solving boundary.

This file adapts the existing Phase-4 equation generators to a direct
prime-field fiber: symbolic pentagon and hexagon equations are combined over
`QQ`, reduced modulo a prime, solved by the existing finite-field Groebner
point enumerator, and wrapped as `FRData{FpElem}`.
"""

"""
    FiniteFieldFRData

Alias for `FRData{FpElem}`, the finite-field F/R data container produced by
`solve_fr_mod_p` and `frdata_from_modp_solution`.
"""
const FiniteFieldFRData = FRData{FpElem}
const _FR_MODP_POINT_CACHE = Dict{Tuple{UInt, Int, Int}, Vector{Vector{Int}}}()

"""
    FRModPSolveFailure

Failure record returned by `solve_fr_mod_p(...; failure=:return)`.
"""
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
    # Fusion rules alone do not determine a braided MTC.  These labels identify
    # ACMG's built-in reference branches used for branch selection and Phase-4
    # fallback, not a mathematical classification from fusion data alone.
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
    one_vec = zeros(Int, fr.rank)
    one_vec[1] = 1
    H, hexagon_eqs, hex_nF, nR = _hexagon_equations_with_symbolic_FR(fr.N, one_vec)
    hex_nF == nF ||
        error("TensorCategories F-variable layout mismatch: pentagon has $nF variables but hexagon reserves $hex_nF")
    nH = hex_nF + 2 * nR
    nvars(H) == nH ||
        error("hexagon ring has $(nvars(H)) variables; expected $nH")
    vars = gens(H)
    pentagon = [_lift_qq_poly_to_vars(eq, vars[1:nF]) for eq in pentagon_eqs]
    hexagon = [_lift_qq_poly_to_vars(eq, vars) for eq in hexagon_eqs]
    return H, filter(e -> !iszero(e), unique(vcat(pentagon, hexagon))), nF, nH
end

function _fr_modp_groebner_points(system::FREquationSystem, p::Int; max_points::Int = 4096)
    eqs = system.equations
    nH = Int(get(system.metadata, :hexagon_variables,
                 get(system.metadata, :f_variables, 0)))
    key = (hash(string(eqs)), nH, p)
    haskey(_FR_MODP_POINT_CACHE, key) && return _FR_MODP_POINT_CACHE[key]
    data = _compute_modular_groebner_data(eqs, nH, [p])
    isempty(data) && error("finite-field Groebner preprocessing failed at p=$p")
    points = _enumerate_modular_triangular_solutions(first(data); max_points = max_points)
    points.complete || error("finite-field solution enumeration at p=$p exceeded max_points=$max_points")
    _FR_MODP_POINT_CACHE[key] = points.points
    return points.points
end

function _fr_modp_groebner_points(fr::FusionRule, p::Int; max_points::Int = 4096)
    return _fr_modp_groebner_points(gauge_fix(fr_equation_system(fr)), p;
                                    max_points = max_points)
end

function _expand_fixed_one_modp_solution(point::AbstractVector{<:Integer},
                                         original_n::Int,
                                         free_indices::AbstractVector{<:Integer},
                                         p::Int)
    length(point) == length(free_indices) ||
        error("reduced solution has length $(length(point)); expected $(length(free_indices))")
    full = ones(Int, original_n)
    for (i, idx) in enumerate(free_indices)
        1 <= idx <= original_n || error("free variable index $idx outside 1:$original_n")
        full[Int(idx)] = mod(Int(point[i]), p)
    end
    return full
end

function _expand_gauge_fixed_modp_points(points::Vector{Vector{Int}},
                                         system::FREquationSystem,
                                         p::Int)
    original_n = Int(get(system.metadata, :original_hexagon_variables,
                         get(system.metadata, :hexagon_variables,
                             get(system.metadata, :f_variables, 0))))
    free_indices = get(system.metadata, :free_hexagon_indices, collect(1:original_n))
    return [_expand_fixed_one_modp_solution(point, original_n, free_indices, p)
            for point in points]
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

function _known_fr_modp_modular_data(info)
    info.name == :semion && return semion_modular_data(CyclotomicContext(info.conductor))
    info.name == :fibonacci && return fibonacci_modular_data(CyclotomicContext(info.conductor))
    info.name == :ising && return ising_modular_data(CyclotomicContext(info.conductor))
    return nothing
end

_modp_inv(a::Integer, p::Int) = invmod(mod(Int(a), p), p)

function _reduce_ctx_mod_p(ctx, x, p::Int)
    return mod(Int(reduce_mod_p(ctx, x, p)), p)
end

function _known_target_twists_mod_p(info, p::Int)
    data = _known_fr_modp_modular_data(info)
    data === nothing && return nothing
    twists = [_reduce_ctx_mod_p(data.context, data.T[i, i], p) for i in 1:length(data.labels)]
    isempty(twists) && return twists
    iszero(twists[1]) && return nothing
    t0inv = _modp_inv(twists[1], p)
    return [_mul_mod_int(t, t0inv, p) for t in twists]
end

function _known_quantum_dimensions_mod_p(info, p::Int)
    data = _known_fr_modp_modular_data(info)
    data === nothing && return nothing
    S11 = _reduce_ctx_mod_p(data.context, data.S[1, 1], p)
    iszero(S11) && return nothing
    invS11 = _modp_inv(S11, p)
    return [_mul_mod_int(_reduce_ctx_mod_p(data.context, data.S[a, 1], p), invS11, p)
            for a in 1:length(data.labels)]
end

function _trace_R_block_mod_p(Rfwd::AbstractVector{<:Integer},
                              positions,
                              Nijk::Array{Int,3},
                              a::Int, b::Int, c::Int,
                              p::Int)
    n = Nijk[a, b, c]
    n == 0 && return 0
    pos = positions[(a, b, c)]
    acc = 0
    for μ in 1:n
        acc = mod(acc + Int(Rfwd[pos[(μ - 1) * n + μ]]), p)
    end
    return acc
end

function _twists_from_solution_mod_p(point::AbstractVector{<:Integer},
                                     fr::FusionRule,
                                     p::Int,
                                     info)
    d = _known_quantum_dimensions_mod_p(info, p)
    d === nothing && return nothing
    _, _, nF = get_pentagon_system(fr.N, fr.rank)
    positions, nR = _braiding_block_positions(fr.N)
    length(point) == nF + 2 * nR || return nothing
    Rfwd = point[(nF + 1):(nF + nR)]
    twists = zeros(Int, fr.rank)
    for a in 1:fr.rank
        iszero(d[a]) && return nothing
        acc = 0
        for c in 1:fr.rank
            fr.N[a, a, c] == 0 && continue
            tr = _trace_R_block_mod_p(Rfwd, positions, fr.N, a, a, c, p)
            acc = mod(acc + _mul_mod_int(d[c], tr, p), p)
        end
        twists[a] = _mul_mod_int(acc, _modp_inv(d[a], p), p)
    end
    iszero(twists[1]) && return nothing
    t0inv = _modp_inv(twists[1], p)
    return [_mul_mod_int(t, t0inv, p) for t in twists]
end

function _twists_match_target_mod_p(twists::AbstractVector{<:Integer},
                                    target::AbstractVector{<:Integer},
                                    fr::FusionRule)
    length(twists) == length(target) == fr.rank || return false
    for perm in fusion_automorphisms(fr)
        all(i -> Int(twists[perm[i]]) == Int(target[i]), 1:fr.rank) &&
            return true
    end
    return false
end

function _select_modp_branch(points, fr::FusionRule, p::Int, branch, info)
    branch isa Integer && return Int(branch)
    branch == :auto || return 1
    target_twists = try
        _known_target_twists_mod_p(info, p)
    catch
        nothing
    end
    if target_twists !== nothing
        idx = findfirst(points) do pt
            twists = _twists_from_solution_mod_p(pt, fr, p, info)
            twists !== nothing && _twists_match_target_mod_p(twists, target_twists, fr)
        end
        idx !== nothing && return idx
    end
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
            d > 0 && (term = _mul_mod_int(term, powermod(mod(vals[i], p), d, p), p))
        end
        total = Int(mod(Int128(total) + Int128(term), p))
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

"""
    frdata_from_modp_solution(rules, solution, p; metadata=Dict())

Wrap a finite-field solution vector in TensorCategories F/R variable order as
`FRData{FpElem}`.  The solution must contain all F coordinates followed by
forward R coordinates and explicit inverse-R coordinates.
"""
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
    solve_fr_mod_p(rules, p; strategy=:toric_snf, branch=:auto, max_points=4096, failure=:throw)

Solve the multiplicity-free pentagon/hexagon equations over `F_p` and return
`FRData{FpElem}`.  The implementation uses the existing Phase-4 equation
generators and finite-field Groebner point enumeration; no hard-coded F/R
coordinates are used.
"""
function solve_fr_mod_p(rules, p::Integer; strategy::Symbol = :toric_snf,
                        branch = :auto, max_points::Int = 4096,
                        solver::Symbol = :auto,
                        failure::Symbol = :throw, kwargs...)
    pp = Int(p)
    pp > 1 && isprime(pp) || error("p must be prime, got $p")
    fr = _fusion_rule(rules)
    require_multiplicity_free(fr)
    try
        info = _known_fr_modp_conductor_and_labels(fr)
        fixed = gauge_fix(fr_equation_system(fr); strategy = strategy)
        reduced = reduce_mod_p(fixed, pp)
        reduced_points = if solver == :phase4
            [_modp_phase4_reference_solution(fr, pp, info; kwargs...)]
        elseif solver == :direct || solver == :auto
            isempty(kwargs) || error("unsupported keyword arguments for direct solver: $(collect(keys(kwargs)))")
            _fr_modp_groebner_points(fixed, pp; max_points = max_points)
        else
            error("unknown solver=$(repr(solver)); expected :auto, :direct, or :phase4")
        end
        points = solver == :phase4 ? reduced_points :
            _expand_gauge_fixed_modp_points(reduced_points, fixed, pp)
        isempty(points) && error("no F/R solution found over F_$pp")
        selected_index = _select_modp_branch(points, fr, pp, branch, info)
        1 <= selected_index <= length(points) ||
            error("branch index $selected_index outside 1:$(length(points))")
        nF = Int(get(fixed.metadata, :original_f_variables,
                     get(fr_equation_system(fr).metadata, :f_variables, 0)))
        verification = Dict(:pentagon => _verify_pentagon_modp(fr, points[selected_index][1:nF], pp),
                            :hexagon => _verify_hexagon_modp(fr, points[selected_index], pp))
        metadata = Dict{Symbol, Any}(:solver_status => :solved,
                                     :solver => solver,
                                     :gauge_fixing => fixed.metadata,
                                     :finite_field_system => reduced.metadata,
                                     :branch => branch,
                                     :branch_index => selected_index,
                                     :solution_count => length(points),
                                     :reduced_solution_count => length(reduced_points),
                                     :verification => verification)
        return frdata_from_modp_solution(fr, points[selected_index], pp; metadata = metadata)
    catch err
        failure == :return || rethrow(err)
        return FRModPSolveFailure(fr, pp, :failed, sprint(showerror, err),
                                  Dict{Symbol, Any}(:strategy => strategy, :branch => branch))
    end
end

"""
    semion_fr_data_mod_p(p; kwargs...)

Solve the semion F/R equations over `F_p` and return `FRData{FpElem}`.
"""
semion_fr_data_mod_p(p::Integer; kwargs...) =
    solve_fr_mod_p(semion_fusion_rules(), p; kwargs...)

"""
    fibonacci_fr_data_mod_p(p; kwargs...)

Solve the Fibonacci F/R equations over `F_p` and return `FRData{FpElem}`.
"""
fibonacci_fr_data_mod_p(p::Integer; kwargs...) =
    solve_fr_mod_p(fibonacci_fusion_rules(), p; kwargs...)

"""
    ising_fr_data_mod_p(p; kwargs...)

Solve the Ising F/R equations over `F_p` and return `FRData{FpElem}`.
"""
ising_fr_data_mod_p(p::Integer; kwargs...) =
    solve_fr_mod_p(ising_fusion_rules(), p; kwargs...)

function _frdata_prime(data::FRData{FpElem})
    !isempty(F_values(data)) && return first(F_values(data)).p
    !isempty(R_values(data)) && return first(R_values(data)).p
    !isempty(R_inverse_values(data)) && return first(R_inverse_values(data)).p
    haskey(fr_metadata(data), :p) && return Int(fr_metadata(data)[:p])
    haskey(fr_metadata(data), :prime) && return Int(fr_metadata(data)[:prime])
    error("cannot infer prime from empty FRData; include metadata[:p]")
end

"""
    verify_pentagon(frdata)

Evaluate the pentagon equations at finite-field F-symbol coordinates.
"""
function verify_pentagon(data::FRData{FpElem})
    return _verify_pentagon_modp(fusion_rule(data),
                                 [x.value for x in F_values(data)],
                                 _frdata_prime(data))
end

"""
    verify_hexagon(frdata)

Evaluate the hexagon equations at finite-field F/R-symbol coordinates.
"""
function verify_hexagon(data::FRData{FpElem})
    p = _frdata_prime(data)
    vals = [x.value for x in vcat(F_values(data), R_values(data), R_inverse_values(data))]
    return _verify_hexagon_modp(fusion_rule(data), vals, p)
end

"""
    verify_FRData(frdata)

Return whether finite-field `FRData` satisfies both pentagon and hexagon
equations.
"""
verify_FRData(data::FRData{FpElem}) = verify_pentagon(data) && verify_hexagon(data)
