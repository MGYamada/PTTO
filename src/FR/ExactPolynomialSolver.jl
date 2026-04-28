"""
Small exact polynomial-system solver used by the Phase-4 F/R routines.

The systems generated here are intentionally small and usually triangular
after a lexicographic Groebner pass.  The solver chooses a deterministic
gauge by setting unconstrained variables to `1`, then solves the remaining
one-variable equations by factorisation over the requested cyclotomic field.
"""

using Oscar

function _default_context_from_kwargs(; context = nothing, conductor = nothing, N = nothing)
    context !== nothing && return context
    n = conductor === nothing ? N : conductor
    n === nothing && error("a CyclotomicContext or conductor N is required for exact Phase 4")
    return CyclotomicContext(n)
end

function _poly_var_degrees(f)
    ds = zeros(Int, nvars(parent(f)))
    iszero(f) && return ds
    for m in monomials(f)
        md = degrees(m)
        for i in eachindex(md)
            ds[i] = max(ds[i], md[i])
        end
    end
    return ds
end

function _lift_qq_poly_to_field(f, vars_K, K)
    iszero(f) && return zero(vars_K[1])
    out = zero(vars_K[1])
    for (c, m) in zip(coefficients(f), monomials(f))
        term = K(c)
        for (i, d) in enumerate(degrees(m))
            d > 0 && (term *= vars_K[i]^d)
        end
        out += term
    end
    return out
end

function _coerce_poly_to_field_ring(f, vars_K, K)
    base_ring(parent(f)) === K && return f
    return _lift_qq_poly_to_field(f, vars_K, K)
end

function _substitute_fixed_one_polys(eqs, n::Int, fixed_indices::Vector{Int};
                                     var_prefix::Symbol = :x)
    isempty(fixed_indices) && return (eqs = eqs, n = n, free_indices = collect(1:n))
    fixed = Set(fixed_indices)
    free_indices = [i for i in 1:n if !(i in fixed)]
    isempty(eqs) && return (eqs = eqs, n = length(free_indices), free_indices = free_indices)

    K = base_ring(parent(eqs[1]))
    R_new, vars_new = polynomial_ring(K, length(free_indices), var_prefix)
    free_pos = Dict(idx => pos for (pos, idx) in enumerate(free_indices))

    function subst_poly(f)
        iszero(f) && return zero(R_new)
        out = zero(R_new)
        for (c, m) in zip(coefficients(f), monomials(f))
            term = R_new(c)
            for (old_idx, d) in enumerate(degrees(m))
                d == 0 && continue
                if haskey(free_pos, old_idx)
                    term *= vars_new[free_pos[old_idx]]^d
                end
            end
            out += term
        end
        return out
    end

    reduced = filter(!iszero, unique([subst_poly(eq) for eq in eqs]))
    return (eqs = reduced, n = length(free_indices), free_indices = free_indices)
end

function _expand_fixed_one_solution(sol, n::Int, free_indices::Vector{Int}, K)
    full = [one(K) for _ in 1:n]
    length(sol) == length(free_indices) ||
        error("reduced solution length $(length(sol)) != free variable count $(length(free_indices))")
    for (pos, old_idx) in enumerate(free_indices)
        full[old_idx] = sol[pos]
    end
    return full
end

function _eval_exact_poly(f, sol::Vector)
    K = parent(sol[1])
    v = zero(K)
    for (c, m) in zip(coefficients(f), monomials(f))
        term = K(c)
        for (i, d) in enumerate(degrees(m))
            d > 0 && (term *= sol[i]^d)
        end
        v += term
    end
    return v
end

function _partially_univariate_coeffs(f, varidx::Int,
                                      assigned::Dict{Int, Any}, K)
    coeffs = Dict{Int, Any}()
    for (c, m) in zip(coefficients(f), monomials(f))
        degs = degrees(m)
        term = K(c)
        pow = degs[varidx]
        for i in eachindex(degs)
            i == varidx && continue
            d = degs[i]
            d == 0 && continue
            if !haskey(assigned, i)
                return nothing
            end
            term *= assigned[i]^d
        end
        coeffs[pow] = get(coeffs, pow, zero(K)) + term
    end
    return Dict(k => v for (k, v) in coeffs if !iszero(v))
end

function _nf_univar_roots(coeffs::AbstractDict{Int}, K)
    isempty(coeffs) && return nothing
    maxpow = maximum(keys(coeffs))
    maxpow == 0 && return iszero(get(coeffs, 0, zero(K))) ? nothing : Any[]

    U, t = polynomial_ring(K, 1, :t)
    g = zero(t[1])
    for (pow, c) in coeffs
        g += c * t[1]^pow
    end

    roots = Any[]
    for (fac, _) in collect(Oscar.factor(g))
        degree(fac, 1) == 1 || continue
        a = zero(K)
        b = zero(K)
        for (c, m) in zip(coefficients(fac), monomials(fac))
            d = degree(m, 1)
            d == 1 && (a += c)
            d == 0 && (b += c)
        end
        iszero(a) && continue
        push!(roots, -b // a)
    end
    isempty(roots) && return Any[]
    return unique(roots)
end

function _candidate_values_for_var(polys, varidx::Int,
                                   assigned::Dict{Int, Any}, K)
    candidates = nothing
    constrained = false
    for f in polys
        coeffs = _partially_univariate_coeffs(f, varidx, assigned, K)
        coeffs === nothing && continue
        has_var = any(pow -> pow > 0, keys(coeffs))
        has_var || continue
        roots = _nf_univar_roots(coeffs, K)
        roots === nothing && continue
        constrained = true
        isempty(roots) && return Any[]
        candidates = candidates === nothing ? roots : [x for x in candidates if any(==(x), roots)]
        isempty(candidates) && return Any[]
    end
    constrained || return Any[one(K)]
    return candidates
end

function _modular_point_equations(gb_data)
    return :equations in keys(gb_data) ? gb_data.equations : collect(gb_data.gb)
end

function _modular_candidate_values_for_var(polys, varidx::Int,
                                           assigned::Dict{Int, Any}, F, p::Int)
    candidates = nothing
    constrained = false
    for f in polys
        coeffs = _fp_partial_univariate_coeffs(f, varidx, assigned, F)
        coeffs === nothing && continue
        any(pow -> pow > 0, keys(coeffs)) || continue
        roots = _fp_roots_from_coeffs(coeffs, F, p)
        roots === nothing && continue
        constrained = true
        isempty(roots) && return Any[]
        candidates = candidates === nothing ? roots : [x for x in candidates if any(==(x), roots)]
        isempty(candidates) && return Any[]
    end
    constrained || return Any[one(F)]
    return candidates
end

function _enumerate_modular_triangular_solutions(gb_data;
                                                 max_points::Int = 1024)
    max_points >= 1 || error("max_points must be positive")
    polys = collect(gb_data.gb)
    any(isone, polys) && return (complete = true, points = Vector{Int}[])

    equations = _modular_point_equations(gb_data)
    F = base_ring(gb_data.ring)
    p = gb_data.p
    n = length(gb_data.vars)
    points = Vector{Int}[]
    truncated = Ref(false)

    function descend(varidx::Int, assigned::Dict{Int, Any})
        if length(points) >= max_points
            truncated[] = true
            return
        end
        if varidx == 0
            vals = [assigned[i] for i in 1:n]
            all(iszero(_eval_fp_poly(f, vals)) for f in equations) || return
            push!(points, [_fp_elem_to_int(v, p) for v in vals])
            return
        end
        vals = _modular_candidate_values_for_var(polys, varidx, assigned, F, p)
        for v in vals
            truncated[] && return
            next_assigned = copy(assigned)
            next_assigned[varidx] = v
            descend(varidx - 1, next_assigned)
        end
    end

    descend(n, Dict{Int, Any}())
    return (complete = !truncated[], points = unique(points))
end

function _filter_exact_solutions_by_modular_points(sols, modular_points,
                                                   ctx::CyclotomicContext,
                                                   p::Int)
    isempty(sols) && return sols
    zeta_Fp = find_zeta_in_Fp(ctx.N, p)
    allowed = Set(modular_points)
    filtered = Vector{eltype(sols)}()
    for sol in sols
        try
            fp_sol = [cyclotomic_to_Fp(x, zeta_Fp, p) for x in sol]
            fp_sol in allowed && push!(filtered, sol)
        catch err
            @warn "exact solution reduction failed at p=$p; skipping modular filter" exception = (err, catch_backtrace())
            return sols
        end
    end
    return filtered
end

function _estimate_modular_tuple_count(modular_filters)
    isempty(modular_filters) && return 0
    total = BigInt(1)
    for filter in modular_filters
        total *= length(filter.points)
    end
    return total
end

function _solve_via_modular_crt(eqs, n::Int, ctx::CyclotomicContext,
                                modular_filters;
                                max_solutions::Int = 32,
                                reconstruction_bound::Int = 4,
                                denominator_bound::Int = 4,
                                max_crt_tuples::Int = 4096,
                                max_ambiguous_crt_coords::Int = 4096)
    n == 0 && return [elem_type(field(ctx))[]]
    length(modular_filters) >= 2 || return nothing
    all(filter -> !isempty(filter.points), modular_filters) || return Vector{elem_type(field(ctx))}[]
    tuple_count = _estimate_modular_tuple_count(modular_filters)
    tuple_count <= max_crt_tuples || return nothing
    well_determined = crt_reconstruction_is_well_determined(ctx, modular_filters;
                                                            reconstruction_bound = reconstruction_bound,
                                                            denominator_bound = denominator_bound)
    (well_determined || tuple_count * n <= max_ambiguous_crt_coords) || return nothing

    sorted_filters = sort(collect(modular_filters); by = f -> length(f.points))
    point_sets = [filter.points for filter in sorted_filters]
    sols = Vector{elem_type(field(ctx))}[]
    for point_tuple in Iterators.product(point_sets...)
        length(sols) >= max_solutions && break
        sol = elem_type(field(ctx))[]
        ok = true
        for j in 1:n
            residues = Dict(sorted_filters[i].p => point_tuple[i][j]
                            for i in eachindex(sorted_filters))
            x = reconstruct_cyclotomic_element_from_residues(residues, ctx;
                                                             coeff_bound = reconstruction_bound,
                                                             denominator_bound = denominator_bound)
            if x === nothing
                ok = false
                break
            end
            push!(sol, x)
        end
        ok || continue
        R, vars = polynomial_ring(field(ctx), n, :y)
        lifted = [_coerce_poly_to_field_ring(eq, vars, field(ctx)) for eq in eqs]
        all(iszero(_eval_exact_poly(f, sol)) for f in lifted) || continue
        push!(sols, sol)
    end

    return unique(sols)
end

function _solve_exact_via_triangular_groebner(eqs, n::Int, K;
                                              var_prefix::Symbol = :x,
                                              max_solutions::Int = 32)
    max_solutions >= 1 || error("max_solutions must be positive")
    n == 0 && return [elem_type(K)[]]
    R, vars = polynomial_ring(K, n, var_prefix)
    lifted = [_coerce_poly_to_field_ring(eq, vars, K) for eq in eqs]
    lifted = filter(!iszero, lifted)
    isempty(lifted) && return [fill(one(K), n)]

    I = ideal(R, lifted)
    G = collect(groebner_basis(I, ordering = lex(R)))
    any(isone, G) && return Vector{elem_type(K)}[]

    sols = Vector{elem_type(K)}[]
    function descend(varidx::Int, assigned::Dict{Int, Any})
        length(sols) >= max_solutions && return
        if varidx == 0
            sol = [assigned[i] for i in 1:n]
            all(iszero(_eval_exact_poly(f, sol)) for f in lifted) && push!(sols, sol)
            return
        end
        vals = _candidate_values_for_var(G, varidx, assigned, K)
        for v in vals
            next_assigned = copy(assigned)
            next_assigned[varidx] = v
            descend(varidx - 1, next_assigned)
        end
    end
    descend(n, Dict{Int, Any}())

    if isempty(sols)
        error("exact Groebner solver found no cyclotomic solution")
    end
    return unique(sols)
end
