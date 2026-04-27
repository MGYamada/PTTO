"""
Small exact polynomial-system solver used by the Phase-4 F/R routines.

The systems generated here are intentionally small and usually triangular
after a lexicographic Groebner pass.  The solver chooses a deterministic
gauge by setting unconstrained variables to `1`, then solves the remaining
one-variable equations by factorisation over the requested cyclotomic field.
"""

using Oscar

const _FR_COEFF_VECTOR_CACHE = Dict{Tuple{Int, Int}, Vector{Vector{Int}}}()
const _FR_MITM_LEFT_CACHE = Dict{Tuple{Int, Int, Tuple{Vararg{Int}}, Tuple{Vararg{Int}}}, Dict{Tuple{Vararg{Int}}, Vector{Int}}}()
const _FR_MITM_RIGHT_CACHE = Dict{Tuple{Int, Int, Int, Tuple{Vararg{Int}}, Tuple{Vararg{Int}}}, Vector{Tuple{Vector{Int}, Vector{Int}}}}()

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

function _eval_power_basis_mod(coeffs::AbstractVector{<:Integer},
                               zeta_Fp::Int,
                               p::Int)
    acc = 0
    zpow = 1
    for c in coeffs
        acc = mod(acc + mod(c, p) * zpow, p)
        zpow = mod(zpow * zeta_Fp, p)
    end
    return acc
end

function _coeff_vectors(len::Int, bound::Int)
    key = (len, bound)
    haskey(_FR_COEFF_VECTOR_CACHE, key) && return _FR_COEFF_VECTOR_CACHE[key]
    len == 0 && return [Int[]]
    ranges = ntuple(_ -> -bound:bound, len)
    vectors = vec([collect(t) for t in Iterators.product(ranges...)])
    sort!(vectors; lt = (a, b) -> isless(_coeff_score(a), _coeff_score(b)))
    _FR_COEFF_VECTOR_CACHE[key] = vectors
    return vectors
end

function _coeff_score(v::Vector{Int})
    abs_vals = [abs(x) for x in v]
    max_abs = isempty(abs_vals) ? 0 : maximum(abs_vals)
    return (sum(abs_vals), max_abs, v)
end

function _reconstruct_power_basis_coeffs_mitm(targets::Vector{Int},
                                              zetas::Vector{Int},
                                              primes::Vector{Int},
                                              degree_K::Int,
                                              coeff_bound::Int)
    left_len = degree_K ÷ 2
    right_len = degree_K - left_len
    left_vectors = _coeff_vectors(left_len, coeff_bound)
    right_vectors = _coeff_vectors(right_len, coeff_bound)

    zetas_key = Tuple(zetas)
    primes_key = Tuple(primes)
    left_key = (left_len, coeff_bound, zetas_key, primes_key)
    table = get!(_FR_MITM_LEFT_CACHE, left_key) do
        built = Dict{Tuple{Vararg{Int}}, Vector{Int}}()
        for c_left in left_vectors
            key = Tuple(_eval_power_basis_mod(c_left, zetas[i], primes[i])
                        for i in eachindex(primes))
            if !haskey(built, key) || _coeff_score(c_left) < _coeff_score(built[key])
                built[key] = c_left
            end
        end
        built
    end

    right_key = (right_len, coeff_bound, left_len, zetas_key, primes_key)
    right_data = get!(_FR_MITM_RIGHT_CACHE, right_key) do
        zeta_offsets = [powermod(zetas[i], left_len, primes[i]) for i in eachindex(primes)]
        data = Tuple{Vector{Int}, Vector{Int}}[]
        for c_right in right_vectors
            vals = Int[]
            for i in eachindex(primes)
                right_val = _eval_power_basis_mod(c_right, zetas[i], primes[i])
                push!(vals, mod(zeta_offsets[i] * right_val, primes[i]))
            end
            push!(data, (vals, c_right))
        end
        data
    end

    best = nothing
    best_score = nothing
    for (vals, c_right) in right_data
        need = Tuple(mod(targets[i] - vals[i], primes[i]) for i in eachindex(primes))
        if haskey(table, need)
            coeffs = vcat(table[need], c_right)
            score = _coeff_score(coeffs)
            if best === nothing || score < best_score
                best = coeffs
                best_score = score
            end
        end
    end
    return best
end

function _cyclotomic_element_from_coeffs(coeffs::Vector{Int}, denom::Int,
                                         ctx::CyclotomicContext)
    K = field(ctx)
    z = zeta(ctx)
    x = zero(K)
    zpow = one(K)
    for c in coeffs
        x += K(c) * zpow
        zpow *= z
    end
    return x // K(denom)
end

function _reconstruct_cyclotomic_element_from_residues(values_by_prime::Dict{Int, Int},
                                                       ctx::CyclotomicContext;
                                                       coeff_bound::Int = 4,
                                                       denominator_bound::Int = 4)
    isempty(values_by_prime) && return nothing
    primes = sort(collect(keys(values_by_prime)))
    zetas = [find_zeta_in_Fp(ctx.N, p) for p in primes]
    degree_K = degree(field(ctx))

    best_coeffs = nothing
    best_denom = 0
    best_score = nothing
    for denom in 1:denominator_bound
        any(p -> denom % p == 0, primes) && continue
        targets = [mod(denom * values_by_prime[p], p) for p in primes]
        coeffs = _reconstruct_power_basis_coeffs_mitm(targets, zetas, primes,
                                                      degree_K, coeff_bound)
        coeffs === nothing && continue
        (l1, max_abs, _) = _coeff_score(coeffs)
        score = (l1, max_abs, denom, coeffs)
        if best_coeffs === nothing || score < best_score
            best_coeffs = coeffs
            best_denom = denom
            best_score = score
        end
    end
    best_coeffs === nothing && return nothing
    return _cyclotomic_element_from_coeffs(best_coeffs, best_denom, ctx)
end

function _estimate_modular_tuple_count(modular_filters)
    isempty(modular_filters) && return 0
    total = BigInt(1)
    for filter in modular_filters
        total *= length(filter.points)
    end
    return total
end

function _crt_reconstruction_is_well_determined(ctx::CyclotomicContext,
                                                modular_filters;
                                                reconstruction_bound::Int,
                                                denominator_bound::Int)
    modulus = BigInt(1)
    for filter in modular_filters
        modulus *= filter.p
    end
    degree_K = degree(field(ctx))
    search_space = BigInt(denominator_bound) * BigInt(2 * reconstruction_bound + 1)^degree_K
    return modulus >= search_space
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
    well_determined = _crt_reconstruction_is_well_determined(ctx, modular_filters;
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
            x = _reconstruct_cyclotomic_element_from_residues(residues, ctx;
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
