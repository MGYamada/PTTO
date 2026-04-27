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

function _solve_exact_via_triangular_groebner(eqs, n::Int, K;
                                              var_prefix::Symbol = :x,
                                              max_solutions::Int = 32)
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
