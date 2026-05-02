"""
    PentagonSolver

Exact pentagon solver over a fixed cyclotomic context.

The solver keeps the old public entry point name
`solve_pentagon_modular_crt`, but it no longer reconstructs floating
representatives.  It caches modular Groebner preprocessing over several
finite fields, uses complete split-prime point enumerations to filter
and reconstruct exact candidates over `Q(ζ_N)`, then falls back to the
exact triangular system when the bounded CRT lift is inconclusive.
"""

using Oscar

const _PENTAGON_MODULAR_GB_CACHE = Dict{Tuple{UInt, Int, Tuple{Vararg{Int}}}, Any}()

function _rational_mod(c, p::Int)
    num = mod(BigInt(numerator(c)), p)
    den = mod(BigInt(denominator(c)), p)
    den == 0 && error("coefficient denominator is zero modulo $p")
    return mod(num * invmod(den, p), p)
end

function _reduce_qq_poly_to_fp(f, vars_fp, p::Int)
    iszero(f) && return zero(vars_fp[1])
    out = zero(vars_fp[1])
    for (c, m) in zip(coefficients(f), monomials(f))
        term = parent(vars_fp[1])(_rational_mod(c, p))
        degs = degrees(m)
        for i in 1:length(degs)
            degs[i] > 0 && (term *= vars_fp[i]^degs[i])
        end
        out += term
    end
    return out
end

function _compute_modular_groebner_data(eqs, n::Int, primes::Vector{Int})
    (n == 0 || isempty(eqs)) && return Any[]
    key = (hash(string(eqs)), n, Tuple(primes))
    haskey(_PENTAGON_MODULAR_GB_CACHE, key) && return _PENTAGON_MODULAR_GB_CACHE[key]

    data = Any[]
    for p in primes
        try
            Fp = GF(p)
            Rfp, vars_fp = polynomial_ring(Fp, n, :x)
            eqs_fp = [_reduce_qq_poly_to_fp(eq, vars_fp, p) for eq in eqs]
            I = ideal(Rfp, eqs_fp)
            G = groebner_basis(I, ordering = lex(Rfp))
            push!(data, (p = p, ring = Rfp, vars = vars_fp,
                         equations = eqs_fp, ideal = I, gb = G))
        catch err
            @warn "pentagon modular Groebner basis failed at p=$p" exception = (err, catch_backtrace())
        end
    end
    _PENTAGON_MODULAR_GB_CACHE[key] = data
    return data
end

function _verify_exact_solution(eqs, sol)
    isempty(eqs) && return true
    K = parent(sol[1])
    for eq in eqs
        v = zero(K)
        for (c, m) in zip(coefficients(eq), monomials(eq))
            term = K(c)
            degs = degrees(m)
            for i in 1:length(degs)
                degs[i] > 0 && (term *= sol[i]^degs[i])
            end
            v += term
        end
        iszero(v) || error("exact pentagon solution failed verification: $v")
    end
    return true
end

"""
    solve_pentagon_modular_crt(eqs, n; Nijk, context, conductor, primes, max_solutions)

Run finite-field Groebner preprocessing and return exact F-symbol
coordinates in the selected cyclotomic field. Complete modular point
enumerations at primes split by the conductor are first lifted by bounded
CRT in the power basis of `Q(ζ_N)`. If this does not reconstruct a
verified exact solution, the solver falls back to triangular Groebner over
`Q(ζ_N)` and uses the modular points as a consistency filter.
"""
function solve_pentagon_modular_crt(eqs, n::Int;
                                    Nijk::Union{Array{Int,3}, Nothing} = nothing,
                                    context = nothing,
                                    conductor = nothing,
                                    N = nothing,
                                    primes::Vector{Int} = [101, 103, 107, 109],
                                    max_solutions::Int = 32,
                                    max_modular_solutions::Int = max(1024, 16 * max_solutions),
                                    reconstruction_bound::Int = 4,
                                    denominator_bound::Int = 4,
                                    max_crt_tuples::Int = 4096,
                                    max_ambiguous_crt_coords::Int = 4096,
                                    exact_fallback::Bool = true,
                                    show_progress::Bool = false,
                                    kwargs...)
    isempty(kwargs) || error("unsupported keyword arguments: $(collect(keys(kwargs)))")
    ctx = _default_context_from_kwargs(context = context, conductor = conductor, N = N)
    gb_data = _compute_modular_groebner_data(eqs, n, primes)
    show_progress && println("  pentagon F_p Groebner: $(length(gb_data)) primes")

    modular_filters = Any[]
    for data in gb_data
        (data.p - 1) % ctx.N == 0 || continue
        points = _enumerate_modular_triangular_solutions(data;
                                                         max_points = max_modular_solutions)
        show_progress && println("  pentagon F_$(data.p) points: $(length(points.points))" *
                                 (points.complete ? "" : " (truncated)"))
        points.complete || continue
        isempty(points.points) && return Vector{elem_type(field(ctx))}[]
        push!(modular_filters, (p = data.p, points = points.points))
    end

    crt_sols = _solve_via_modular_crt(eqs, n, ctx, modular_filters;
                                      max_solutions = max_solutions,
                                      reconstruction_bound = reconstruction_bound,
                                      denominator_bound = denominator_bound,
                                      max_crt_tuples = max_crt_tuples,
                                      max_ambiguous_crt_coords = max_ambiguous_crt_coords)
    if crt_sols !== nothing && !isempty(crt_sols)
        show_progress && println("  pentagon CRT lift: $(length(crt_sols)) exact solution(s)")
        for sol in crt_sols
            length(sol) == n || error("pentagon solution length $(length(sol)) != variable count $n")
            _verify_exact_solution(eqs, sol)
        end
        return crt_sols
    end
    exact_fallback || error("pentagon modular CRT did not reconstruct a verified exact solution")

    exact_max_solutions = isempty(modular_filters) ? max_solutions :
                          min(max_solutions, minimum(length(f.points) for f in modular_filters))
    sols = _solve_exact_via_triangular_groebner(eqs, n, field(ctx);
                                                var_prefix = :x,
                                                max_solutions = exact_max_solutions)
    for filter in modular_filters
        sols = _filter_exact_solutions_by_modular_points(sols, filter.points, ctx, filter.p)
        isempty(sols) && return sols
    end
    for sol in sols
        length(sol) == n || error("pentagon solution length $(length(sol)) != variable count $n")
        _verify_exact_solution(eqs, sol)
    end
    return sols
end

solve_pentagon_homotopy(eqs, n::Int; kwargs...) =
    solve_pentagon_modular_crt(eqs, n; kwargs...)

solve_pentagon_newton(eqs, n, zeta; kwargs...) =
    solve_pentagon_modular_crt(eqs, n; kwargs...)
