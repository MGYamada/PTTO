"""
    HexagonSolver

Exact Phase-4 braiding solver over `Q(ζ_N)`.
"""

using Oscar

const _HEXAGON_MODULAR_GB_CACHE = Dict{Tuple{UInt, Int, Tuple{Vararg{Int}}}, Any}()

function _reduce_nf_elem_to_fp(c, zeta_Fp::Int, p::Int)
    return cyclotomic_to_Fp(c, zeta_Fp, p)
end

function _reduce_nf_poly_to_fp(f, vars_fp, zeta_Fp::Int, p::Int)
    iszero(f) && return zero(vars_fp[1])
    out = zero(vars_fp[1])
    for (c, m) in zip(coefficients(f), monomials(f))
        term = parent(vars_fp[1])(_reduce_nf_elem_to_fp(c, zeta_Fp, p))
        degs = degrees(m)
        for i in 1:length(degs)
            degs[i] > 0 && (term *= vars_fp[i]^degs[i])
        end
        out += term
    end
    return out
end

function _hexagon_modular_groebner_data(eqs, n::Int, ctx::CyclotomicContext, primes::Vector{Int})
    key = (hash(string(eqs)), n, Tuple(primes))
    haskey(_HEXAGON_MODULAR_GB_CACHE, key) && return _HEXAGON_MODULAR_GB_CACHE[key]
    data = Any[]
    for p in primes
        try
            zeta_Fp = find_zeta_in_Fp(ctx.N, p)
            Fp = GF(p)
            Rfp, vars_fp = polynomial_ring(Fp, n, :r)
            eqs_fp = [_reduce_nf_poly_to_fp(eq, vars_fp, zeta_Fp, p) for eq in eqs]
            I = ideal(Rfp, eqs_fp)
            G = groebner_basis(I, ordering = lex(Rfp))
            push!(data, (p = p, ring = Rfp, vars = vars_fp,
                         equations = eqs_fp, ideal = I, gb = G))
        catch err
            @warn "hexagon modular Groebner basis failed at p=$p" exception = (err, catch_backtrace())
        end
    end
    _HEXAGON_MODULAR_GB_CACHE[key] = data
    return data
end

function _verify_hexagon_solution(eqs, sol)
    K = parent(sol[1])
    for eq in eqs
        v = zero(K)
        for (c, m) in zip(coefficients(eq), monomials(eq))
            term = c
            degs = degrees(m)
            for i in 1:length(degs)
                degs[i] > 0 && (term *= sol[i]^degs[i])
            end
            v += term
        end
        iszero(v) || error("exact hexagon solution failed verification: $v")
    end
    return true
end

function solve_hexagon_modular_crt(eqs, n::Int;
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
    gb_data = _hexagon_modular_groebner_data(eqs, n, ctx, primes)
    show_progress && println("  hexagon F_p Groebner: $(length(gb_data)) primes")
    modular_filters = Any[]
    for data in gb_data
        points = _enumerate_modular_triangular_solutions(data;
                                                         max_points = max_modular_solutions)
        show_progress && println("  hexagon F_$(data.p) points: $(length(points.points))" *
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
        show_progress && println("  hexagon CRT lift: $(length(crt_sols)) exact solution(s)")
        for sol in crt_sols
            length(sol) == n || error("hexagon solution length $(length(sol)) != variable count $n")
            _verify_hexagon_solution(eqs, sol)
        end
        return crt_sols
    end
    exact_fallback || error("hexagon modular CRT did not reconstruct a verified exact solution")

    exact_max_solutions = isempty(modular_filters) ? max_solutions :
                          min(max_solutions, minimum(length(f.points) for f in modular_filters))
    sols = _solve_exact_via_triangular_groebner(eqs, n, field(ctx);
                                                var_prefix = :r,
                                                max_solutions = exact_max_solutions)
    for filter in modular_filters
        sols = _filter_exact_solutions_by_modular_points(sols, filter.points, ctx, filter.p)
        isempty(sols) && return sols
    end
    for sol in sols
        length(sol) == n || error("hexagon solution length $(length(sol)) != variable count $n")
        _verify_hexagon_solution(eqs, sol)
    end
    return sols
end

solve_hexagon_homotopy(eqs, n::Int; kwargs...) =
    solve_hexagon_modular_crt(eqs, n; kwargs...)
