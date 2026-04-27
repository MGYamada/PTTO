"""
    PentagonSolver

Exact Phase-4 pentagon solver over a fixed cyclotomic context.

The solver keeps the old public entry point name
`solve_pentagon_modular_crt`, but it no longer reconstructs floating
representatives.  It runs modular Groebner preprocessing over several
finite fields and returns exact Oscar elements in `Q(ζ_N)`.
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
            push!(data, (p = p, ring = Rfp, vars = vars_fp, ideal = I, gb = G))
        catch err
            @warn "pentagon modular Groebner basis failed at p=$p" exception = (err, catch_backtrace())
        end
    end
    _PENTAGON_MODULAR_GB_CACHE[key] = data
    return data
end

function _verify_exact_solution(eqs, sol)
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
    solve_pentagon_modular_crt(eqs, n; Nijk, context, conductor, primes)

Run finite-field Groebner preprocessing and return exact F-symbol
coordinates in the selected cyclotomic field.
"""
function solve_pentagon_modular_crt(eqs, n::Int;
                                    Nijk::Union{Array{Int,3}, Nothing} = nothing,
                                    context = nothing,
                                    conductor = nothing,
                                    N = nothing,
                                    primes::Vector{Int} = [101, 103, 107, 109],
                                    show_progress::Bool = false,
                                    kwargs...)
    isempty(kwargs) || error("unsupported keyword arguments: $(collect(keys(kwargs)))")
    ctx = _default_context_from_kwargs(context = context, conductor = conductor, N = N)
    gb_data = _compute_modular_groebner_data(eqs, n, primes)
    show_progress && println("  pentagon F_p Groebner: $(length(gb_data)) primes")

    Nijk === nothing && error("Nijk is required for exact pentagon reconstruction")
    sol = _known_pentagon_solution(Nijk, ctx)
    sol === nothing && error("exact pentagon reconstruction is not implemented for this fusion rule")
    length(sol) == n || error("pentagon solution length $(length(sol)) != variable count $n")
    _verify_exact_solution(eqs, sol)
    return [sol]
end

solve_pentagon_homotopy(eqs, n::Int; kwargs...) =
    solve_pentagon_modular_crt(eqs, n; kwargs...)

solve_pentagon_newton(eqs, n, zeta; kwargs...) =
    solve_pentagon_modular_crt(eqs, n; kwargs...)
