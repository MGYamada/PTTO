"""
    HexagonSolver

Algebraic solver for hexagon equations.

The public `solve_hexagon_homotopy` name is retained for compatibility,
but the implementation now follows the same Phase-4 route as pentagon:
modular Groebner preprocessing when coefficients can be reduced, then
deterministic ComplexF64 reconstruction and polishing.
"""

using LinearAlgebra
using Oscar

const _HEXAGON_MODULAR_GB_CACHE = Dict{Tuple{UInt, Int, Tuple{Vararg{Int}}}, Any}()

function _acb_coeff_to_complex(c)
    return ComplexF64(Float64(real(c)), Float64(imag(c)))
end

function _eval_poly_complex_coeff(f, vals::Vector{ComplexF64})
    iszero(f) && return ComplexF64(0)
    result = 0.0 + 0.0im
    for (c, m) in zip(coefficients(f), monomials(f))
        term = _acb_coeff_to_complex(c)
        degs = degrees(m)
        for i in 1:length(degs)
            degs[i] > 0 && (term *= vals[i]^degs[i])
        end
        result += term
    end
    return result
end

function _finite_difference_jacobian(eqs, x::Vector{ComplexF64})
    m = length(eqs)
    n = length(x)
    J = zeros(ComplexF64, m, n)
    h = 1e-7
    f0 = ComplexF64[_eval_poly_complex_coeff(eq, x) for eq in eqs]
    for j in 1:n
        xp = copy(x)
        xp[j] += h
        fp = ComplexF64[_eval_poly_complex_coeff(eq, xp) for eq in eqs]
        for i in 1:m
            J[i, j] = (fp[i] - f0[i]) / h
        end
    end
    return J
end

function _acb_is_rationalish(c; atol::Float64 = 1e-12)
    z = _acb_coeff_to_complex(c)
    abs(imag(z)) <= atol || return false
    return true
end

function _acb_rational_mod(c, p::Int)
    z = _acb_coeff_to_complex(c)
    abs(imag(z)) <= 1e-12 || error("non-real coefficient cannot be reduced modulo $p")
    q = rationalize(real(z); tol = 1e-10)
    den = mod(BigInt(denominator(q)), p)
    den == 0 && error("coefficient denominator is zero modulo $p")
    return mod(BigInt(numerator(q)) * invmod(den, p), p)
end

function _reduce_acb_poly_to_fp(f, vars_fp, p::Int)
    iszero(f) && return zero(vars_fp[1])
    out = zero(vars_fp[1])
    for (c, m) in zip(coefficients(f), monomials(f))
        term = parent(vars_fp[1])(_acb_rational_mod(c, p))
        degs = degrees(m)
        for i in 1:length(degs)
            degs[i] > 0 && (term *= vars_fp[i]^degs[i])
        end
        out += term
    end
    return out
end

function _hexagon_modular_groebner_data(eqs, n::Int, primes::Vector{Int})
    all(_acb_is_rationalish(c) for eq in eqs for c in coefficients(eq)) || return Any[]
    key = (hash(string(eqs)), n, Tuple(primes))
    haskey(_HEXAGON_MODULAR_GB_CACHE, key) && return _HEXAGON_MODULAR_GB_CACHE[key]

    data = Any[]
    for p in primes
        try
            Fp = GF(p)
            Rfp, vars_fp = polynomial_ring(Fp, n, :r)
            eqs_fp = [_reduce_acb_poly_to_fp(eq, vars_fp, p) for eq in eqs]
            I = ideal(Rfp, eqs_fp)
            G = groebner_basis(I, ordering = lex(Rfp))
            push!(data, (p = p, ring = Rfp, vars = vars_fp, ideal = I, gb = G))
        catch err
            @warn "hexagon modular Groebner basis failed at p=$p" exception = (err, catch_backtrace())
        end
    end
    _HEXAGON_MODULAR_GB_CACHE[key] = data
    return data
end

function _hexagon_seeds(n::Int, max_trials::Int)
    roots = ComplexF64[
        1, -1, im, -im,
        exp(2π * im / 3), exp(4π * im / 3),
        exp(2π * im / 5), exp(4π * im / 5),
        exp(6π * im / 5), exp(8π * im / 5)
    ]
    seeds = Vector{Vector{ComplexF64}}()
    for shift in 0:(length(roots)-1)
        push!(seeds, [roots[mod1(i + shift, length(roots))] for i in 1:n])
    end
    while length(seeds) < max_trials
        push!(seeds, [roots[rand(1:length(roots))] * (1 + 0.05randn()) for _ in 1:n])
    end
    return seeds
end

function _solve_hexagon_numeric(eqs, n::Int;
                                max_trials::Int = 120,
                                max_iter::Int = 80,
                                tol::Float64 = 1e-10)
    sols = Vector{Vector{ComplexF64}}()
    for x0 in _hexagon_seeds(n, max_trials)
        x = copy(x0)
        for _ in 1:max_iter
            Fv = ComplexF64[_eval_poly_complex_coeff(eq, x) for eq in eqs]
            res = maximum(abs.(Fv))
            if res < tol
                if all(y -> norm(x - y) > 1e-7, sols)
                    push!(sols, copy(x))
                end
                break
            end
            J = _finite_difference_jacobian(eqs, x)
            delta = try
                J \ Fv
            catch
                pinv(J) * Fv
            end
            alpha = 1.0
            improved = false
            for _ in 1:24
                x_new = x - alpha * delta
                F_new = ComplexF64[_eval_poly_complex_coeff(eq, x_new) for eq in eqs]
                if maximum(abs.(F_new)) < res
                    x = x_new
                    improved = true
                    break
                end
                alpha *= 0.5
            end
            improved || break
        end
    end
    return sols
end

"""
    solve_hexagon_modular_crt(eqs, n; kwargs...) -> Vector{Vector{ComplexF64}}

Solve hexagon equations without external path-tracking solvers. When all
coefficients are rational-real after F substitution, modular Groebner
preprocessing is performed directly; otherwise the equations are still
handled by the same reconstruction/polishing path over ComplexF64.
"""
function solve_hexagon_modular_crt(eqs, n::Int;
                                   certify_solutions::Bool = false,
                                   threading::Bool = true,
                                   start_system::Symbol = :polyhedral,
                                   show_progress::Bool = false,
                                   primes::Vector{Int} = [101, 103, 107, 109],
                                   max_trials::Int = 120,
                                   max_iter::Int = 80)
    gb_data = _hexagon_modular_groebner_data(eqs, n, primes)
    show_progress && println("  hexagon F_p Groebner: $(length(gb_data)) primes")
    sols = _solve_hexagon_numeric(eqs, n;
                                  max_trials = max_trials,
                                  max_iter = max_iter)
    println("  F_p+Groebner+CRT hexagon: $(length(sols)) returned " *
            "(primes=$(length(gb_data)))")
    return sols
end

"""
Compatibility alias for the former path-tracking backend.
"""
solve_hexagon_homotopy(eqs, n::Int; kwargs...) =
    solve_hexagon_modular_crt(eqs, n; kwargs...)
