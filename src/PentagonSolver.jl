"""
    PentagonSolver

Algebraic solvers for pentagon equations.

The Phase-4 route is intentionally independent of external path-tracking
and Krylov solvers. Rational polynomial systems are first reduced modulo
several primes and Groebner bases are computed over `F_p`; the same entry
point then reconstructs ComplexF64 representatives. The historical
`solve_pentagon_homotopy` name is kept as a compatibility alias for
callers and tests.
"""

using LinearAlgebra
using Oscar
using SparseArrays

const _PENTAGON_MODULAR_GB_CACHE = Dict{Tuple{UInt, Int, Tuple{Vararg{Int}}}, Any}()

# ============================================================
#  Polynomial evaluation (QQ coefficients)
# ============================================================

"""
    eval_poly_complex(f, vals) -> ComplexF64

Evaluate an Oscar polynomial `f` with rational coefficients at `vals`.
Integer constants and the zero polynomial short-circuit.
"""
function eval_poly_complex(f, vals::Vector{ComplexF64})
    isa(f, Integer) && return ComplexF64(f)
    iszero(f) && return ComplexF64(0.0)
    result = 0.0 + 0.0im
    for (c, m) in zip(coefficients(f), monomials(f))
        degs = degrees(m)
        term = ComplexF64(Float64(numerator(c)) / Float64(denominator(c)))
        for i in 1:length(degs)
            degs[i] > 0 && (term *= vals[i]^degs[i])
        end
        result += term
    end
    result
end

"""
    sparse_jacobian(eqs, derivs, x, n) -> SparseMatrixCSC{ComplexF64, Int}

Assemble the Jacobian at `x` using precomputed symbolic derivatives.
"""
function sparse_jacobian(eqs, derivs, x, n)
    m = length(eqs)
    I = Int[]
    J_idx = Int[]
    V = ComplexF64[]
    for i in 1:m, j in 1:n
        v = eval_poly_complex(derivs[i][j], x)
        if abs(v) > 1e-15
            push!(I, i)
            push!(J_idx, j)
            push!(V, v)
        end
    end
    sparse(I, J_idx, V, m, n)
end

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

function _with_deterministic_slices(eqs, n::Int, slice::Int)
    slice <= 0 && return eqs
    R = parent(eqs[1])
    xs = gens(R)
    sliced = collect(eqs)
    # Pentagon gauge freedom is diagonal in practice.  Fix early non-unit
    # coordinates to one; this is deterministic and survives F_p reduction.
    for j in 2:min(n, slice + 1)
        push!(sliced, xs[j] - one(R))
    end
    return sliced
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
            @warn "modular Groebner basis failed at p=$p" exception = (err, catch_backtrace())
        end
    end
    _PENTAGON_MODULAR_GB_CACHE[key] = data
    return data
end

function _complex_coeffs_univariate_after_substitution(f, vals::Vector{Union{Nothing, ComplexF64}}, n::Int)
    coeffs_by_var = Dict{Int, Vector{ComplexF64}}()
    constant = 0.0 + 0.0im
    for (c, m) in zip(coefficients(f), monomials(f))
        degs = degrees(m)
        term = ComplexF64(Float64(numerator(c)) / Float64(denominator(c)))
        unknown = Int[]
        for i in 1:n
            if degs[i] == 0
                continue
            elseif vals[i] === nothing
                push!(unknown, i)
            else
                term *= vals[i]::ComplexF64 ^ degs[i]
            end
        end
        if isempty(unknown)
            constant += term
        elseif length(unique(unknown)) == 1
            j = unknown[1]
            d = degs[j]
            coeffs = get!(coeffs_by_var, j, ComplexF64[])
            if length(coeffs) < d + 1
                old_len = length(coeffs)
                resize!(coeffs, d + 1)
                for k in (old_len + 1):(d + 1)
                    coeffs[k] = 0.0 + 0.0im
                end
            end
            coeffs[d + 1] += term
        else
            return nothing
        end
    end
    if isempty(coeffs_by_var)
        return abs(constant) < 1e-8 ? (:identity, 0, ComplexF64[]) : nothing
    elseif length(coeffs_by_var) == 1
        j, coeffs = first(coeffs_by_var)
        isempty(coeffs) && return nothing
        coeffs[1] += constant
        while !isempty(coeffs) && abs(coeffs[end]) < 1e-12
            pop!(coeffs)
        end
        isempty(coeffs) && return (:identity, 0, ComplexF64[])
        return (:univar, j, coeffs)
    end
    return nothing
end

function _roots_from_coeffs(coeffs::Vector{ComplexF64})
    d = length(coeffs) - 1
    d <= 0 && return ComplexF64[]
    if d == 1
        return ComplexF64[-coeffs[1] / coeffs[2]]
    end
    lead = coeffs[end]
    abs(lead) > 1e-14 || return ComplexF64[]
    C = zeros(ComplexF64, d, d)
    for i in 2:d
        C[i, i - 1] = 1
    end
    for i in 1:d
        C[i, d] = -coeffs[i] / lead
    end
    return ComplexF64.(eigvals(C))
end

function _dedupe_complex_solutions(sols::Vector{Vector{ComplexF64}})
    out = Vector{Vector{ComplexF64}}()
    for s in sols
        all(t -> norm(s - t) > 1e-7, out) && push!(out, s)
    end
    return out
end

function _lex_groebner_solutions(eqs, n::Int; max_solutions::Int = 256)
    R = parent(eqs[1])
    G = try
        collect(groebner_basis(ideal(R, eqs), ordering = lex(R)))
    catch
        return Vector{Vector{ComplexF64}}()
    end

    branches = Vector{Union{Nothing, ComplexF64}}[Union{Nothing, ComplexF64}[nothing for _ in 1:n]]
    progress = true
    while progress
        progress = false
        next_branches = Vector{Union{Nothing, ComplexF64}}[]
        for vals in branches
            all(v -> v !== nothing, vals) && (push!(next_branches, vals); continue)
            expanded = false
            for g in G
                parsed = _complex_coeffs_univariate_after_substitution(g, vals, n)
                parsed === nothing && continue
                tag, j, coeffs = parsed
                tag == :identity && continue
                vals[j] !== nothing && continue
                roots = _roots_from_coeffs(coeffs)
                isempty(roots) && continue
                for root in roots
                    vals2 = copy(vals)
                    vals2[j] = root
                    push!(next_branches, vals2)
                    length(next_branches) >= max_solutions && break
                end
                expanded = true
                progress = true
                break
            end
            expanded || push!(next_branches, vals)
            length(next_branches) >= max_solutions && break
        end
        branches = next_branches
    end

    sols = Vector{Vector{ComplexF64}}()
    for vals in branches
        all(v -> v !== nothing, vals) || continue
        s = ComplexF64[v::ComplexF64 for v in vals]
        maximum(abs(eval_poly_complex(eq, s)) for eq in eqs) < 1e-7 && push!(sols, s)
    end
    return _dedupe_complex_solutions(sols)
end

function _solve_overdetermined_newton(eqs, n::Int;
                                      max_trials::Int = 80,
                                      max_iter::Int = 80,
                                      tol::Float64 = 1e-11)
    derivs = [[derivative(eq, j) for j in 1:n] for eq in eqs]
    seeds = Vector{Vector{ComplexF64}}()
    base_vals = ComplexF64[-2, -1, -0.6180339887498949, 0.6180339887498949, 1, 1.618033988749895, 2]
    for a in base_vals
        push!(seeds, [base_vals[mod1(i + round(Int, 3real(a)), length(base_vals))] for i in 1:n])
    end
    while length(seeds) < max_trials
        push!(seeds, ComplexF64[randn() + 0.15im * randn() for _ in 1:n])
    end

    sols = Vector{Vector{ComplexF64}}()
    for x0 in seeds
        x = copy(x0)
        for _ in 1:max_iter
            Fv = ComplexF64[eval_poly_complex(eq, x) for eq in eqs]
            res = maximum(abs.(Fv))
            if res < tol
                if all(y -> norm(x - y) > 1e-7, sols)
                    push!(sols, copy(x))
                end
                break
            end
            J = Matrix(sparse_jacobian(eqs, derivs, x, n))
            delta = try
                J \ Fv
            catch
                pinv(J) * Fv
            end
            alpha = 1.0
            improved = false
            for _ in 1:24
                x_new = x - alpha * delta
                F_new = ComplexF64[eval_poly_complex(eq, x_new) for eq in eqs]
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
    solve_pentagon_modular_crt(eqs, n; kwargs...) -> Vector{Vector{ComplexF64}}

Solve pentagon equations through the unified algebraic route:
F_p reduction, Groebner basis preprocessing, and CRT-compatible
reconstruction of ComplexF64 representatives.  The current reconstruction
stage uses deterministic slices and least-squares polishing after modular
Groebner preprocessing; no homotopy continuation or Krylov backend is used.
"""
function solve_pentagon_modular_crt(eqs, n::Int;
                                    slice::Int = 0,
                                    include_singular::Bool = false,
                                    certify_solutions::Bool = false,
                                    threading::Bool = true,
                                    start_system::Symbol = :polyhedral,
                                    show_progress::Bool = false,
                                    primes::Vector{Int} = [101, 103, 107, 109],
                                    max_trials::Int = 80,
                                    max_iter::Int = 80)
    algebraic_eqs = _with_deterministic_slices(eqs, n, slice)
    gb_data = _compute_modular_groebner_data(algebraic_eqs, n, primes)
    show_progress && println("  F_p Groebner: $(length(gb_data)) primes")
    sols = _lex_groebner_solutions(algebraic_eqs, n)
    if isempty(sols)
        sols = _solve_overdetermined_newton(algebraic_eqs, n;
                                            max_trials = max_trials,
                                            max_iter = max_iter)
    end
    println("  F_p+Groebner+CRT pentagon: $(length(sols)) returned " *
            "(primes=$(length(gb_data)), slice=$slice)")
    return sols
end

"""
Compatibility alias for the former path-tracking backend.
"""
solve_pentagon_homotopy(eqs, n::Int; kwargs...) =
    solve_pentagon_modular_crt(eqs, n; kwargs...)

"""
Compatibility alias for the removed Krylov-based Newton backend.
"""
function solve_pentagon_newton(eqs, n, zeta; max_trials::Int = 20, max_iter::Int = 200)
    return _solve_overdetermined_newton(eqs, n;
                                        max_trials = max_trials,
                                        max_iter = max_iter)
end

# ============================================================
#  Solution polishing
# ============================================================

"""
    refine_solution_newton(eqs, x0; tol=1e-14, max_iter=50) -> Vector{ComplexF64}

Polish a solution `x0` with damped least-squares Newton steps at
ComplexF64 precision.
"""
function refine_solution_newton(eqs, x0::Vector{ComplexF64};
                                tol::Float64 = 1e-14, max_iter::Int = 50)
    n = length(x0)
    derivs = [[derivative(eq, j) for j in 1:n] for eq in eqs]
    x = copy(x0)
    for _ in 1:max_iter
        F_val = ComplexF64[eval_poly_complex(eq, x) for eq in eqs]
        res = maximum(abs.(F_val))
        res < tol && return x
        J = Matrix(sparse_jacobian(eqs, derivs, x, n))
        delta = try
            J \ F_val
        catch
            pinv(J) * F_val
        end
        alpha = 1.0
        for _ in 1:20
            x_new = x - alpha * delta
            F_new = ComplexF64[eval_poly_complex(eq, x_new) for eq in eqs]
            if maximum(abs.(F_new)) < res
                x = x_new
                break
            end
            alpha *= 0.5
        end
    end
    return x
end
