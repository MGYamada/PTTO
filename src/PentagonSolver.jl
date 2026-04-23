"""
    PentagonSolver

Numerical solvers for pentagon equations.

Provides:
- `solve_pentagon_homotopy` (MAIN): polynomial homotopy continuation via
  HomotopyContinuation.jl. Near-exhaustive coverage of solutions;
  supports random linear slicing to break pentagon's gauge symmetry.
- `solve_pentagon_newton`: legacy damped Newton + KrylovKit. Kept for
  comparison and for small problems where HC overhead is unjustified.
- `refine_solution_newton`: polish an HC solution with a few Newton
  steps at ComplexF64 precision; useful before algebraic recognition.

All routines take Oscar polynomials in variables `x_1, …, x_n` with
rational coefficients (as produced by `get_pentagon_system`).

Depends on: Oscar, SparseArrays, KrylovKit, HomotopyContinuation.
"""

using Oscar
using SparseArrays
using KrylovKit
import HomotopyContinuation
const HC = HomotopyContinuation


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

Assemble the Jacobian at `x` using precomputed symbolic derivatives
`derivs[i][j] = ∂eqs[i]/∂x_j`. Entries below 1e-15 in absolute value
are dropped for sparsity.
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

# ============================================================
#  Damped Newton solver (legacy)
# ============================================================

"""
    solve_pentagon_newton(eqs, n, zeta; max_trials=20, max_iter=200)
        -> Vector{Vector{ComplexF64}}

Legacy damped-Newton solver. Draws random starting points from
`ℤ[ζ] + ζ·ℤ[ζ]` (small coefficients) and applies Newton with KrylovKit
for the normal-equations linear solve and a backtracking line search.

This solver is kept because ℂ has the right topology (line search,
damping) that 𝔽_p lacks, but it gets trapped in local minima on all but
the smallest problems. Use `solve_pentagon_homotopy` as default.
"""
function solve_pentagon_newton(eqs, n, zeta; max_trials::Int = 20, max_iter::Int = 200)
    derivs = [[derivative(eq, j) for j in 1:n] for eq in eqs]

    solutions = Vector{Vector{ComplexF64}}()
    for trial in 1:max_trials
        x = ComplexF64[rand(-3:3) + rand(-3:3) * zeta for _ in 1:n]
        for i in 1:n
            abs(x[i]) < 0.1 && (x[i] = 1.0 + zeta)
        end

        for iter in 1:max_iter
            F_val = ComplexF64[eval_poly_complex(eq, x) for eq in eqs]
            res = maximum(abs.(F_val))
            if res < 1e-12
                push!(solutions, copy(x))
                break
            end

            J = sparse_jacobian(eqs, derivs, x, n)
            delta, _ = linsolve(v -> J' * (J * v), J' * F_val;
                                ishermitian = true, isposdef = true, verbosity = 0)

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
    end
    solutions
end

# ============================================================
#  HomotopyContinuation backend
# ============================================================

"""
    oscar_poly_to_hc(f, hc_vars) -> HC.Expression

Convert an Oscar multivariate polynomial `f` (QQ coefficients) into a
HomotopyContinuation.jl Expression in the provided `hc_vars`.
"""
function oscar_poly_to_hc(f, hc_vars::Vector{HC.Variable})
    isa(f, Integer) && return HC.Expression(f)
    iszero(f) && return HC.Expression(0)
    expr = HC.Expression(0)
    for (c, m) in zip(coefficients(f), monomials(f))
        num = BigInt(numerator(c))
        den = BigInt(denominator(c))
        coef_expr = HC.Expression(num // den)
        degs = degrees(m)
        term = coef_expr
        for i in 1:length(degs)
            if degs[i] > 0
                term *= hc_vars[i]^degs[i]
            end
        end
        expr += term
    end
    expr
end

"""
    build_hc_system(eqs, n) -> (HC.System, Vector{HC.Variable})

Build a HomotopyContinuation `System` in `n` variables (named `x_1`…`x_n`)
from Oscar polynomial equations `eqs`. Filters out trivially-zero
equations.
"""
function build_hc_system(eqs, n::Int)
    hc_vars = [HC.Variable(Symbol("x", i)) for i in 1:n]
    hc_exprs = HC.Expression[]
    for eq in eqs
        (isa(eq, Integer) || iszero(eq)) && continue
        push!(hc_exprs, oscar_poly_to_hc(eq, hc_vars))
    end
    return HC.System(hc_exprs; variables = hc_vars), hc_vars
end

"""
    solve_pentagon_homotopy(eqs, n; kwargs...) -> Vector{Vector{ComplexF64}}

Solve pentagon equations via polynomial homotopy continuation.

Pentagon carries a diagonal-basis **gauge symmetry**, so the raw variety
is positive-dimensional and HC tends to flag all endpoints as singular.
The `slice` option appends random complex-linear forms `aᵀx = b` to
cut gauge directions and recover a zero-dimensional system.

Keyword arguments:
- `slice = 0`:                   number of random linear slices to append
- `include_singular = false`:    if true, return singular endpoints too
- `certify_solutions = false`:   rigorous certification via HC.certify
- `threading = true`:            parallel path tracking
- `start_system = :polyhedral`:  :polyhedral or :total_degree
- `show_progress = false`:       HC progress bar
"""
function solve_pentagon_homotopy(eqs, n::Int;
                                 slice::Int = 0,
                                 include_singular::Bool = false,
                                 certify_solutions::Bool = false,
                                 threading::Bool = true,
                                 start_system::Symbol = :polyhedral,
                                 show_progress::Bool = false)
    sys, hc_vars = build_hc_system(eqs, n)

    if slice > 0
        extra = HC.Expression[]
        for _ in 1:slice
            a = randn(ComplexF64, n)
            b = randn(ComplexF64)
            lf = sum(HC.Expression(a[i]) * hc_vars[i] for i in 1:n) - HC.Expression(b)
            push!(extra, lf)
        end
        all_eqs = vcat(HC.expressions(sys), extra)
        sys = HC.System(all_eqs; variables = hc_vars)
    end

    result = HC.solve(sys;
                      start_system = start_system,
                      threading = threading,
                      show_progress = show_progress)

    sols_raw = if include_singular
        HC.solutions(result)
    else
        HC.solutions(result; only_nonsingular = true)
    end

    if certify_solutions && !isempty(sols_raw)
        cert = HC.certify(sys, sols_raw)
        certs = HC.certificates(cert)
        sols_raw = [HC.solution_approximation(c) for c in certs if HC.is_certified(c)]
        println("  Certified: $(length(sols_raw)) / $(length(certs)) solutions")
    end

    nsols_total = HC.nsolutions(result)
    nsing = HC.nsingular(result)
    natinf = HC.nat_infinity(result)
    nfail = HC.nfailed(result)
    slicestr = slice > 0 ? " [+$slice slice]" : ""
    println("  HC$slicestr: $(length(sols_raw)) returned " *
            "(total=$nsols_total, singular=$nsing, at_infinity=$natinf, failed=$nfail)")

    return [ComplexF64.(s) for s in sols_raw]
end

# ============================================================
#  Solution polishing
# ============================================================

"""
    refine_solution_newton(eqs, x0; tol=1e-14, max_iter=50) -> Vector{ComplexF64}

Polish a solution `x0` with damped Newton steps at ComplexF64 precision.
Useful to tighten HC's double-precision output before PSLQ / algebraic
recognition downstream.
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
        J = sparse_jacobian(eqs, derivs, x, n)
        delta, _ = linsolve(v -> J' * (J * v), J' * F_val;
                            ishermitian = true, isposdef = true, verbosity = 0)
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

