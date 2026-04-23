"""
    HexagonSolver

HomotopyContinuation solver for hexagon equations.

Unlike pentagon, hexagon has NO continuous gauge freedom (pentagon has
already exhausted the F-gauge; hexagon adds discrete R-data on top of a
fixed F). The system should therefore be zero-dimensional with all
endpoints nonsingular.

The hexagon polynomials live over `AcbField` (Arb complex field) because
the F-values are ComplexF64 constants absorbed into the coefficients.
HC requires ComplexF64 coefficients, so conversion goes through
`oscar_poly_to_hc_complex`.

Depends on: Oscar, HomotopyContinuation.
"""

using Oscar
import HomotopyContinuation
# `const HC = HomotopyContinuation` is defined in PentagonSolver.jl; it is
# shared across ACMG after flattening and we do not re-declare it here.


"""
    oscar_poly_to_hc_complex(f, hc_vars) -> HC.Expression

Convert an Oscar polynomial with **AcbField (complex)** coefficients into
a HomotopyContinuation Expression. Coefficients are converted via
ComplexF64 (midpoints of Arb intervals).
"""
function oscar_poly_to_hc_complex(f, hc_vars::Vector{HC.Variable})
    iszero(f) && return HC.Expression(0)
    expr = HC.Expression(0)
    for (c, m) in zip(coefficients(f), monomials(f))
        # c is an AcbFieldElem; convert to ComplexF64 via real/imag midpoints
        cre = Float64(real(c))
        cim = Float64(imag(c))
        coef = ComplexF64(cre, cim)
        term = HC.Expression(coef)
        degs = degrees(m)
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
    build_hc_system_complex(eqs, n) -> (HC.System, Vector{HC.Variable})

Build an HC `System` in `n` variables (named `r_1`…`r_n`) from Oscar
polynomial equations with AcbField coefficients. Filters trivially-zero
equations.
"""
function build_hc_system_complex(eqs, n::Int)
    hc_vars = [HC.Variable(Symbol("r", i)) for i in 1:n]
    hc_exprs = HC.Expression[]
    for eq in eqs
        iszero(eq) && continue
        push!(hc_exprs, oscar_poly_to_hc_complex(eq, hc_vars))
    end
    return HC.System(hc_exprs; variables = hc_vars), hc_vars
end

"""
    solve_hexagon_homotopy(eqs, n; kwargs...) -> Vector{Vector{ComplexF64}}

Solve hexagon equations via polynomial homotopy continuation.

Hexagon has NO gauge freedom, so this routine does not offer a `slice`
option (unlike `solve_pentagon_homotopy`).

Keyword arguments:
- `certify_solutions = false`: run `HC.certify` for rigorous validation
- `threading = true`:          parallel path tracking
- `start_system = :polyhedral`: :polyhedral (sparse) or :total_degree
- `show_progress = false`
"""
function solve_hexagon_homotopy(eqs, n::Int;
                                certify_solutions::Bool = false,
                                threading::Bool = true,
                                start_system::Symbol = :polyhedral,
                                show_progress::Bool = false)
    sys, hc_vars = build_hc_system_complex(eqs, n)

    result = HC.solve(sys;
                      start_system = start_system,
                      threading = threading,
                      show_progress = show_progress)
    sols_raw = HC.solutions(result; only_nonsingular = true)

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
    println("  HC hexagon: $(length(sols_raw)) nonsingular " *
            "(total=$nsols_total, singular=$nsing, at_infinity=$natinf, failed=$nfail)")
    return [ComplexF64.(s) for s in sols_raw]
end

