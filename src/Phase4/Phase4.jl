"""
    Phase4

Phase 4 of ACMG: Pentagon and hexagon solvers for MTC classification.

This module wraps the four migrated submodules:
- `PentagonEquations`:   generates pentagon polynomials via TensorCategories.
- `PentagonSolver`:      HomotopyContinuation + damped Newton solvers.
- `HexagonEquations`:    generates hexagon polynomials with F fixed.
- `HexagonSolver`:       HomotopyContinuation solver for R.

Usage:

    using ACMG.Phase4

    # Pentagon: F-symbol classification for a given fusion rule
    R, eqs, n = Phase4.get_pentagon_system(Nijk, r)
    F_sols = Phase4.solve_pentagon_homotopy(eqs, n; slice = 1)

    # Hexagon: R-symbol classification for a given F
    R_ring, hex_eqs, n_r = Phase4.get_hexagon_system(Nijk, r, F_sols[1])
    R_sols = Phase4.solve_hexagon_homotopy(hex_eqs, n_r)

Each F_sol in `F_sols` is a `Vector{ComplexF64}` of pentagon-variable
values (ordering determined by TensorCategories' internal traversal).
Each R_sol in `R_sols` is a `Vector{ComplexF64}` of R-symbol values,
split as [forward braidings; reverse braidings].

Design notes:
- Pentagon/Hexagon classification is over ℂ (ComplexF64). Algebraic
  lift to ℚ(ζ_N) is left to downstream code (PSLQ or LLL-based
  recognition).
- Phase 4 is currently **standalone**: not wired to ACMG's Phase 3
  output (which produces F_p / ℤ[√d] modular data). A future
  integration layer will bridge Phase 3 → Phase 4.
"""
module Phase4

include("PentagonEquations.jl")
include("PentagonSolver.jl")
include("HexagonEquations.jl")
include("HexagonSolver.jl")

using .PentagonEquations
using .PentagonSolver
using .HexagonEquations
using .HexagonSolver

# Re-export the user-facing API
export get_pentagon_system
export solve_pentagon_newton, solve_pentagon_homotopy, refine_solution_newton
export hexagon_equations, get_hexagon_system
export solve_hexagon_homotopy
export assign_F_to_associator!, invert_associator_numeric

end # module Phase4
