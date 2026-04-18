"""
    Phase4

Phase 4 of ACMG: Pentagon and hexagon solvers for MTC classification.

This module wraps five migrated / new submodules:
- `PentagonEquations`:   generates pentagon polynomials via TensorCategories.
- `PentagonSolver`:      HomotopyContinuation + damped Newton solvers.
- `HexagonEquations`:    generates hexagon polynomials with F fixed.
- `HexagonSolver`:       HomotopyContinuation solver for R.
- `ModularDataLift`:     bridge from v0.2 Phase 3 output (F_p + ℤ[√d])
                         to Phase 4 input (ℂ).

Usage:

    using ACMG.Phase4

    # Pentagon: F-symbol classification for a given fusion rule
    R, eqs, n = Phase4.get_pentagon_system(Nijk, r)
    F_sols = Phase4.solve_pentagon_homotopy(eqs, n; slice = 1)

    # Hexagon: R-symbol classification for a given F
    R_ring, hex_eqs, n_r = Phase4.get_hexagon_system(Nijk, r, F_sols[1])
    R_sols = Phase4.solve_hexagon_homotopy(hex_eqs, n_r)

    # Bridge from v0.2 Phase 3:
    zeta_Fp = ACMG.find_zeta_in_Fp(N, candidate.p)
    S_ℂ, T_ℂ, Nijk = Phase4.lift_mtc_candidate(candidate, recon_S;
                                                d = 3, N = 24,
                                                zeta_Fp = zeta_Fp)

Each F_sol in `F_sols` is a `Vector{ComplexF64}` of pentagon-variable
values (ordering determined by TensorCategories' internal traversal).
Each R_sol in `R_sols` is a `Vector{ComplexF64}` of R-symbol values,
split as [forward braidings; reverse braidings].

Design notes:
- Pentagon/Hexagon classification is over ℂ (ComplexF64). Algebraic
  lift to ℚ(ζ_N) is left to downstream code (PSLQ or LLL-based
  recognition).
- ModularDataLift accepts Phase 3 output from v0.2 (MTCCandidate and
  `reconstruct_S_matrix` result) and returns (S, T, Nijk) in ℂ.
  Verification of a Phase 4 solution against the original (S, T) is
  delegated to Verify.jl (forthcoming).
"""
module Phase4

include("PentagonEquations.jl")
include("PentagonSolver.jl")
include("HexagonEquations.jl")
include("HexagonSolver.jl")
include("ModularDataLift.jl")
include("Verify.jl")
include("KitaevComplex.jl")

using .PentagonEquations
using .PentagonSolver
using .HexagonEquations
using .HexagonSolver
using .ModularDataLift
using .Verify
using .KitaevComplex

# Re-export the user-facing API
export get_pentagon_system
export solve_pentagon_newton, solve_pentagon_homotopy, refine_solution_newton
export hexagon_equations, get_hexagon_system
export solve_hexagon_homotopy
export assign_F_to_associator!, invert_associator_numeric
# ModularDataLift
export DiscreteLogTable, lift_T_Fp_to_complex, lift_S_sqrtd_to_complex
export lift_mtc_candidate
# Verify
export pentagon_residuals, hexagon_residuals
export extract_R_block, block_positions_R
export ribbon_residuals, VerifyReport, verify_mtc

# KitaevComplex (Kitaev 2006 App E.6 tangent cohomology + slice)
export KitaevCx, build_kitaev_complex, build_complex
export Cn_dim, Cn_basis
export delta_matrix, chi_matrix
export verify_homotopy, verify_complex
export slice_constraint, kernel_chi
export tangent_cohomology_dims

end # module Phase4
