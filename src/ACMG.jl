"""
    ACMG — Arithmetic Condensed Matter Geometry

Modular Tensor Category classification by fixing the conductor `N`,
running an SL(2, ℤ/N) stratum + block-U enumeration, multi-prime F_p
sweep + CRT reconstruction, and pentagon/hexagon `(F, R)` solution.

Top-level pipeline:

    using ACMG
    auto = classify_mtcs_auto(24; skip_FR = true)
    classified = auto.classified

returns a `Vector{ClassifiedMTC}`, one per Galois sector per stratum.

Module organisation (all at ACMG top level — no submodules):
- Types:               Core layer (data types + F_p arithmetic +
                       modular-data validation + Verlinde extraction)
- SL2Reps:             SL(2, ℤ/N) irreducible representation catalog
                       (Oscar + GAP/SL2Reps)                  [Phase 0]
- StratumEnum:         combinatorial partition `Σ m_λ d_λ = r`[Phase 1]
- BlockU:              general O(n) Cayley + single-prime driver
                       `find_mtcs_at_prime`                   [Phase 2]
- CRT:                 multi-prime CRT reconstruction + Galois-aware
                       grouping                               [Phase 3]
- PentagonEquations:   pentagon equations from fusion rule
                       (TensorCategories wrapper)             [Phase 4]
- PentagonSolver:      HC + damped Newton pentagon solvers    [Phase 4]
- HexagonEquations:    hexagon equations with F baked in      [Phase 4]
- HexagonSolver:       HC solver for R                        [Phase 4]
- ModularDataLift:     F_p / ℤ[√d] → ℂ lift                   [Phase 4]
- Verify:              pentagon / hexagon / ribbon residuals,
                       `VerifyReport`, `verify_mtc`           [Phase 4]
- Pipeline:            end-to-end driver `classify_mtcs_at_conductor`,
                       auto wrapper `classify_mtcs_auto`,
                       `compute_FR_from_ST`, `classify_from_group`,
                       and the `ClassifiedMTC` output type    [Phase 5]
"""
module ACMG

using LinearAlgebra
using Primes

# Core types + arithmetic + modular-data/fusion helpers
include("Types.jl")


# SL(2, ℤ/N) irrep catalog (Oscar + GAP/SL2Reps) — Phase 0
include("SL2Reps.jl")

# Stratum enumeration (combinatorial partition of rank by irrep dimensions) — Phase 1
include("StratumEnum.jl")

# Block-U parametrisation and MTC reconstruction — Phase 2
include("BlockU.jl")

# Multi-prime CRT reconstruction — Phase 3
include("CRT.jl")

# Phase 4: Pentagon/Hexagon (F, R) solver + verify, all at top level
include("PentagonEquations.jl")
include("PentagonSolver.jl")
include("HexagonEquations.jl")
include("HexagonSolver.jl")
include("ModularDataLift.jl")
include("Verify.jl")

# Prime selection helpers
include("PrimeSelection.jl")

# End-to-end pipeline driver
include("Pipeline.jl")

# ============================================================
#  Exports
# ============================================================

# Core types and F_p arithmetic
export ModularDatumFp, FusionRule
export validate_modular_data, build_modular_datum, compute_alpha, compute_charge_conjugation
export extract_fusion_rule_Fp, lift_fusion_to_Z, extract_and_lift
export verlinde_coefficient
export is_square, sqrt_mod, primitive_root, root_of_unity, roots_of_unity
export matmul_mod, matpow_mod, diagmul_right, diagmul_left, lift_symmetric
export fusion_isomorphic, fusion_matrix, validate

# Phase 0: atomic catalog
export AtomicIrrep, build_atomic_catalog, all_divisors

# Phase 1: strata
export Stratum, enumerate_strata, count_strata, describe_stratum, find_unit_indices

# Phase 2: block-U
export MTCCandidate, build_block_diagonal, reduce_matrix_to_Fp, reduce_vector_to_Fp
export find_zeta_in_Fp, cyclotomic_to_Fp
export t_eigenspace_decomposition, parameter_dim
export o2_circle_points, apply_o2_block, verlinde_find_unit
export cayley_so_n, inverse_mod_p, enumerate_so_n_Fp, enumerate_o_n_Fp, apply_block_U
export find_mtcs_at_prime, signed_Fp

# Phase 3: CRT reconstruction
export acmg_crt, crt2, rational_reconstruct, compute_sqrt_d_mod_p
export compute_sqrt3_cyclotomic_mod_p, compute_sqrt2_cyclotomic_mod_p
export compute_sqrt5_cyclotomic_mod_p
export fusion_signature, group_mtcs_by_fusion, group_mtcs_galois_aware
export reconstruct_rational, reconstruct_in_Z_sqrt_d
export reconstruct_matrix_in_Z_sqrt_d, reconstruct_S_matrix
export verify_reconstruction, describe_matrix

# Phase 4: (F, R) classification
export get_pentagon_system
export solve_pentagon_newton, solve_pentagon_homotopy, refine_solution_newton
export eval_poly_complex, sparse_jacobian
export oscar_poly_to_hc, build_hc_system
export hexagon_equations, get_hexagon_system
export number_of_variables_in_hexagon_equations
export coerce_complex, invert_associator_numeric, assign_F_to_associator!
export solve_hexagon_homotopy
export oscar_poly_to_hc_complex, build_hc_system_complex
export DiscreteLogTable, lift_T_Fp_to_complex, lift_S_sqrtd_to_complex
export lift_mtc_candidate
export pentagon_residuals, hexagon_residuals
export extract_R_block, block_positions_R
export ribbon_residuals, VerifyReport, verify_mtc

# End-to-end pipeline
export ClassifiedMTC
export select_admissible_primes
export compute_FR_from_ST, classify_from_group
export classify_mtcs_at_conductor, classify_mtcs_auto

end # module ACMG
