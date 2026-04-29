"""
    ACMG — Arithmetic Condensed Matter Geometry

Modular Tensor Category classification by fixing the conductor `N`,
running an SL(2, ℤ/N) stratum + block-U enumeration, multi-prime F_p
sweep + CRT reconstruction, and exact lift to `Q(ζ_N)`.

Top-level pipeline:

    using ACMG
    auto = classify_mtcs_auto(24; skip_FR = true)
    classified = auto.classified

returns a `Vector{ClassifiedMTC}`, one per Galois sector per stratum.

Module organisation (all at ACMG top level — no submodules):
- Core:            shared data types and finite-field arithmetic
- Cyclotomics:     exact `Q(ζ_N)` context, Galois/Frobenius action,
                   and finite-field reduction hooks
- ModularData:     target home for S/T modular-data validation and
                   Verlinde extraction
- SL2:             SL(2, ℤ/N) irreducible representation catalog
- Search:          stratum enumeration and finite-field block-U search
- Reconstruction:  CRT, Galois-aware grouping, and exact cyclotomic lift
- FR:              exact F/R equation generation and solvers
- Gauge:           gauge representatives, transforms, and fixing helpers
- Pipeline:        conductor-first orchestration and output records
- IO:              repository-facing serialization and reporting hooks
- Experimental:    incubating code outside stable public APIs
"""
module ACMG

using LinearAlgebra
using Primes

# Core: shared types and finite-field arithmetic.
include("Core/Types.jl")
include("Core/FpArith.jl")

# Cyclotomics: exact cyclotomic context, reduction, and Galois actions.
include("Cyclotomics/Context.jl")
include("Cyclotomics/Reduction.jl")

# ModularData: exact containers, finite-field S/T validation, and Verlinde extraction.
include("ModularData/Exact.jl")
include("ModularData/Reduction.jl")
include("Cyclotomics/Galois.jl")
include("Cyclotomics/Frobenius.jl")
include("ModularData/Validation.jl")
include("ModularData/Verlinde.jl")
include("ModularData/ExactValidation.jl")

# SL2: SL(2, ℤ/N) irreducible representation catalog.
include("SL2/SL2Reps.jl")

# Search: stratum enumeration and block-U finite-field search.
include("Search/StratumEnum.jl")
include("Search/BlockU.jl")

# Reconstruction: multi-prime CRT reconstruction and exact lifting.
include("Reconstruction/CRT.jl")
include("Reconstruction/CyclotomicCRT.jl")

# FR: exact F/R equations and solvers over cyclotomic fields.
include("FR/PentagonEquations.jl")
include("FR/HexagonEquations.jl")
include("FR/ExactPolynomialSolver.jl")
include("FR/PentagonSolver.jl")
include("FR/HexagonSolver.jl")
include("Reconstruction/ModularDataLift.jl")

# Gauge: public gauge API surface and gauge-fixing helpers.
include("Gauge/Gauge.jl")

# Pipeline: result records, prime selection, FR layer, and conductor-first orchestration.
include("Pipeline/Types.jl")
include("Pipeline/PrimeSelection.jl")
include("Pipeline/Strategy.jl")
include("Pipeline/Auto.jl")
include("Pipeline/FRLayer.jl")
include("Pipeline/Pipeline.jl")

# Gauss sums: exact cyclotomic Gauss sums and higher central charges.
include("GaussSums/GaussSums.jl")
include("GaussSums/HigherCentralCharge.jl")

# Experimental: finite-field F/R prototypes.
include("Experimental/FiniteFieldHigherCentralCharge.jl")

# IO: JSON export/import and Markdown reports for classification outputs.
include("IO/Serialization.jl")

# ============================================================
#  Exports
# ============================================================

# Core types and F_p arithmetic
export ModularDatumFp, FusionRule
export CyclotomicContext, ModularData
export field, zeta, conductor, cond_S, cond_T, cond_F
export semion_modular_data, fibonacci_modular_data, ising_modular_data
export toric_code_modular_data, modular_data
export galois_action, galois_orbit, frobenius, reduce_mod_p
export quantum_dimensions
export total_quantum_dimension_squared, gauss_sum_plus, gauss_sum_minus
export normalized_gauss_sum, HigherCentralChargeResult
export higher_central_charge, higher_central_charges, central_charge
export FRSolutionModP, HigherCentralChargeModPResult
export solve_FR_mod_p, lift_higher_central_charge
export validate_modular_data, build_modular_datum, compute_alpha, compute_charge_conjugation
export check_modular_relations, check_unitarity, check_verlinde_integrality
export check_twist_balance, check_vafa_constraints, check_galois_symmetry
export validate_exact_modular_data, validate_exact_mtc
export fusion_automorphisms, is_fusion_automorphism
export modular_data_automorphisms, is_modular_data_automorphism
export galois_anyon_action, galois_anyon_orbits
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
export passes_unit_axiom, solve_cayley_unit_filtered_blocks
export cayley_so_n, inverse_mod_p, enumerate_so_n_Fp, enumerate_o_n_Fp
export enumerate_o_n_Fp_groebner, enumerate_block_candidates, apply_block_U
export build_verlinde_unit_equations
export build_cayley_link_equations
export is_orthogonal_mod_p
export validate_search_mode
export find_mtcs_at_prime, signed_Fp

# Phase 3: CRT reconstruction
export acmg_crt, crt2, rational_reconstruct
export fusion_signature, canonical_rule, group_mtcs_by_fusion
export reconstruct_rational, reconstruct_cyclotomic_element_from_residues

# Phase 4: (F, R) classification
export get_pentagon_system
export solve_pentagon_modular_crt, solve_pentagon_homotopy, solve_pentagon_newton
export assign_F_to_associator!, hexagon_equations, get_hexagon_system
export number_of_variables_in_hexagon_equations
export solve_hexagon_modular_crt, solve_hexagon_homotopy
export DiscreteLogTable, lift_T_Fp_to_cyclotomic
export canonical_gauge, gauge_equivalent, gauge_transform

# End-to-end pipeline
export ClassifiedMTC, FRRoundtripReport, FRStatus
export FRSkipped, FRSolved, FRNoSolutionFound, FRTimeoutLikeFailure
export FRReconstructionFailed, FRVerificationFailed, fr_status
export select_admissible_primes
export estimate_search_complexity, estimate_phase4_complexity
export recommend_primes, recommend_skip_FR
export compute_FR_from_ST, classify_from_group
export classify_mtcs_at_conductor, classify_mtcs_auto
export save_classification, load_classification
export export_modular_data, export_fusion_rule, export_FR, write_report
export conductor_sanity_table, conductor_sanity_markdown

end # module ACMG
