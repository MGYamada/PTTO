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
include("Core/ExperimentalWarnings.jl")

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
include("FR/EquationInfrastructure.jl")
include("FR/PentagonEquations.jl")
include("FR/HexagonEquations.jl")
include("FR/ExactPolynomialSolver.jl")
include("FR/PentagonSolver.jl")
include("FR/HexagonSolver.jl")
include("FR/FRData.jl")
include("FR/FiniteFieldFRData.jl")
include("Reconstruction/ModularDataLift.jl")

# Braid representations consume FRData from the FR layer.
include("BraidRepresentations/BraidRepresentations.jl")

# Gauge: public gauge API surface and gauge-fixing helpers over FRData accessors.
include("Gauge/GaugeWeights.jl")
include("Gauge/GaugeTypes.jl")
include("Gauge/ToricGauge.jl")
include("Gauge/FiniteFieldGauge.jl")
include("Gauge/Gauge.jl")
include("Gauge/GaugeActions.jl")
include("Gauge/GaugeConstraints.jl")
include("Gauge/GaugeValidation.jl")
include("Gauge/GaugeNormalForms.jl")
include("Gauge/Stabilizers.jl")

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

# Experimental finite-field Zariski diagnostics.
include("ZariskiDiagnostics/ZariskiDiagnostics.jl")

# IO: JSON export/import and Markdown reports for classification outputs.
include("IO/Serialization.jl")

# Documentation-only attachments for exported APIs whose definitions live in
# compact implementation files.
include("Documentation.jl")

# ============================================================
#  Public API
# ============================================================
#
# Exported names are the thin public surface intended for ordinary docs and
# examples.  Julia 1.11 `public` declarations below are reserved for stable
# names that should remain qualified as `ACMG.name`.

public HigherCentralChargeResult, HCCGeneratingFunction, HCCLocalFactor
public higher_central_charge_result
public FRStatus, FRSkipped, FRSolved, FRNoSolutionFound
public FRTimeoutLikeFailure, FRReconstructionFailed, FRVerificationFailed

export CyclotomicContext, ModularData, FusionRule, FRData
export BraidRepresentation, ClassifiedMTC

export field, zeta, conductor, cond_S, cond_T, cond_F
export modular_data
export semion_modular_data, fibonacci_modular_data, ising_modular_data
export toric_code_modular_data

export validate_exact_modular_data, validate_modular_data
export check_modular_relations, check_unitarity, check_verlinde_integrality
export check_twist_balance, check_vafa_constraints, check_galois_symmetry
export quantum_dimensions, total_quantum_dimension_squared
export extract_fusion_rule_Fp, lift_fusion_to_Z, extract_and_lift
export verlinde_coefficient

export galois_action, galois_orbit, frobenius, reduce_mod_p
export modular_data_automorphisms, is_modular_data_automorphism
export galois_anyon_action, galois_anyon_orbits

export gauss_sum_plus, gauss_sum_minus, normalized_gauss_sum, central_charge
export higher_central_charge, higher_central_charges
export higher_central_charge_period, higher_central_charge_sequence
export higher_central_charge_generating_function
export higher_central_charge_modp, higher_central_charge_sequence_modp
export hcc_local_factor, hcc_local_factors

export simple_objects, fusion_product, fusion_channels, is_admissible
export simples, fusion_coeff, hom_basis
export F_symbol, R_symbol, has_F_symbol, has_R_symbol
export fusion_rule, F_values, R_values, R_inverse_values, fr_metadata
export fr_scalar_type, fr_value_one, frdata_from_vectors, frdata_from_namedtuple
export is_multiplicity_free, require_multiplicity_free
export semion_fusion_rules, fibonacci_fusion_rules, toric_code_fusion_rules
export ising_fusion_rules

export GaugeAction, GaugeParameters, GaugeChoice, GaugeFixingResult
export identity_gauge, apply_gauge, compose_gauge, inverse_gauge
export gauge_normal_form, validate_gauge_fixed
export gauge_parameters, gauge_fixing_plan, is_gauge_fixed, gauge_fix
export StabilizerProblem, StabilizerEquations, StabilizerResult
export stabilizer, stabilizer_equations, stabilizer_order
export automorphisms, is_trivial_stabilizer, stabilizer_metadata

export FusionPath, FusionTreeBasis
export fusion_paths, fusion_basis, dim
export braid_representation, braid_generator, braid_generators
export braid_generators_B3, check_braid_relations

export classify_mtcs_at_conductor, classify_mtcs_auto
export fr_status
export recommend_primes, recommend_skip_FR
export save_classification, load_classification
export export_modular_data, export_fusion_rule, export_FR, write_report
export conductor_sanity_table, conductor_sanity_markdown

# Experimental entry points used by docs/examples.  They are exported for
# ergonomics but remain outside the stable surface described in api_stability.
export solve_fr_mod_p, semion_fr_data_mod_p, fibonacci_fr_data_mod_p
export ising_fr_data_mod_p, verify_FRData
export generated_subgroup, finite_group_diagnostics
export generated_matrix_algebra, commutant, zariski_closure_diagnostics

end # module ACMG
