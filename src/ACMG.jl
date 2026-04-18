"""
    ACMG — Arithmetic Condensed Matter Geometry

Modular Tensor Category classification via NRWW block-U enumeration
over Z[ζ_N], with F_p filtering/validation layer.

Design philosophy (v0.2):
- Primary search: Z[ζ_N], using SL(2, ℤ/N) irrep catalog + block-U structure
- Validation: F_p layer (from v0.1), already validated against Fibonacci & Ising
- Fusion rule, F, R symbols: derived via Verlinde (F_p) then lifted

Module organization:
- FpArith:       F_p arithmetic primitives
- Types:         ModularDatumFp, FusionRule with axiom validation
- Dimensions:    Quantum dimension enumeration (stub)
- ModularData:   (S, T) axiom checking over F_p
- FusionExtract: Verlinde-based fusion rule extraction
- Enumerator:    top-level driver (stub)
- SL2Reps:       SL(2, ℤ/N) irreducible representation catalog (via GAP/SL2Reps)
"""
module ACMG

using LinearAlgebra
using Primes

# Core arithmetic layer
include("FpArith.jl")

# Modular data and fusion rule types
include("Types.jl")

# Quantum dimension enumeration
include("Dimensions.jl")

# (S, T) enumeration over F_p
include("ModularData.jl")

# Fusion rule extraction
include("FusionExtract.jl")

# Top-level enumeration driver
include("Enumerator.jl")

# SL(2, ℤ/N) irrep catalog (Oscar + GAP/SL2Reps)
include("SL2Reps.jl")

# Stratum enumeration (combinatorial partition of rank by irrep dimensions)
include("StratumEnum.jl")

# Block-U parametrisation and MTC reconstruction (Phase 2)
include("BlockU.jl")

# Multi-prime CRT reconstruction (Phase 3)
include("CRT.jl")

# Phase 4: Pentagon/Hexagon solver for (F, R) classification
include("Phase4/Phase4.jl")

# Exports
export ModularDatumFp, FusionRule
export validate_modular_data, build_modular_datum, compute_alpha, compute_charge_conjugation
export extract_fusion_rule_Fp, lift_fusion_to_Z, extract_and_lift
export verlinde_coefficient
export is_square, sqrt_mod, primitive_root, root_of_unity, roots_of_unity
export matmul_mod, matpow_mod, diagmul_right, diagmul_left, lift_symmetric
export fusion_isomorphic, fusion_matrix, validate
export AtomicIrrep, build_atomic_catalog, all_divisors
export Stratum, enumerate_strata, count_strata, describe_stratum, find_unit_indices
export MTCCandidate, build_block_diagonal, reduce_matrix_to_Fp, reduce_vector_to_Fp
export find_zeta_in_Fp, cyclotomic_to_Fp
export t_eigenspace_decomposition, parameter_dim
export o2_circle_points, apply_o2_block, verlinde_find_unit
export cayley_so_n, inverse_mod_p, enumerate_so_n_Fp, enumerate_o_n_Fp, apply_block_U
export find_mtcs_at_prime, signed_Fp
# Phase 3: CRT reconstruction
export acmg_crt, crt2, rational_reconstruct, compute_sqrt_d_mod_p
export compute_sqrt3_cyclotomic_mod_p, compute_sqrt2_cyclotomic_mod_p
export fusion_signature, group_mtcs_by_fusion, group_mtcs_galois_aware
export reconstruct_rational, reconstruct_in_Z_sqrt_d
export reconstruct_matrix_in_Z_sqrt_d, reconstruct_S_matrix
export verify_reconstruction, describe_matrix

# Phase 4 exports (namespace access via ACMG.Phase4.*)
# The Phase4 submodule is accessible as ACMG.Phase4 without re-export.
# Users should call e.g. ACMG.Phase4.get_pentagon_system(...)

end # module ACMG
