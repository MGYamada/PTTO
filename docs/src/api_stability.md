# API Stability

ACMG.jl is a research computational framework.  Stable names are deliberately
thin; lower-level search, reconstruction, and diagnostic hooks may change as
the mathematical and computational interfaces settle.

## Stable/public API

Stable public API means documented names with Julia 1.11 public status.
Ordinary docs and examples use exported names via `using ACMG`; stable names
that should stay qualified are marked with the `public` keyword instead.

The stable/public surface is intentionally thin and centered on normal user
workflows:

- conductor and cyclotomic utilities
- exact modular-data construction and validation
- Verlinde checks and fusion-rule extraction
- Galois action utilities
- higher central charge utilities over exact modular data
- basic F/R data containers and accessors
- basic braid representation constructors and braid-relation checks
- conductor-level classification pipeline entry points

### Public names

Cyclotomic and modular data:

- `CyclotomicContext`, `ModularData`, `field`, `zeta`, `conductor`,
  `cond_S`, `cond_T`, `cond_F`
- `modular_data`, `semion_modular_data`, `fibonacci_modular_data`,
  `ising_modular_data`, `toric_code_modular_data`
- `validate_exact_modular_data`, `validate_modular_data`,
  `check_modular_relations`, `check_unitarity`,
  `check_verlinde_integrality`, `check_twist_balance`,
  `check_vafa_constraints`, `check_galois_symmetry`
- `quantum_dimensions`, `total_quantum_dimension_squared`,
  `extract_fusion_rule_Fp`, `lift_fusion_to_Z`, `extract_and_lift`,
  `verlinde_coefficient`

Galois and higher central charge:

- `galois_action`, `galois_orbit`, `modular_data_automorphisms`,
  `is_modular_data_automorphism`, `galois_anyon_action`,
  `galois_anyon_orbits`
- `gauss_sum_plus`, `gauss_sum_minus`, `normalized_gauss_sum`,
  `central_charge`
- `HigherCentralChargeResult`, `HCCGeneratingFunction`, `HCCLocalFactor`,
  `higher_central_charge`, `higher_central_charges`,
  `higher_central_charge_result`, `higher_central_charge_period`,
  `higher_central_charge_sequence`,
  `higher_central_charge_generating_function`

F/R data, gauge, and braids:

- `FusionRule`, `FRData`, `simple_objects`, `fusion_product`,
  `fusion_channels`, `is_admissible`, `simples`, `fusion_coeff`,
  `hom_basis`
- `F_symbol`, `R_symbol`, `has_F_symbol`, `has_R_symbol`, `fusion_rule`,
  `F_values`, `R_values`, `R_inverse_values`, `fr_metadata`,
  `fr_scalar_type`, `fr_value_one`, `frdata_from_vectors`,
  `frdata_from_namedtuple`
- `is_multiplicity_free`, `require_multiplicity_free`,
  `semion_fusion_rules`, `fibonacci_fusion_rules`,
  `toric_code_fusion_rules`, `ising_fusion_rules`
- `GaugeAction`, `GaugeParameters`, `GaugeChoice`, `GaugeFixingResult`,
  `identity_gauge`, `apply_gauge`, `compose_gauge`, `inverse_gauge`,
  `gauge_normal_form`, `validate_gauge_fixed`, `gauge_parameters`,
  `gauge_fixing_plan`, `is_gauge_fixed`, `gauge_fix`
- `FusionPath`, `FusionTreeBasis`, `BraidRepresentation`, `fusion_paths`,
  `fusion_basis`, `dim`, `braid_representation`, `braid_generator`,
  `braid_generators`, `braid_generators_B3`, `check_braid_relations`

Pipeline and IO:

- `ClassifiedMTC`, `classify_mtcs_at_conductor`, `classify_mtcs_auto`,
  `recommend_primes`, `recommend_skip_FR`
- `save_classification`, `load_classification`, `export_modular_data`,
  `export_fusion_rule`, `export_FR`, `write_report`,
  `conductor_sanity_table`, `conductor_sanity_markdown`

## Experimental API

Experimental APIs may receive breaking changes.  Pin the ACMG.jl revision when
using them in research scripts or generated data pipelines.  Some experimental
entry points are exported because they appear in docs or examples, but they are
not part of the stable surface listed above.

Diagnostics return computational evidence, not mathematical theorems.  In
particular, Zariski and finite-field image diagnostics should not be cited as
proofs of characteristic-zero density, irreducibility, or exact classification
without separate mathematical verification.

### Experimental API list

`solve_fr_mod_p`

- status: public finite-field F/R entry point for fusion rules and
  `FRData{FpElem}` construction
- mathematical caveat: finite-field F/R output is a fiber of the
  pentagon/hexagon scheme, not a cyclotomic lift
- compatibility note: the old uppercase `solve_FR_mod_p` name and the
  built-in modular-data prototype overload have been removed

Finite-field higher central charge helpers

- status: experimental for `higher_central_charge_modp`,
  `higher_central_charge_sequence_modp`, `hcc_local_factor`, and
  `hcc_local_factors`
- mathematical caveat: prime-field reductions are implemented; extension-field
  norm factors are future work
- stable alternative: exact `higher_central_charge`, `higher_central_charges`,
  `higher_central_charge_period`, `higher_central_charge_sequence`, and
  `higher_central_charge_generating_function`

`solve_finite_field`

- status: experimental
- input/output may change
- mathematical caveat: general finite-field F/R solving is incomplete; future
  solutions require exact verification
- stable alternative: use `fr_equation_system`, exact F/R reconstruction helpers, or
  built-in exact F/R data

`cyclotomic_reconstruct`

- status: experimental
- input/output may change
- mathematical caveat: interface validation only; modular residues do not
  constitute a reconstruction proof
- stable alternative: use exact cyclotomic constructors and validators

`generated_subgroup`

- status: experimental
- input/output may change
- mathematical caveat: bounded finite-field enumeration is computational
  evidence and may truncate
- stable alternative: none for full image diagnostics; use
  `check_braid_relations` for stable relation checking

`finite_group_diagnostics`

- status: experimental
- input/output may change
- mathematical caveat: finite-field image data and sampled traces are
  diagnostic only
- stable alternative: `check_braid_relations`

`generated_matrix_algebra`

- status: experimental
- input/output may change
- mathematical caveat: finite-field matrix algebra dimension is not a full
  characteristic-zero Zariski closure computation
- stable alternative: none

`commutant`

- status: experimental
- input/output may change
- mathematical caveat: finite-field commutant data are diagnostic evidence
- stable alternative: none

`zariski_closure_diagnostics`

- status: experimental
- input/output may change
- mathematical caveat: does not compute a full Zariski closure or classify
  algebraic groups
- stable alternative: none

Finite-field braid representation helpers

- status: experimental
- input/output may change
- mathematical caveat: finite-field reductions require split-prime conditions
  and do not replace exact mathematical verification
- stable alternative: finite-field `FRData` plus `braid_representation` and
  `check_braid_relations` over the selected prime field

GaugeFixing low-level normal form helpers

- status: experimental/internal
- input/output may change
- mathematical caveat: low-level normal forms encode implementation choices
  rather than canonical mathematical data
- stable alternative: `GaugeAction`, `apply_gauge`, `identity_gauge`,
  `gauge_normal_form`, `validate_gauge_fixed`, and the compatibility wrappers
  `gauge_parameters`, `gauge_fixing_plan`, `is_gauge_fixed`, and `gauge_fix`

Search/BlockU low-level enumeration helpers

- status: experimental/internal
- input/output may change
- mathematical caveat: enumeration order, pruning heuristics, and diagnostic
  metadata are implementation details
- stable alternative: `classify_mtcs_at_conductor`, `classify_mtcs_auto`,
  `find_mtcs_at_prime`

## Internal API

Internal APIs are implementation details and may change without deprecation.
This includes cache structures, polynomial equation internals, search
heuristics, low-level Block-U enumeration state, and private helper functions
whose names begin with `_`.

Some internal or experimental helpers remain exported for backward
compatibility.  They are classified here rather than removed to avoid a
breaking release.

## Experimental warnings

Experimental warnings are opt-in:

```julia
ENV["ACMG_WARN_EXPERIMENTAL"] = "true"
```

When enabled, experimental APIs call `warn_experimental(name)` and warn at
most once per API name in a Julia session.
