# API Stability

ACMG.jl is a research computational framework.  Until v1.0, some APIs may
change as the mathematical and computational interfaces settle.

## Stable/public API

Stable public API means documented and exported functions intended for normal
users.  It also includes the basic APIs used in examples and tests, especially:

- conductor and cyclotomic utilities
- exact modular-data construction and validation
- Verlinde checks and fusion-rule extraction
- Galois action utilities
- higher central charge utilities over exact modular data
- basic F/R data containers and accessors
- basic braid representation constructors and braid-relation checks
- conductor-level classification pipeline entry points

## Experimental API

Experimental APIs may receive breaking changes before v1.0.  Pin the ACMG.jl
version when using them in research scripts or generated data pipelines.

Diagnostics return computational evidence, not mathematical theorems.  In
particular, Zariski and finite-field image diagnostics should not be cited as
proofs of characteristic-zero density, irreducibility, or exact classification
without separate mathematical verification.

### Experimental API list

`solve_fr_mod_p`

- status: public finite-field F/R entry point; the built-in modular-data
  prototype overload remains experimental
- mathematical caveat: finite-field F/R output is a fiber of the
  pentagon/hexagon scheme, not a cyclotomic lift; modular-data prototype
  residues likewise require exact comparison or lifting before use as
  cyclotomic data
- compatibility note: the old uppercase `solve_FR_mod_p` name is deprecated
  and no longer exported

`lift_higher_central_charge`

- status: experimental
- input/output may change before v1.0
- mathematical caveat: placeholder reconstruction helper; it does not prove a
  cyclotomic lift
- stable alternative: compute exact `higher_central_charge` and compare by
  `reduce_mod_p`

`solve_finite_field`

- status: experimental
- input/output may change before v1.0
- mathematical caveat: general finite-field F/R solving is incomplete; future
  solutions require exact verification
- stable alternative: use `fr_equation_system`, exact Phase-4 helpers, or
  built-in exact F/R data

`cyclotomic_reconstruct`

- status: experimental
- input/output may change before v1.0
- mathematical caveat: interface validation only; modular residues do not
  constitute a reconstruction proof
- stable alternative: use exact cyclotomic constructors and validators

`generated_subgroup`

- status: experimental
- input/output may change before v1.0
- mathematical caveat: bounded finite-field enumeration is computational
  evidence and may truncate
- stable alternative: none for full image diagnostics; use
  `check_braid_relations` for stable relation checking

`finite_group_diagnostics`

- status: experimental
- input/output may change before v1.0
- mathematical caveat: finite-field image data and sampled traces are
  diagnostic only
- stable alternative: `check_braid_relations`

`generated_matrix_algebra`

- status: experimental
- input/output may change before v1.0
- mathematical caveat: finite-field matrix algebra dimension is not a full
  characteristic-zero Zariski closure computation
- stable alternative: none

`commutant`

- status: experimental
- input/output may change before v1.0
- mathematical caveat: finite-field commutant data are diagnostic evidence
- stable alternative: none

`zariski_closure_diagnostics`

- status: experimental
- input/output may change before v1.0
- mathematical caveat: does not compute a full Zariski closure or classify
  algebraic groups
- stable alternative: none

Finite-field braid representation helpers

- status: experimental
- input/output may change before v1.0
- mathematical caveat: finite-field reductions require split-prime conditions
  and do not replace exact mathematical verification
- stable alternative: finite-field `FRData` plus `braid_representation` and
  `check_braid_relations` over the selected prime field

GaugeFixing low-level normal form helpers

- status: experimental/internal
- input/output may change before v1.0
- mathematical caveat: low-level normal forms encode implementation choices
  rather than canonical mathematical data
- stable alternative: `GaugeAction`, `apply_gauge`, `identity_gauge`,
  `gauge_normal_form`, `validate_gauge_fixed`, and the compatibility wrappers
  `gauge_parameters`, `gauge_fixing_plan`, `is_gauge_fixed`, and `gauge_fix`

Search/BlockU low-level enumeration helpers

- status: experimental/internal
- input/output may change before v1.0
- mathematical caveat: enumeration order, pruning heuristics, and diagnostic
  metadata are implementation details
- stable alternative: `classify_mtcs_at_conductor`, `classify_mtcs_auto`,
  `find_mtcs_at_prime`

## Internal API

Internal APIs are implementation details and may change without deprecation
before v1.0.  This includes cache structures, polynomial equation internals,
search heuristics, low-level Block-U enumeration state, and private helper
functions whose names begin with `_`.

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
