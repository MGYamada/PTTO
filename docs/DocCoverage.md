# ACMG.jl Documentation Coverage

This file records the documentation coverage inventory.  It is a maintenance
log rather than a rendered manual page.

## Scope

The documentation work concentrates on public docstrings, Documenter pages,
README orientation, runnable examples, and clearer stable vs experimental
boundaries.

## Public Inventory

`src/ACMG.jl` exposes one top-level module; ACMG currently has no public
submodules.  Names are grouped below by stability.  Ordinary examples use a
thin exported surface; stable qualified-only names are marked with Julia
1.11's `public` keyword.

| Area | Stable public API | Semi-public low-level API | Experimental/internal API | Notes |
| --- | --- | --- | --- | --- |
| Core and cyclotomics | `CyclotomicContext`, `field`, `zeta`, `conductor`, `galois_action`, `galois_orbit`, `frobenius`, `reduce_mod_p` | finite-field arithmetic helpers such as `matmul_mod`, `sqrt_mod`, `root_of_unity` | root-selection and coercion internals when used outside documented inputs | conductor-first semantics documented |
| Modular data | `ModularData`, built-in constructors, `validate_modular_data`, `validate_exact_modular_data`, Verlinde checks | automorphism and extraction helpers | low-level fixed-stratum lift helpers | built-ins and conductor metadata documented |
| Gauss sums | `gauss_sum_plus`, `gauss_sum_minus`, `normalized_gauss_sum`, `central_charge`, `higher_central_charge`, `higher_central_charges`, `higher_central_charge_period`, `higher_central_charge_sequence`, `HCCGeneratingFunction` | `higher_central_charge_result`, `HigherCentralChargeResult` | finite-field HCC reductions and local factors | `D^2` moment normalization documented |
| SL2/Search/Reconstruction | `classify_mtcs_at_conductor`, `classify_mtcs_auto`, prime recommendation helpers | `AtomicIrrep`, `Stratum`, stratum/block-U helpers | most low-level block-U enumeration, CRT reconstruction hooks, Galois grouping internals | stable entry points separated from search machinery |
| F/R infrastructure | `FusionRule`, `FRData`, `F_symbol`, `R_symbol`, fusion-rule constructors, `fr_equation_system`, `pentagon_equations`, `hexagon_equations`, `validate_fr_system` | equation containers and variable records | general finite-field solving and cyclotomic reconstruction | multiplicity-free scope documented |
| Gauge | `GaugeAction`, `GaugeParameters`, `GaugeChoice`, `GaugeFixingResult`, `identity_gauge`, `apply_gauge`, `compose_gauge`, `inverse_gauge`, `gauge_normal_form`, `validate_gauge_fixed`, compatibility wrappers | `GaugeTransform`, scalar transform helpers, explicit constraint records | quotient diagnostics, finite-field stacky weights, future solver-specific normal forms | algebraic gauge-action role documented |
| Braid representations | `FusionPath`, `FusionTreeBasis`, `braid_representation`, `braid_generator`, `braid_generators`, `check_braid_relations` | exact fusion-space basis helpers | finite-field braid reduction and image diagnostics | exact construction separated from diagnostics |
| Zariski diagnostics | none as theorem-level stable API | none | `finite_group_diagnostics`, `generated_subgroup`, `generated_matrix_algebra`, `commutant`, `zariski_closure_diagnostics` | diagnostic-only status documented |
| IO/export | `save_classification`, `load_classification`, `export_modular_data`, `export_fusion_rule`, `export_FR`, `write_report` | sanity-table helpers | JSON schema internals | textual exact-value serialization explained |

## Docstring Status

The docstring inventory checks public bindings with
`Docs.hasdoc(Docs.Binding(ACMG, sym))`.  The public surface should keep direct
documentation attachments, including compact docstrings for low-level
compatibility entry points.

## Existing Manual Coverage

| Page | Original coverage | Current hardening |
| --- | --- | --- |
| `docs/src/index.md` | project overview | linked to stable/experimental boundaries |
| `docs/src/getting_started.md` | quick examples | kept lightweight |
| `docs/src/concepts.md` | conductor-first sketch | expanded to concept-page structure |
| `docs/src/modular_data.md` | modular data overview | clarified exact cyclotomic metadata |
| `docs/src/finite_fields.md` | finite-field notes | clarified shared conductor context |
| `docs/src/fr_symbols.md` | F/R overview | expanded multiplicity-free equation layer |
| `docs/src/gauge_fixing.md` | gauge overview | expanded to algebraic action, constraints, normal forms, finite-field relation, and migration notes |
| `docs/src/braid_representations.md` | braid construction | separated exact representation from diagnostics |
| `docs/src/search_blocku.md` | Block-U overview | marked low-level helpers experimental |
| `docs/src/zariski_diagnostics.md` | diagnostic caveats | strengthened non-proof language |
| `docs/src/api_stability.md` | stable/experimental policy | refreshed for current API boundaries |
| `docs/src/api_reference.md` | hand-picked `@docs` blocks | supplemented by `docs/src/api.md` `@autodocs` |

## Remaining Intentional Gaps

- Low-level search and reconstruction helpers remain implementation details or
  experimental entry points rather than promoted stable API.
- Zariski and finite-field image routines remain diagnostics.  Documentation
  avoids presenting them as classification proofs.
- General finite-field F/R solving and cyclotomic reconstruction are interface
  placeholders or experimental paths; examples avoid depending on them as
  mathematical certificates.
