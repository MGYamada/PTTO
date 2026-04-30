# ACMG.jl v0.8.6 Documentation Coverage

This file records the documentation hardening inventory for v0.8.6.  It is a
maintenance log rather than a rendered manual page.

## Scope

v0.8.6 is documentation-focused.  No new public API or intentional functional
change is introduced.  The work concentrates on exported docstrings,
Documenter pages, README orientation, runnable examples, and clearer stable vs
experimental boundaries.

## Export Inventory

`src/ACMG.jl` exports one top-level API surface; ACMG currently has no public
submodules.  Exported names are grouped below by stability.

| Area | Stable public API | Semi-public low-level API | Experimental/internal exported API | v0.8.6 notes |
| --- | --- | --- | --- | --- |
| Core and cyclotomics | `CyclotomicContext`, `field`, `zeta`, `conductor`, `galois_action`, `galois_orbit`, `frobenius`, `reduce_mod_p` | finite-field arithmetic helpers such as `matmul_mod`, `sqrt_mod`, `root_of_unity` | root-selection and coercion internals when used outside documented inputs | conductor-first semantics documented |
| Modular data | `ModularData`, built-in constructors, `validate_modular_data`, `validate_exact_modular_data`, Verlinde checks | automorphism and extraction helpers | low-level fixed-stratum lift helpers | built-ins and conductor metadata documented |
| Gauss sums | `gauss_sum_plus`, `gauss_sum_minus`, `normalized_gauss_sum`, `central_charge`, `higher_central_charge`, `higher_central_charges`, `HigherCentralChargeResult` | none | finite-field higher-central-charge prototypes | exact normalization caveats documented |
| SL2/Search/Reconstruction | `classify_mtcs_at_conductor`, `classify_mtcs_auto`, prime recommendation helpers | `AtomicIrrep`, `Stratum`, stratum/block-U helpers | most low-level block-U enumeration, CRT reconstruction hooks, Galois grouping internals | stable entry points separated from search machinery |
| F/R infrastructure | `FusionRule`, `FRData`, `F_symbol`, `R_symbol`, fusion-rule constructors, `fr_equation_system`, `pentagon_equations`, `hexagon_equations`, `validate_fr_system` | equation containers and variable records | general finite-field solving and cyclotomic reconstruction | multiplicity-free scope documented |
| Gauge | `GaugeParameters`, `GaugeChoice`, `GaugeFixingResult`, `gauge_parameters`, `gauge_fixing_plan`, `is_gauge_fixed`, `gauge_fix` | `GaugeTransform`, scalar transform helpers | canonical low-level normal-form internals, finite-field stacky weights | safe/unit-channel gauge role documented |
| Braid representations | `FusionPath`, `FusionTreeBasis`, `braid_representation`, `braid_generator`, `braid_generators`, `check_braid_relations` | exact fusion-space basis helpers | finite-field braid reduction and image diagnostics | exact construction separated from diagnostics |
| Zariski diagnostics | none as theorem-level stable API | none | `finite_group_diagnostics`, `generated_subgroup`, `generated_matrix_algebra`, `commutant`, `zariski_closure_diagnostics` | diagnostic-only status documented |
| IO/export | `save_classification`, `load_classification`, `export_modular_data`, `export_fusion_rule`, `export_FR`, `write_report` | sanity-table helpers | JSON schema internals | textual exact-value serialization explained |

## Docstring Status

Before v0.8.6 hardening, `Docs.hasdoc(Docs.Binding(ACMG, sym))` reported 91
exported symbols without direct docstrings.  v0.8.6 adds documentation
attachments for the remaining exported surface, including compact docstrings
for low-level compatibility exports.  The local inventory command now reports
zero missing exported doc bindings.

## Existing Manual Coverage

| Page | Covered before v0.8.6 | v0.8.6 hardening |
| --- | --- | --- |
| `docs/src/index.md` | project overview | linked to stable/experimental boundaries |
| `docs/src/getting_started.md` | quick examples | kept lightweight |
| `docs/src/concepts.md` | conductor-first sketch | expanded to concept-page structure |
| `docs/src/modular_data.md` | modular data overview | clarified exact cyclotomic metadata |
| `docs/src/finite_fields.md` | finite-field notes | clarified shared conductor context |
| `docs/src/fr_symbols.md` | F/R overview | expanded multiplicity-free equation layer |
| `docs/src/gauge_fixing.md` | gauge overview | clarified solver preprocessing and normalization role |
| `docs/src/braid_representations.md` | braid construction | separated exact representation from diagnostics |
| `docs/src/search_blocku.md` | Block-U overview | marked low-level helpers experimental |
| `docs/src/zariski_diagnostics.md` | diagnostic caveats | strengthened non-proof language |
| `docs/src/api_stability.md` | stable/experimental policy | refreshed for v0.8.6 |
| `docs/src/api_reference.md` | hand-picked `@docs` blocks | supplemented by `docs/src/api.md` `@autodocs` |

## Remaining Intentional Gaps

- Low-level search and reconstruction helpers remain exported for v0.8
  compatibility, but are documented as semi-public or experimental rather than
  promoted as stable API.
- Zariski and finite-field image routines remain diagnostics.  Documentation
  avoids presenting them as classification proofs.
- General finite-field F/R solving and cyclotomic reconstruction are interface
  placeholders or experimental paths; examples avoid depending on them as
  mathematical certificates.
