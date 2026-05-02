# Changelog

## v0.9.1 - Gauge Stabilizers and API Boundary Refinements

### Added
- Added experimental gauge stabilizer infrastructure for computing
  `Aut(x) = Stab_G(x)` for already available solution data:
  `StabilizerProblem`, `StabilizerEquations`, `StabilizerResult`,
  `stabilizer`, `stabilizer_equations`, `stabilizer_order`,
  `automorphisms`, `is_trivial_stabilizer`, and `stabilizer_metadata`.
- Added brute-force stabilizer enumeration for small finite gauge groups,
  with exact structural comparison and no floating-point conversion.
- Added finite-field toric gauge-group support for multiplicity-free
  `FRData{FpElem}` and deterministic character-equation generation for
  active F/R coordinates.
- Added `docs/src/gauge_stabilizers.md`, explaining gauge stabilizers,
  their role as future stacky denominators `1 / |Aut(x)|`, and the v0.9.1
  non-goals around zeta functions, Lefschetz sums, trace formulas, and
  quotient-stack classification.
- Added `examples/15_gauge_stabilizers.jl`, a small finite-field toy example
  computing stabilizer orders and future stacky weights.
- Added stabilizer tests covering the toy action `F_q^*` on `F_q`,
  problem/result accessors, equation generation, identity-gauge behavior, and
  explicit errors for too-large finite gauge groups.

### Changed
- Bumped the package version to `0.9.1`.
- Exported `fr_status` so documentation examples can use
  `fr_status(m)` directly.
- Marked `FRStatus`, `FRSkipped`, `FRSolved`, `FRNoSolutionFound`,
  `FRTimeoutLikeFailure`, `FRReconstructionFailed`, and
  `FRVerificationFailed` as Julia `public` names.
- Aligned `central_charge(data)` with the D²-normalized convention used by
  `higher_central_charge(data, 1)`.  The legacy Galois-normalized structured
  result remains available through
  `higher_central_charge_result(data; normalization = :galois)`.
- Clarified that `solve_fr_mod_p` is an exported experimental finite-field
  F/R entry point rather than part of the stable public surface.
- Added the gauge stabilizer page to the Documenter navigation and updated
  API stability docs/tests for the new experimental exports.

### Non-goals
- v0.9.1 does not implement zeta functions, Euler factors, Lefschetz trace
  sums, higher-central-charge trace sums, étale cohomology, full stacky point
  counting, or automatic quotient-stack classification.

### Tests
- Verified the full package test suite with
  `julia --project=. -e 'using Pkg; Pkg.test()'`.
- Verified the documentation build with
  `julia --project=docs -e 'using Pkg; Pkg.instantiate(); include("docs/make.jl")'`.

## v0.9.0 - Higher Central Charges, Public API Boundary, and Classification Docs

### Added
- Added a dedicated higher-central-charge API using the canonical
  D²-normalized moment convention:
  `higher_central_charge`, `higher_central_charges`,
  `higher_central_charge_period`, `higher_central_charge_sequence`, and
  `higher_central_charge_generating_function`.
- Added `HCCGeneratingFunction` as a structured exact generating-function
  object storing twist-spectrum weights and twists.
- Added experimental finite-field higher-central-charge helpers:
  `higher_central_charge_modp`, `higher_central_charge_sequence_modp`,
  `hcc_local_factor`, and `hcc_local_factors`.
- Added regression coverage for Semion, Fibonacci, and Ising higher central
  charges, including the Fibonacci `p = 11`, `ζ_5 ↦ 3` sequence
  `[1, 6, 7, 5, 9]` and local factor checks.
- Added `docs/src/invariants/higher_central_charge.md` with the exact moment
  definition, examples, finite-field reduction notes, and local Euler factors.
- Added `docs/src/classify_mtcs.md`, a mathematically focused classification
  guide explaining conductor-first search, `SL(2, Z/N)` representation strata,
  Block-U search, exact reconstruction, and F/R reconstruction.
- Added API stability tests for the thin exported/public surface under
  Julia 1.11.

### Changed
- Updated `Project.toml` compatibility to require Julia `1.11` and
  `LinearAlgebra` `1.11`, enabling use of Julia's `public` keyword.
- Reworked the top-level API boundary:
  docs/examples-level names are exported for `using ACMG`, while qualified-only
  stable records such as HCC helper records are marked `public` without export.
- Expanded docstrings for classification records, strata, conductor-first
  pipeline entry points, SL2 atomic catalogs, exact F/R reconstruction helpers,
  finite-field diagnostics, and higher-central-charge utilities.
- Removed version-specific prose from the docs and public-facing docstrings.
- Replaced internal stage terminology in docs, docstrings, error messages,
  verbose output, and tests with mathematical or functional names such as
  Block-U search, CRT modular-data reconstruction, and exact F/R reconstruction.
- Renamed the experimental finite-field F/R reference-solution selector from
  `solver = :phase4` to `solver = :reference`.
- Renamed the coarse F/R cost helper to
  `estimate_fr_reconstruction_complexity`.
- Refined `classify_mtcs_auto` and `classify_mtcs_at_conductor` docs to state
  that the input conductor is fixed and the modular representation factors
  through `SL(2, Z/N)`.
- Updated `docs/src/api_stability.md` to describe the thin stable surface,
  exported example-level names, qualified-only public records, experimental
  APIs, and current limitations.

### Removed
- Removed duplicate higher-central-charge finite-field prototype plumbing from
  `Experimental/FiniteFieldHigherCentralCharge.jl`.
- Removed the old uppercase `solve_FR_mod_p` public prototype path and the
  obsolete `FRSolutionModP`, `HigherCentralChargeModPResult`, and
  `lift_higher_central_charge` experimental records.
- Removed public exposure of internals from
  `Experimental/FiniteFieldHigherCentralCharge.jl`; finite-field HCC behavior
  now routes through the main HCC API.

### Compatibility
- `higher_central_charge_result` remains available as the structured legacy
  Gauss-sum compatibility result.
- Exact F/R reconstruction and finite-field F/R APIs remain experimental
  outside the documented stable workflow boundary.
- `conductor_mode` accepts only `:full_mtc`; unsupported modes now report a
  direct unsupported-mode error.
- Low-level Block-U, reconstruction, and diagnostic helpers remain outside the
  stable public API even when reachable for research workflows.

### Documentation
- Added the HCC documentation page to the manual navigation.
- Added the classification guide to the manual navigation and linked the
  lower-level Block-U page to it.
- Improved the home page documentation map and stability guidance.
- Kept examples in the `using ACMG` style so documented workflows match the
  exported API surface.

### Tests
- Verified the full package test suite with `julia --project=. test/runtests.jl`.
- Verified the documentation build with `julia --project=docs docs/make.jl`.

## v0.8.9 - Algebraic Gauge Fixing Infrastructure

### Added
- Added `GaugeAction` as the FRData-centered gauge action record, with
  `identity_gauge`, `apply_gauge`, `compose_gauge`, `inverse_gauge`, and
  `validate_gauge_action`.
- Added explicit gauge-fixing constraint records:
  `GaugeConstraint`, `FixUnitConstraints`, `FixSelectedFSymbols`,
  `FixSelectedRSymbols`, and `NormalizationConstraint`.
- Added staged normal-form helpers:
  `gauge_degrees_of_freedom`, `build_gauge_constraints`,
  `solve_gauge_constraints`, `gauge_normal_form`, and
  `validate_gauge_fixed`.
- Added FRData equality on the mathematical F/R payload, ignoring metadata, so
  identity and inverse gauge-action tests can compare representatives directly.
- Expanded gauge-fixing documentation with the algebraic action formula,
  finite-field notes, examples, and migration guidance.

### Changed
- Split gauge infrastructure into focused files for type definitions, actions,
  constraints, validation, and normal forms while preserving legacy wrappers.
- Corrected the R-symbol gauge character to match the TensorCategories
  hexagon convention and explicit inverse-R variables:
  `R'^{ab}_c = u[a,b,c] / u[b,a,c] R^{ab}_c`.
- Updated finite-field and toric gauge helpers to consume `GaugeAction`.
- Gauge tests now avoid invoking the expensive `ising_fr_data_mod_p(17)`
  solver path; the Ising hexagon-preservation regression uses a precomputed
  finite-field fixture.
- Fixed finite-field F/R solver issues carried forward from v0.8.8:
  `solver=:phase4` now accepts Phase-4 keyword arguments, empty-F pentagon
  verification no longer throws, known-family branch selection checks target
  twists, and Ising finite-field fixture tests are included in the default
  test suite without invoking the expensive solver path.
- Unified finite-field solver naming on `solve_fr_mod_p`; the old uppercase
  `solve_FR_mod_p` compatibility wrapper is deprecated and no longer exported.
- Implemented equation-layer `gauge_fix` end to end: selected F-symbol pivot
  coordinates are substituted by `1`, solver variables are reduced, and
  finite-field solutions are expanded back to TensorCategories F/R order.
- Added Smith-normal-form F-symbol toric gauge slices before finite-field
  F/R solving, with stabilizer and stacky metadata reported on solved
  `FRData`.  R-symbols are not gauge-fixed by this pre-solver slice.

### Compatibility
- Legacy `GaugeTransform`, `GaugeParameters`, `GaugeChoice`,
  `gauge_transform`, `canonical_gauge`, `gauge_equivalent`,
  `gauge_fixing_plan`, and `is_gauge_fixed` remain available.
- Equation-layer `gauge_fix(system; strategy=:safe)` keeps the public entry
  point while now returning a genuinely reduced equation system with metadata
  needed to recover the original coordinate order.

## v0.8.8 - Finite-Field FRData Source of Truth

### Added
- Added `FpElem` and `FRData{FpElem}` finite-field F/R data produced by
  `solve_fr_mod_p(rules, p)`.
- Added `semion_fr_data_mod_p`, `fibonacci_fr_data_mod_p`, and
  `ising_fr_data_mod_p` convenience constructors that solve Phase-4
  pentagon/hexagon fibers over `F_p`.
- Added finite-field F/R verification helpers: `verify_pentagon`,
  `verify_hexagon`, and `verify_FRData`.
- Added finite-field braid tests for Semion, Fibonacci, and Ising, including
  B3 braid-relation checks over `F_p`.
- Added finite-field FRData documentation.

### Changed
- Braid construction now consumes finite-field `FRData{FpElem}` directly.
- Removed legacy hard-coded small-example F/R constructors; new braid examples
  use solved finite-field FRData and integer simple-object labels.

## v0.8.6 - Documentation Hardening

### Added
- Added a documentation coverage inventory at `docs/DocCoverage.md`.
- Added an `@autodocs` API reference page at `docs/src/api.md`.
- Added runnable documentation examples for conductor-first modular data,
  F/R-gauge-braid usage, and export helpers.

### Changed
- Expanded docstrings for exported public API groups across cyclotomics,
  modular data, higher central charges, F/R infrastructure, gauge helpers,
  braid representations, pipeline records, diagnostics, and IO/export helpers.
- Reworked concept pages around the standard structure: scope, minimal
  example, mathematical meaning, API overview, common pitfalls, and stability
  notes.
- Shortened README and moved detailed explanations into the manual pages.
- Clarified stable, semi-public, experimental, and internal API boundaries.

### Compatibility
- No intentional functional changes.
- Existing exported names, type names, and return values are unchanged.

### Tests
- Documentation examples were run locally with `julia --project=.`
  invocations.
- Full package tests and docs build were checked for this documentation-only
  release.

## v0.8.5 - API Stability Documentation and Experimental API Marking

### Added
- Added a Documenter-based docs tree under `docs/src/` with pages for getting
  started, concepts, modular data, F/R symbols, gauge fixing, braid
  representations, finite fields, Block-U search, Zariski diagnostics, API
  stability, and API reference.
- Added `docs/src/api_stability.md`, defining stable/public, experimental,
  and internal API boundaries for pre-1.0 releases.
- Added an experimental API list covering finite-field F/R prototypes,
  cyclotomic reconstruction hooks, finite-field braid diagnostics, Zariski
  diagnostics, gauge-fixing internals, and low-level Block-U search helpers.
- Added `warn_experimental(name)` as an opt-in warning helper.  Warnings are
  emitted only when `ENV["ACMG_WARN_EXPERIMENTAL"] == "true"` and are
  de-duplicated per API name.
- Added GitHub Actions CI with separate test and docs build jobs.
- Added tests for API stability docs existence, experimental warning behavior,
  and documented experimental symbol availability.

### Changed
- Reorganized legacy docs pages into the new `docs/src/` structure.
- Updated README with the v0.8.5 research focus, API stability policy link,
  and experimental API warning.
- Added "Experimental API" docstrings and mathematical caveats to selected
  exported experimental functions without changing their behavior.
- Classified potentially internal exported helpers in docs/TODO comments
  instead of removing exports, preserving v0.8 compatibility.

### Removed
- Removed legacy top-level docs pages after folding their content into the new
  docs tree: `docs/BraidRepresentations.md`, `docs/FRInfrastructure.md`, and
  `docs/HigherCentralCharges.md`.

### Tests
- Verified the full test suite with
  `julia --project=. -e 'using Pkg; Pkg.test()'`.
- Verified docs build with `julia --project=docs docs/make.jl`.

## v0.8.4 - Block-U Refactor and Strict Phase 4 Roundtrips

### Changed
- Split the large `src/Search/BlockU.jl` implementation into focused files
  under `src/Search/BlockU/` for types, finite-field helpers, block assembly,
  constraints, enumeration, canonicalization, and search orchestration.
- Kept `src/Search/BlockU.jl` as a thin aggregator so existing
  `include("Search/BlockU.jl")` users continue to work.
- Tightened exact `(F, R)` attachment: Phase 4 candidates are now attachable
  only after strict exact `(S, T)` roundtrip verification.
- When Phase 3 carries a projective or SL2-representation `T` branch that
  disagrees with a valid exact `(F, R)` solution, the pipeline now replaces it
  with the exact `T` reconstructed from the selected braiding data and then
  rechecks strict roundtrip consistency.

### Fixed
- Fixed the built-in Ising reference `FRData` braiding signs so it roundtrips
  to the target Ising `T` matrix.
- Fixed conductor-16 Ising classification so rank-3 Ising sectors are retained
  with exact `(F, R) roundtrip=✓` instead of being returned as `S✓/T✗` or
  discarded.

### Tests
- Added coverage for the Ising reference `FRData` exact `T` roundtrip.
- Verified the full test suite with `julia --project=. -e 'using Pkg; Pkg.test()'`.

## v0.8.3 - FRData-Centered Gauge Refactor

### Added
- `src/FR/FRData.jl` as the canonical home for `FRData` and F/R accessors.
- FRData accessors: `simples`, `fusion_coeff`, `hom_basis`, `F_symbol`,
  `R_symbol`, `has_F_symbol`, `has_R_symbol`, `gauge_basis_indices`, and
  `validate_frdata_for_gauge`.
- FRData-backed gauge APIs: `gauge_transform(data::FRData, ...)`,
  `canonical_gauge(data::FRData)`, and `gauge_equivalent(data1, data2)`.
- Typed gauge containers `GaugeParameters`, `GaugeChoice`, and
  `GaugeFixingResult`.
- Gauge tests covering object-label access, Hom-basis indices, identity
  transforms, scalar transforms, typed gauge parameters, and finite-field
  gauge application.

### Changed
- Gauge internals now consume fusion/object data through FRData/FusionRule
  accessors instead of assuming ad hoc F/R dictionary shapes.
- Gauge-specific FRData validation is owned by the Gauge layer, while FRData
  owns storage and F/R/object/Hom accessors.
- Braid representations now consume `FRData` from the FR layer rather than
  defining it locally.
- `FRData` is now parameterized as `FRData{T}` over the F/R scalar type, and
  `BraidRepresentation{T,M}` stores concrete generator matrix vectors.
- Built-in Semion, Fibonacci, and Ising FRData include simple-object labels
  for label-stable accessor calls and concrete scalar vectors instead of
  `Vector{Any}`.
- Phase 4 candidate selection now prioritizes exact S-roundtrip matches before
  comparing T-branch errors, and the conductor pipeline keeps exact `(F, R)`
  data when only the projective T branch disagrees with the selected target.

### Compatibility
- Legacy vector APIs such as `gauge_transform(F, R, gauge; Nijk=...)`,
  `canonical_gauge(F, R, Nijk)`, and tuple-keyed `GaugeTransform` remain
  supported as multiplicity-free wrappers.
- Multiplicityful vector-backed gauge transforms remain intentionally
  unsupported; explicit Hom-basis-aware metadata is accepted at the FRData
  boundary for future TensorCategories integration.

## v0.8.2 - Braid Representations and Finite-Field Diagnostics

### Added
- Fusion tree bases for multiplicity-free fusion spaces.
- Braid group representation construction from F/R symbols.
- Braid relation checks.
- Finite-field reduction of braid representations.
- Finite group diagnostics over F_p.
- Generated matrix algebra diagnostics.
- Commutant and irreducibility diagnostics.
- Experimental Zariski closure diagnostics.
- Examples for Semion, Fibonacci, and Ising braid representations.

### Limitations
- Multiplicity > 1 fusion rules are not yet supported.
- Zariski closure computation is diagnostic only, not a complete algebraic group computation.
- Finite-field closure data should be interpreted as arithmetic evidence, not as a proof of characteristic-zero Zariski density.
- Finite-field extension support is limited.
- Large generated groups are not exhaustively enumerated by default.

## v0.8.1 - F/R Equation Infrastructure

### Added

- Multiplicity-free F-symbol variable generation.
- Pentagon equation generation.
- Multiplicity-free R-symbol variable generation.
- Left/right hexagon equation generation.
- `FREquationSystem` abstraction.
- Gauge variable and gauge transformation infrastructure.
- Safe gauge fixing for small multiplicity-free examples.
- Finite-field reduction interface for FR equation systems.
- Experimental cyclotomic reconstruction interface.
- Examples for Semion, Toric code, Fibonacci, and Ising.

### Limitations

- Multiplicity > 1 fusion rules are not yet supported.
- Gauge fixing is safe but not canonical.
- Finite-field solving is limited.
- Cyclotomic reconstruction is experimental.
- This release provides infrastructure, not a complete MTC classifier.
