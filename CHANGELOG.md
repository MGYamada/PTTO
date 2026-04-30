# Changelog

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
