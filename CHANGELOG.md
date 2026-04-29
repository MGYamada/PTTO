# Changelog

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
