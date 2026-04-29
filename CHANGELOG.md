# Changelog

## v0.8.0 - F/R Equation Infrastructure

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
