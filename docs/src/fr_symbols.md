# F/R Symbols

## What this page covers

The F/R layer describes associators and braidings for multiplicity-free fusion
rules.  It provides equation containers, finite-field solved `FRData`, F/R
accessors, and gauge-facing metadata.

## Minimal example

```julia
using ACMG

rules = fibonacci_fusion_rules()
system = fr_equation_system(rules)

@assert is_multiplicity_free(rules)
@assert validate_fr_system(system)
@assert system.metadata[:rank] == 2
```

## Mathematical meaning

For a multiplicity-free fusion rule, F-symbols and R-symbols become scalar
coordinates satisfying pentagon and hexagon equations.  ACMG records these
equations without committing the public boundary to one solver backend.

The current equation infrastructure targets multiplicity-free rules.  It is
the equation layer for pentagon/hexagon constraints, not a general
TensorCategories replacement.

## API overview

- Fusion rules: `semion_fusion_rules`, `fibonacci_fusion_rules`,
  `ising_fusion_rules`, `toric_code_fusion_rules`
- Equation layer: `FREquationSystem`, `fr_equation_system`,
  `pentagon_equations`, `hexagon_equations`, `validate_fr_system`
- Finite-field data: `solve_fr_mod_p`, `semion_fr_data_mod_p`,
  `fibonacci_fr_data_mod_p`, `ising_fr_data_mod_p`, `verify_FRData`
- Data container: `FRData`
- Accessors: `F_symbol`, `R_symbol`, `F_values`, `R_values`,
  `R_inverse_values`, `fr_metadata`

## Common pitfalls

- Multiplicity-free support does not imply support for arbitrary Hom-space
  multiplicities.
- Vector-backed F/R coordinates follow ACMG/TensorCategories ordering; prefer
  accessors unless you are working on internals.
- Small-example F/R data should be produced by the finite-field constructors
  or `solve_fr_mod_p`, not by hard-coded symbolic object labels.

## Stability notes

`FRData`, F/R accessors, fusion-rule constructors, and equation-system
construction are stable public API in v0.8.6.  Finite-field F/R solving is the
preferred small-example computation path, while cyclotomic reconstruction
remains experimental.
