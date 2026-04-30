# F/R Symbols

## What this page covers

The F/R layer describes associators and braidings for multiplicity-free fusion
rules.  It provides equation containers, exact built-in examples, F/R
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
- Exact data: `FRData`, `semion_fr_data`, `fibonacci_fr_data`,
  `ising_fr_data`
- Accessors: `F_symbol`, `R_symbol`, `F_values`, `R_values`,
  `R_inverse_values`, `fr_metadata`

## Common pitfalls

- Multiplicity-free support does not imply support for arbitrary Hom-space
  multiplicities.
- Vector-backed F/R coordinates follow ACMG/TensorCategories ordering; prefer
  accessors unless you are working on internals.
- Finite-field residues are not exact cyclotomic F/R data until reconstructed
  and verified.

## Stability notes

`FRData`, F/R accessors, fusion-rule constructors, and equation-system
construction are stable public API in v0.8.6.  General finite-field F/R solving
and cyclotomic reconstruction are experimental.
