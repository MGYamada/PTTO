# Modular Data

## What this page covers

The modular-data layer constructs and validates exact S/T data over a
`CyclotomicContext`.

## Minimal example

```julia
using ACMG

data = ising_modular_data(CyclotomicContext(16))

@assert conductor(data) == 16
@assert cond_S(data) == 16
@assert check_unitarity(data.S).valid
```

## Mathematical meaning

`ModularData` stores exact matrices in the cyclotomic field selected by its
context.  The conductor metadata distinguishes the conductor needed for `S`,
for `T`, and optionally for an attached F/R layer.

## API overview

- Constructors: `modular_data`, `semion_modular_data`,
  `fibonacci_modular_data`, `ising_modular_data`,
  `toric_code_modular_data`
- Metadata: `field`, `conductor`, `cond_S`, `cond_T`, `cond_F`
- Validation: `validate_modular_data`, `validate_exact_modular_data`,
  `check_modular_relations`, `check_unitarity`,
  `check_verlinde_integrality`
- Extraction: `verlinde_coefficient`, `extract_fusion_rule_Fp`,
  `extract_and_lift`

## Common pitfalls

- Built-in examples require their intended conductor.
- Finite-field extraction helpers are not a replacement for exact validation.
- `cond_T(data)` can be smaller than the ambient cyclotomic field conductor.

## Stability notes

Exact built-in constructors, validation functions, and conductor metadata are
stable public API.  Low-level finite-field extraction and reconstruction hooks
are semi-public or experimental depending on the input type.
