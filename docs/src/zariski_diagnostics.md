# Zariski Diagnostics

## What this page covers

Zariski diagnostics summarize finite-field evidence from braid representation
matrices: generated algebra dimensions, commutants, bounded image samples, and
projective-order data.

## Minimal example

```julia
using ACMG

data = semion_fr_data_mod_p(17)
br = braid_representation(data, [2, 2, 2], 2)
brp = reduce_mod_p(br)
ok = check_braid_relations(brp)

@assert ok.ok
@assert zariski_closure_diagnostics(brp).diagnostic_only
```

## Mathematical meaning

These routines inspect algebraic behavior suggested by finite-field braid
matrices.  They do not compute a full characteristic-zero Zariski closure and
should not be cited as classification proofs without separate exact arguments.

## API overview

- `generated_subgroup`
- `finite_group_diagnostics`
- `generated_matrix_algebra`
- `commutant`
- `zariski_closure_diagnostics`

## Common pitfalls

- Bounded enumeration may truncate.
- Finite-field irreducibility or image size is evidence, not a theorem about
  the exact representation.
- Diagnostic return schemas may change before v1.0.

## Stability notes

All Zariski and finite-field image diagnostics are experimental in v0.8.6.
The stable neighbor API is finite-field `FRData` construction plus exact
`check_braid_relations` over the chosen prime field.
