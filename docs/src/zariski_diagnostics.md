# Zariski Diagnostics

## What this page covers

Zariski diagnostics summarize finite-field evidence from braid representation
matrices: generated algebra dimensions, commutants, bounded image samples, and
projective-order data.

## Minimal example

```julia
using ACMG

data = semion_fr_data()
br = braid_representation(data, [:s, :s, :s], :s)
ok = check_braid_relations(br)

@assert ok.ok
```

## Mathematical meaning

These routines help inspect algebraic behavior suggested by finite-field
reductions.  They do not compute a full characteristic-zero Zariski closure and
should not be cited as classification proofs without separate exact
arguments.

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
The stable neighbor API is exact `braid_representation` plus
`check_braid_relations`.
