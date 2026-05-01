# Braid Representations

## What this page covers

ACMG constructs braid group representation matrices from `FRData` on
multiplicity-free fusion spaces, including finite-field solved `FRData`.

## Minimal example

```julia
using ACMG

data = fibonacci_fr_data_mod_p(101)
br = braid_representation(data, [2, 2, 2], 2)

@assert length(braid_generators(br)) == 2
@assert check_braid_relations(br).ok
```

## Mathematical meaning

A fusion-tree basis expresses the state space for a sequence of simple
objects.  F-moves and R-symbols produce the standard braid generators on that
basis.  When the underlying data are `FRData{FpElem}`, the matrices are exact
over the prime field.

## API overview

- Basis objects: `FusionPath`, `FusionTreeBasis`, `fusion_paths`,
  `fusion_basis`, `dim`
- Representations: `BraidRepresentation`, `braid_representation`
- Matrix access: `braid_generator`, `braid_generators`,
  `braid_generators_B3`
- Checks: `check_braid_relations`

## Common pitfalls

- New finite-field braid computations should consume `solve_fr_mod_p` output,
  not legacy hard-coded small-example F/R data.
- A finite-field image sample is not a proof about a characteristic-zero image.
- Multiplicity support follows the current `FRData` accessor limitations.

## Stability notes

Exact multiplicity-free braid construction and relation checking are stable in
v0.8.6.  Finite-field reductions, generated subgroup diagnostics, matrix
algebra diagnostics, commutants, and Zariski-style reports are experimental.
