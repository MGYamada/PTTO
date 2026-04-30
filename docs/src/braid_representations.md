# Braid Representations

## What this page covers

ACMG constructs braid group representation matrices from exact F/R data on
multiplicity-free fusion spaces.

## Minimal example

```julia
using ACMG

data = fibonacci_fr_data()
br = braid_representation(data, [:tau, :tau, :tau], :tau)

@assert length(braid_generators(br)) == 2
@assert check_braid_relations(br).ok
```

## Mathematical meaning

A fusion-tree basis expresses the state space for a sequence of simple
objects.  F-moves and R-symbols produce the standard braid generators on that
basis.  ACMG keeps this exact when the underlying F/R data are exact.

## API overview

- Basis objects: `FusionPath`, `FusionTreeBasis`, `fusion_paths`,
  `fusion_basis`, `dim`
- Representations: `BraidRepresentation`, `braid_representation`
- Matrix access: `braid_generator`, `braid_generators`
- Checks: `check_braid_relations`

## Common pitfalls

- Exact braid construction is separate from finite-field diagnostics.
- A finite-field image sample is not a proof about a characteristic-zero image.
- Multiplicity support follows the current `FRData` accessor limitations.

## Stability notes

Exact multiplicity-free braid construction and relation checking are stable in
v0.8.6.  Finite-field reductions, generated subgroup diagnostics, matrix
algebra diagnostics, commutants, and Zariski-style reports are experimental.
