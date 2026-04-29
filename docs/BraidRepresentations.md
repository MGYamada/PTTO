# Braid representations and finite-field diagnostics

ACMG v0.8.2 adds an experimental braid-representation layer for
multiplicity-free F/R data.

```julia
using ACMG

fr = fibonacci_fr_data()
br = braid_representation(fr, [:τ, :τ, :τ], :τ)

check_braid_relations(br).ok

brp = reduce_mod_p(br, 41; conductor = 20)
alg = generated_matrix_algebra(brp)
comm = commutant(brp)
zdiag = zariski_closure_diagnostics(brp; max_words = 200, max_degree = 2)
```

`FRData` uses the same TensorCategories-backed coordinate order as the FR
equation APIs:

- `F_values` follows `get_pentagon_system`.
- `R_values` and `R_inverse_values` follow `get_hexagon_fr_system`.

The built-in `semion_fr_data`, `fibonacci_fr_data`, and `ising_fr_data` are
hard-coded exact solutions in that order, and tests substitute them into
`pentagon_equations(rules)` and `hexagon_equations(rules)`.

The standard basis is the deterministic left-associated fusion-tree basis
`(((a1 a2 -> b2) a3 -> b3) ... an -> x)`.  The current implementation supports
only multiplicity-free fusion rules and scalar F/R symbols.  Higher
multiplicity needs matrix-valued local channels and is intentionally rejected
with an explicit error.

Finite-field reduction maps braid matrices to `F_p` when all entries can be
reduced and, for cyclotomic entries, when the requested conductor divides
`p - 1`.  Finite-field extensions are not implemented yet.

The Zariski closure API is diagnostic only.  It reports generated matrix
algebra dimension, commutant dimension, determinant and projective-order data,
and related small finite-field evidence.  It does not compute a full
characteristic-zero Zariski closure, does not classify algebraic groups, and
finite-field image data should be treated as arithmetic evidence rather than a
proof of density.
