# Classifying MTCs

ACMG's classification pipeline is conductor-first.  The input is a positive
integer `N`, interpreted as the cyclotomic conductor through which the modular
representation of a modular tensor category is expected to factor:

```math
\rho : SL(2,\mathbb Z) \longrightarrow GL_r(K),
\qquad
K = \mathbb Q(\zeta_N),
```

with `ρ` factoring through the finite group

```math
G_N = SL(2,\mathbb Z/N\mathbb Z).
```

The pipeline searches for exact modular data `(S,T)` over `Q(ζ_N)`, validates
the Verlinde fusion rules, and optionally reconstructs exact `(F,R)` data.
This page describes the mathematical organization of the search.

## Strata

Let

```math
\operatorname{Irr}(G_N) = \{\rho_\lambda\}_{\lambda \in \Lambda_N}
```

be the finite catalog of irreducible `G_N`-representations used by ACMG,
realized over `Q(ζ_N)`.  A rank-`r` stratum is a finite-support multiplicity
function

```math
m : \Lambda_N \to \mathbb Z_{\ge 0},
\qquad
\sum_{\lambda \in \Lambda_N} m_\lambda \dim(\rho_\lambda) = r.
```

It determines the semisimple representation type

```math
\rho_m \cong
\bigoplus_{\lambda \in \Lambda_N}
\rho_\lambda^{\oplus m_\lambda}.
```

Equivalently, inside the representation variety of `G_N` in dimension `r`, a
stratum is the locally closed locus with this fixed irreducible decomposition.
The stratum is only a representation-theoretic layer.  It is not yet a modular
tensor category and not yet modular data in the simple-object basis.

This distinction matters for the tensor unit.  The unit object is a
distinguished basis vector after a candidate `(S,T)` has been assembled.  It
need not appear as a separate trivial one-dimensional summand of `ρ_m`; it may
lie inside a higher-dimensional irreducible summand whose `T`-spectrum contains
`1`.  For this reason, ACMG does not require a trivial summand during stratum
enumeration by default.

In Julia:

```julia
using ACMG

catalog = ACMG.build_atomic_catalog(8; max_rank = 3, verbose = false)
strata = ACMG.enumerate_strata(catalog, 2)
map(s -> ACMG.describe_stratum(s, catalog), strata)
```

`enumerate_strata` is a low-level combinatorial tool.  Most users should start
from `classify_mtcs_auto` or `classify_mtcs_at_conductor`.

## Block-U Search

For a fixed stratum, ACMG first builds a block-diagonal representation
realizing `ρ_m`.  The simple-object basis for a modular datum is not known a
priori.  The Block-U step searches for allowable changes of basis, especially
inside degenerate `T`-eigenspaces, and tests whether the transformed matrices
can satisfy the modular-data constraints.

The tests include:

- the modular relations for the chosen `S` and diagonal `T`,
- a distinguished unit object,
- Verlinde integrality and nonnegativity of the fusion coefficients,
- compatibility checks used later for exact cyclotomic reconstruction.

The finite-field search is arithmetic rather than numerical.  ACMG chooses
primes `p` with

```math
N \mid (p - 1),
```

so that the `N`-th roots of unity split in `F_p`.  Candidate matrices are
searched over such prime fields and grouped by their fusion tensors and Galois
data.

## Exact Reconstruction

Finite-field candidates are not final mathematical output.  They are residues
that guide exact reconstruction.  ACMG uses several admissible primes:

- a subset of primes is used to reconstruct cyclotomic entries in `Q(ζ_N)`,
- remaining fresh primes are used as cross-checks,
- candidates failing exact validation are discarded.

The result is exact modular data over `Q(ζ_N)`, not an approximate complex
matrix.

## F/R Reconstruction

After exact modular data are found, ACMG can attempt to solve the pentagon and
hexagon equations for `(F,R)`.  This is substantially more expensive than the
modular-data search.  Use `skip_FR = true` when you want a fast
modular-data-only exploration, and inspect `fr_status(m)` on each result when
F/R reconstruction is enabled.

Immediately before this F/R reconstruction step, ACMG v0.9.2 performs toric
gauge preconditioning for multiplicity-free fusion rules.  The pipeline checks
the fusion tensor only: all multiplicities must be `0` or `1`.  In that case,
it records the torus of trivalent-channel scalars and fixes a deterministic
Smith-normal-form slice of F-symbol coordinates to `1` before solving.  For
fusion rules with multiplicity greater than one, the step is skipped by
default.

Controls:

- `gauge_fixing = :auto`: apply toric preconditioning only when
  multiplicity-free.
- `gauge_fixing = :toric`: require the multiplicity-free toric step and error
  otherwise.
- `gauge_fixing = :none` or `toric_gauge_fixing = false`: preserve the old
  ungauge-fixed reconstruction path.

Gauge-fixed F/R representatives can differ from older output, but accepted
results must still pass the exact pentagon, hexagon, and modular-data
roundtrip checks.

A classification result is represented by `ClassifiedMTC`.  It stores:

- conductor and rank,
- the representation stratum,
- integer fusion tensor,
- primes used for reconstruction and fresh verification,
- exact cyclotomic `S` and `T`,
- optional exact `F_values` and `R_values`,
- an exact roundtrip verification report when available.

## Public Workflow

The convenience entry point is `classify_mtcs_auto`:

```julia
using ACMG

result = classify_mtcs_auto(8;
    max_rank_candidates = [2],
    min_primes = 2,
    prime_start = 17,
    prime_max = 100,
    skip_FR = true,
    verbose = false,
)

result.classified
```

For reproducible runs, choose the primes and rank bound explicitly:

```julia
using ACMG

classified = classify_mtcs_at_conductor(8;
    max_rank = 2,
    primes = [17, 41],
    skip_FR = true,
    verbose = false,
)
```

`classify_mtcs_auto` returns a named tuple with `classified`, selected primes,
rank cutoff, attempt count, and a history of stages.  `classify_mtcs_at_conductor`
returns the vector of `ClassifiedMTC` records directly.

## Interpreting Results

The modular-data layer of a `ClassifiedMTC` is exact:

```julia
using ACMG

m = first(classified)
m.N
m.rank
m.stratum
m.Nijk
m.S_cyclotomic
m.T_cyclotomic
fr_status(m)
```

When `fr_status(m) == ACMG.FRSolved`, ACMG also found exact `(F,R)` data and
verified that the reconstructed modular data match the target.  When
`skip_FR = true`, the result is intentionally modular-data-only.

## Scope

The conductor-level classification API is stable at the workflow boundary:

- `classify_mtcs_auto`
- `classify_mtcs_at_conductor`
- `ClassifiedMTC`
- `fr_status`
- `recommend_primes`
- `recommend_skip_FR`

Lower-level Block-U enumeration, Gröbner extraction hooks, reconstruction
internals, and search diagnostics are implementation details or experimental
interfaces.  They are useful for development and research experiments, but the
mathematical object to cite from the public pipeline is the exact validated
`ClassifiedMTC` output together with its recorded provenance.
