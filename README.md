# ACMG — Arithmetic Condensed Matter Geometry

ACMG.jl is a conductor-first experimental library for modular tensor
category arithmetic.  The base field is chosen up front as `Q(ζ_N)`,
and modular data, Galois action, Frobenius, and finite-field reduction
all keep that same conductor context.

```julia
using ACMG

ctx = CyclotomicContext(20)
fib = fibonacci_modular_data(ctx)

field(fib)        # Q(ζ_20)
conductor(fib)    # 20
cond_S(fib)       # 20
cond_T(fib)       # 5

galois_orbit(fib)
frobenius(fib, 41)
reduce_mod_p(fib, 41)
```

The built-in exact examples are:

```julia
semion = semion_modular_data(CyclotomicContext(8))
fib    = fibonacci_modular_data(CyclotomicContext(20))
ising  = ising_modular_data(CyclotomicContext(16))
```

## Design

The main design rule is:

```text
N -> K = Q(ζ_N) -> (F, R, S, T) over K -> Galois / F_p reduction
```

`N` is not discovered after constructing modular data.  It specifies
the cyclotomic ground field for the computation.

`CyclotomicContext` stores the conductor, Oscar field, chosen primitive
root, and related context:

```julia
struct CyclotomicContext
    N::Int
    field
    zeta
    real_subfield
    conductor::Int
end
```

`ModularData` always carries its context:

```julia
struct ModularData
    context::CyclotomicContext
    labels::Vector{Symbol}
    S
    T
    cond_S::Int
    cond_T::Int
    cond_F::Union{Int, Nothing}
end
```

## Pipeline

The conductor pipeline searches modular data exactly through finite
fields and CRT reconstruction:

```julia
mtcs = classify_mtcs_at_conductor(
    20;
    max_rank = 2,
    primes = [41, 61],
    skip_FR = false,
)
```

The current supported path is the arithmetic modular-data pipeline:

```text
N
  -> SL(2, Z/N) atomic irreps
  -> stratum enumeration
  -> block-U sweep over F_p
  -> CRT and Galois-aware grouping
  -> exact lift to Q(ζ_N)
  -> exact Phase 4 F/R roundtrip over Q(ζ_N)
```

Phase 4 computes exact cyclotomic `(F, R)` data for the built-in small
rank examples.  It caches finite-field Groebner preprocessing data,
enumerates complete split-prime point sets when they are small enough,
reconstructs F/R coordinates by bounded CRT in the power basis of
`Q(ζ_N)`, and verifies the lifted candidates exactly.  If the bounded CRT
lift is inconclusive, it falls back to exact triangular Groebner over
`Q(ζ_N)` and filters those exact candidates by their modular reductions.

### Sample conductor searches

The following small examples start only from the conductor `N`.  Modular
data are reconstructed directly in `Q(ζ_N)` by the same bounded
power-basis CRT layer used for exact `(F, R)` reconstruction.

```julia
using ACMG

samples = [
    (name = "Semion",    N = 8,  max_rank = 2, primes = [17, 41]),
    (name = "Ising",     N = 16, max_rank = 3, primes = [17, 97]),
    (name = "Fibonacci", N = 20, max_rank = 2, primes = [41, 61]),
]

for sample in samples
    classified = classify_mtcs_at_conductor(sample.N;
                                            max_rank = sample.max_rank,
                                            primes = sample.primes,
                                            verbose = false)

    @assert !isempty(classified)
    @assert any(m -> m.rank == sample.max_rank, classified)
    @assert all(m -> m.N == sample.N, classified)

    println(sample.name, ": ", length(classified), " classified MTC(s)")
end
```

For quick modular-data-only exploration, set `skip_FR = true`.  For a
more automatic search over rank and reconstruction-field candidates, use
`classify_mtcs_auto`:

```julia
semion_auto = classify_mtcs_auto(8;
                                 max_rank_candidates = [2],
                                 min_primes = 2,
                                 prime_start = 11,
                                 prime_max = 250,
                                 skip_FR = true,
                                 verbose = false)

@assert any(m -> m.rank == 2, semion_auto.classified)
```

## Galois And Frobenius

For `a` in `(Z/NZ)^*`, ACMG applies

```text
σ_a(ζ_N) = ζ_N^a
```

entrywise to context-carrying modular data:

```julia
galois_action(fib, 3)
```

For a prime `p` coprime to `N`, Frobenius is represented by the same
arithmetic action:

```julia
frobenius(fib, p)
```

Finite-field reduction uses the same context and the selected
primitive `N`-th root in `F_p`:

```julia
reduce_mod_p(fib, 41)
```

## Installation

```julia
julia> ]
pkg> activate .
pkg> instantiate
```

Dependencies:

- Oscar.jl for cyclotomic fields and GAP/SL2Reps access
- TensorCategories.jl for pentagon equation generation
- Primes.jl and LinearAlgebra

Julia 1.9 or newer is required.

## Repository Layout

```text
src/
  ACMG.jl                 single module root, include order, and exports
  Core/                   shared types and finite-field arithmetic
  Cyclotomics/            Q(ζ_N) context, Galois/Frobenius, F_p reduction
  ModularData/            exact S/T containers, validation, Verlinde helpers
  SL2/                    SL(2, Z/N) atomic representation catalog
  Search/                 stratum enumeration and block-U search
  Reconstruction/         CRT, Galois-aware grouping, cyclotomic lifting
  FR/                     exact pentagon/hexagon systems and F/R solvers
  Pipeline/               conductor-first orchestration and result records
  IO/                     reserved for serialization/reporting utilities
  Experimental/           reserved for incubating non-stable code
```

## Tests

Run the test suite with:

```julia
julia --project=. -e 'using Pkg; Pkg.test()'
```

The tests now emphasize static arithmetic unit checks and fast rank
2-3 pipeline smoke tests:

- semion at `N = 8`
- Fibonacci at `N = 20`
- Ising at `N = 16`

The Phase 4 roundtrip test checks semion, Fibonacci, and Ising over
their cyclotomic fields.

## References

- Bruillard–Ng–Rowell–Wang, *On classification of modular categories by rank*,
  arXiv:1507.05139.
- Ng–Rowell–Wang–Wen, *Reconstruction of modular data from SL_2(Z)
  representations*, arXiv:2203.14829.
- Ng–Rowell–Wen, *Classification of modular data up to rank 12*,
  arXiv:2308.09670.
