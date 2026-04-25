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
rank examples using finite-field Groebner preprocessing and exact lift
back to `Q(ζ_N)`.

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
  ACMG.jl                 module root and exports
  Types.jl                core F_p arithmetic and Verlinde helpers
  SL2Reps.jl              SL(2, Z/N) atomic representation catalog
  StratumEnum.jl          stratum enumeration
  BlockU.jl               block-U search and F_p modular data
  CyclotomicContext.jl    exact Q(ζ_N) context and ModularData
  CRT.jl                  CRT and Galois-aware grouping
  PentagonEquations.jl    TensorCategories pentagon wrapper
  ModularDataLift.jl      exact lift from F_p / Z[sqrt(d)] to Q(ζ_N)
  Pipeline.jl             conductor-first search pipeline
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
