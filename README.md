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

Exact Gauss sums and higher central charges use the same cyclotomic
backend:

```julia
ising = ising_modular_data(CyclotomicContext(16))

gauss_sum_plus(ising)                         # raw τ_1^+
normalized_gauss_sum(ising; normalization=:D) # τ_1^+ / D = ζ_16
central_charge(ising)                         # ordinary central charge

higher_central_charge(ising; n=3)             # τ_3^+ / σ_3(D)
higher_central_charges(ising, [1, 3, 5])
```

`higher_central_charge` returns a `HigherCentralChargeResult` so cases
without a Galois normalization, such as non-units modulo the conductor,
can be reported structurally.  Core values remain exact elements of
`Q(ζ_N)`; numerical evaluation is only a display concern.

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

v0.8 also exposes a lower-level, backend-neutral equation layer for
multiplicity-free fusion rules:

```julia
rules = fibonacci_fusion_rules()
system = fr_equation_system(rules)
fixed = gauge_fix(system; strategy = :safe)
fp_system = reduce_mod_p(fixed, 11)
```

This layer generates F-symbol variables, R-symbol variables, pentagon
equations, left/right hexagon equations, safe unit-channel gauge fixing,
and finite-field reductions without requiring a solver backend at the
representation boundary.  See `docs/FRInfrastructure.md`.

`ClassifiedMTC` records carry an explicit `fr_status` field.  Use
`fr_status(m)` to distinguish skipped, solved, reconstruction-failed, and
verification-failed Phase 4 outcomes.

### Exporting results

Classification outputs can be archived as JSON and summarized as Markdown:

```julia
classified = classify_mtcs_at_conductor(8; max_rank = 2, verbose = false)

save_classification("N8.json", classified)
loaded = load_classification("N8.json")

export_modular_data(classified[1])
export_fusion_rule(classified[1])
export_FR(classified[1])

write_report("N8.md", classified)
```

JSON export stores exact cyclotomic values in stable textual form.  The
loader returns the saved payload as a `Dict`.

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
galois_action(fib.S, fib.T, 3; context = fib.context)
```

Fusion-rule and modular-data automorphisms are exposed as unit-fixing
1-indexed permutations:

```julia
Nijk = check_verlinde_integrality(semion.S).Nijk
fusion_automorphisms(Nijk)
modular_data_automorphisms(semion)
```

The Galois action on anyon labels is recovered from projective rows of
`S`:

```julia
galois_anyon_action(ising, 3)
galois_anyon_orbits(ising)
```

For conductor-wise smoke checks, use the compact sanity table:

```julia
conductor_sanity_table([semion, ising])
println(conductor_sanity_markdown([semion, ising]))
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
  Gauge/                  gauge representatives, transforms, and fixing helpers
  Pipeline/               conductor-first orchestration and result records
  IO/                     JSON serialization and Markdown reporting
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
