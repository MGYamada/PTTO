# ACMG.jl

ACMG.jl is a research-oriented Julia library for conductor-first experiments
with exact modular data, F/R symbols, gauge choices, braid representations, and
finite-field diagnostics for modular tensor categories.

## Quick start

```julia
using ACMG

ctx = CyclotomicContext(20)
fib = fibonacci_modular_data(ctx)

@assert conductor(fib) == 20
@assert validate_exact_modular_data(fib).ok

hcc = higher_central_charge(fib; n = 1)
@assert hcc.ok
```

## Main design principle

ACMG.jl is conductor-first, and new F/R braid computations are equation-first:

```text
N -> Q(zeta_N) -> modular data/F/R -> Galois and F_p reduction

fusion rules -> pentagon/hexagon scheme -> F_p FRData -> braid matrices over F_p
```

The conductor `N` is not only an invariant inferred after the fact.  It selects
the cyclotomic ground field and keeps exact arithmetic, Frobenius/Galois
actions, finite-field reduction, and reconstruction in the same context.

## Stable in v0.8.6

- `CyclotomicContext` and conductor metadata helpers
- built-in exact modular data for semion, Fibonacci, Ising, and toric code
- exact modular-data validation, Verlinde checks, and Gauss sums
- higher central charge utilities over exact modular data
- multiplicity-free F/R equation containers and `FRData` accessors
- finite-field Phase-4 F/R solving via `solve_fr_mod_p`
- safe gauge records and gauge-fixing entry points
- exact multiplicity-free braid representation construction
- conductor-level pipeline entry points and classification result records
- JSON/Markdown export helpers for classification outputs

## Experimental features

The following remain experimental or semi-public compatibility exports:

- finite-field F/R solving and cyclotomic reconstruction hooks
- low-level Block-U search and fixed-stratum reconstruction helpers
- finite-field braid image diagnostics
- Zariski-style diagnostics
- low-level gauge normal-form and finite-field gauge internals

Diagnostic routines produce computational evidence, not classification
theorems.  Pin ACMG.jl when depending on experimental APIs.

## Examples

Runnable examples live in `examples/`.  The v0.8.6 examples are intentionally
small and use explicit conductors, small ranks, `skip_FR=true` where useful,
and sanity-check assertions.

```bash
julia --project=. examples/12_conductor_first.jl
julia --project=. examples/13_fr_gauge_braid.jl
julia --project=. examples/14_exports.jl
```

## Documentation

Build the local manual with:

```bash
julia --project=docs -e 'using Pkg; Pkg.instantiate(); include("docs/make.jl")'
```

Important pages:

- `docs/src/concepts.md`
- `docs/src/fr_symbols.md`
- `docs/src/finite_field_fr_data.md`
- `docs/src/api_stability.md`
- `docs/src/api.md`
- `docs/DocCoverage.md`

## References

ACMG.jl uses standard modular tensor category, cyclotomic field, and
finite-field methods.  The documentation describes the package conventions and
marks diagnostic routines carefully where they should not be read as proofs.

## Developer notes

Run tests before publishing documentation-focused changes:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
julia --project=docs -e 'using Pkg; Pkg.instantiate(); include("docs/make.jl")'
```

v0.8.6 is a documentation hardening release with no intentional functional
changes.
