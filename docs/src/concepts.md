# Concepts

## What this page covers

ACMG.jl is organized around a conductor-first pipeline:

```text
N -> Q(zeta_N) -> modular data/F/R -> Galois and F_p reduction
```

The conductor `N` is a computational input.  It chooses the cyclotomic ground
field and therefore fixes the arithmetic context used by exact modular data,
F/R data, Galois action, Frobenius, and finite-field reduction.

## Minimal example

```julia
using ACMG

ctx = CyclotomicContext(20)
fib = fibonacci_modular_data(ctx)

@assert conductor(ctx) == 20
@assert conductor(fib) == 20
@assert cond_T(fib) == 5
@assert validate_exact_modular_data(fib).ok
```

## Mathematical meaning

For a modular tensor category whose modular data lie in a cyclotomic field, the
choice of `Q(zeta_N)` controls exact arithmetic.  ACMG keeps this choice
explicit rather than reconstructing it implicitly from floating-point or
symbolic output.

Exact cyclotomic arithmetic and finite-field reduction are meant to share the
same conductor context.  A finite-field calculation is useful only after the
prime has been interpreted relative to the conductor, for example through
split-prime or Frobenius metadata.

## API overview

The most common entry points are:

- `CyclotomicContext`, `field`, `zeta`, `conductor`
- `ModularData`, `semion_modular_data`, `fibonacci_modular_data`,
  `ising_modular_data`, `toric_code_modular_data`
- `galois_action`, `galois_orbit`, `frobenius`, `reduce_mod_p`
- `fr_equation_system`, `FRData`, `gauge_fix`, `braid_representation`
- `classify_mtcs_at_conductor`, `classify_mtcs_auto`

## Common pitfalls

- Do not treat `N` as merely a value discovered after building `S` and `T`.
  In ACMG, `N` selects the exact ground field.
- Do not compare finite-field residues from unrelated conductors.
- Do not treat diagnostic routines as classification proofs without exact
  verification.
- Keep heavy searches out of doctests and examples; use `skip_FR=true`, small
  ranks, and explicit small prime lists for documentation snippets.

## Stability notes

Conductor contexts, exact built-in modular data, Galois/Frobenius utilities,
basic F/R containers, gauge accessors, exact braid representations, and
pipeline entry points are the stable public surface.

Low-level Block-U search helpers, finite-field F/R solving, cyclotomic
reconstruction hooks, finite-field braid image diagnostics, and Zariski
diagnostics are experimental implementation APIs.
