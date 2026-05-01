# Finite Fields

## What this page covers

Finite-field support reduces conductor-aware exact data modulo primes and
provides lightweight diagnostic hooks.

## Minimal example

```julia
using ACMG

data = semion_modular_data()
fp = reduce_mod_p(data, 17)

@assert fp.p == 17
@assert fp.rank == 2
```

## Mathematical meaning

Finite-field calculations are interpreted relative to the same conductor used
for exact cyclotomic arithmetic.  Frobenius action, split-prime behavior, and
residue comparisons all depend on that conductor.

## API overview

- Reduction: `reduce_mod_p`
- Frobenius metadata: `frobenius`, `frobenius_metadata`
- Root utilities: `find_zeta_in_Fp`, `cyclotomic_to_Fp`, `root_of_unity`
- Experimental F/R hooks: `solve_fr_mod_p`, `solve_finite_field`,
  `cyclotomic_reconstruct`

## Common pitfalls

- A residue modulo `p` is not an exact cyclotomic value.
- Use primes compatible with the conductor when comparing reductions.
- Keep docs examples small; finite-field searches can grow quickly.

## Stability notes

Exact modular-data reduction is part of the main workflow.  General
finite-field F/R solving, finite-field higher-central-charge lifting, and
cyclotomic reconstruction hooks are experimental.
