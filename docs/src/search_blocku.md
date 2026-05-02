# Block-U Search

## What this page covers

The search layer enumerates SL(2, Z/N) strata and performs finite-field
Block-U searches used by the conductor-first pipeline.

For the mathematical definition of strata and the public classification
workflow, see [Classifying MTCs](@ref).  This page is a short orientation to
the lower-level search layer.

## Minimal example

```julia
using ACMG

result = classify_mtcs_auto(8;
    max_rank_candidates = [2],
    min_primes = 1,
    prime_start = 17,
    prime_max = 50,
    skip_FR = true,
    verbose = false,
)

@assert all(m -> m.N_input == 8, result.classified)
```

## Mathematical meaning

The conductor-level pipeline starts with `N`, builds arithmetic representation
data modulo `N`, searches finite-field candidates at admissible primes, groups
Galois-related data, and attempts exact cyclotomic reconstruction.

Block-U search is the fixed-stratum step: after choosing a semisimple
`SL(2, Z/N)` representation type, ACMG changes basis inside the allowed
`T`-eigenspace blocks and tests modular-data constraints over finite fields.

## API overview

- Stable pipeline: `classify_mtcs_at_conductor`, `classify_mtcs_auto`
- Strategy helpers: `select_admissible_primes`, `recommend_primes`,
  `recommend_skip_FR`
- Records: `ClassifiedMTC`, `FRRoundtripReport`, `fr_status`
- Lower-level search: `enumerate_strata`, `find_mtcs_at_prime`, Block-U
  enumeration helpers

## Common pitfalls

- Low-level enumeration order and pruning are implementation details.
- Heavy searches do not belong in docs examples.
- `skip_FR=true` is appropriate for fast modular-data-only exploration.

## Stability notes

The conductor-level classification entry points and result records are the
stable API.  Low-level Block-U helpers, fixed-stratum reconstruction hooks, and
internal search diagnostics remain experimental or semi-public compatibility
exports.
