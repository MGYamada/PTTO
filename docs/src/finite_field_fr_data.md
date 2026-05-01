# Finite-Field FRData

## Computational path

ACMG treats finite-field F/R data as a solved fiber of the F/R equation
scheme, not as a reduction of hand-written complex matrices:

```text
fusion rules
  -> pentagon/hexagon scheme
  -> equation-level gauge fixing
  -> finite-field reduction
  -> F_p solution
  -> expansion to TensorCategories F/R order
  -> FRData{FpElem}
  -> finite-field braid representation
```

The high-level entry point is:

```julia
using ACMG

fib = solve_fr_mod_p(fibonacci_fusion_rules(), 101)
@assert verify_FRData(fib)

σ1, σ2 = braid_generators_B3(fib, [2, 2, 2]; total_charge = 2)
@assert σ1 * σ2 * σ1 == σ2 * σ1 * σ2
```

## Convenience constructors

Small examples are available through finite-field constructors:

```julia
semion = semion_fr_data_mod_p(17)
fib = fibonacci_fr_data_mod_p(101)
```

These call the finite-field solving path.  The legacy no-argument
`semion_fr_data`, `fibonacci_fr_data`, and `ising_fr_data` constructors have
been removed; new examples should use finite-field constructors or
`solve_fr_mod_p` directly.  `ising_fr_data_mod_p(p)` is available for
admissible primes such as `17`; the default test suite uses a precomputed
finite-field Ising fixture for braid and verification regressions so routine
tests do not pay the full solver cost.

`solve_fr_mod_p` is the finite-field F/R solver.  The older uppercase
`solve_FR_mod_p` name is deprecated and no longer exported.

## Verification

`verify_pentagon`, `verify_hexagon`, and `verify_FRData` reevaluate the
Phase-4 equations at the finite-field coordinates stored in `FRData{FpElem}`.
Braid examples should verify F/R data before constructing generators.

The solver records metadata such as `:p`, `:base_field`, `:solver_status`,
`:branch`, `:solution_count`, gauge-fixing metadata, and verification results
in `fr_metadata(fr)`.
