# Getting Started

This page shows the normal public workflow for constructing exact modular
data and running basic checks.

Representative functions:

- `CyclotomicContext`
- `fibonacci_modular_data`, `ising_modular_data`, `semion_modular_data`
- `conductor`, `cond_S`, `cond_T`
- `galois_orbit`, `frobenius`, `reduce_mod_p`

```julia
using ACMG

ctx = CyclotomicContext(20)
fib = fibonacci_modular_data(ctx)

conductor(fib)
galois_orbit(fib)
reduce_mod_p(fib, 41)
```

Notes: the conductor is part of the data model.  Choose it before constructing
modular data, and prefer documented public functions for ordinary use.
