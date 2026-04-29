# Higher central charges

ACMG computes the higher central charge of modular data by the exact
cyclotomic Gauss sum

```julia
xi_n = (1 / D) * sum(d[a]^2 * theta[a]^n for a in simples)
```

with `d[a] = S[1,a] / S[1,1]`, `D = 1 / S[1,1]`, and `theta[a] = T[a,a]`.
The built-in modular data use pure twists in `T`: the tensor unit has
`T[1,1] == 1`; no global `exp(-2*pi*i*c/24)` prefactor is included.

`xi_1` is the ordinary topological central charge phase
`exp(2*pi*i*c/8)`.

```julia
using ACMG

data = fibonacci_modular_data()
xi1 = higher_central_charge(data, 1; normalization = :D)
xis = higher_central_charges(data, 1:5; normalization = :D)
```

The historical exact API returns a `HigherCentralChargeResult`.  Its default
normalization is `:galois`, which computes `tau_n^+ / sigma_n(D)` and is
defined for `gcd(n, conductor) == 1`.  Use `normalization = :D` for the raw
higher central charge formula above at arbitrary `n`.

# Finite-field computation and cyclotomic reconstruction

The experimental finite-field path lives in `src/Experimental`.  It supports
the built-in sanity examples `:semion`, `:fibonacci`, `:ising`, and
`:toric_code`.

```julia
sol = solve_FR_mod_p(:fibonacci, 41)
xi = higher_central_charge(sol, 3)
```

The prototype reduces exact cyclotomic modular data to `F_p`, and for
Fibonacci can also store a reduced reference Phase-4 F/R solution.  Higher
central charges are evaluated in `F_p` from the gauge-invariant Gauss sum, so
tests compare final values rather than individual F- or R-symbol entries.

Cyclotomic reconstruction is intentionally minimal in this prototype:
`lift_higher_central_charge` is a placeholder, while exact comparisons are
done by reducing the cyclotomic answer with `reduce_mod_p`.

# Fibonacci worked example

```julia
using ACMG

data = fibonacci_modular_data()
sol = solve_FR_mod_p(:fibonacci, 41)

for n in 1:5
    exact = normalized_gauss_sum(data; n = n, normalization = :D)
    reduced_exact = reduce_mod_p(data.context, exact, 41)
    finite_field = higher_central_charge(sol, n).value
    @assert finite_field == reduced_exact
end
```

This checks the modular-data path and the finite-field F/R prototype on the
gauge-invariant quantity `xi_n`.
