# Higher Central Charge

ACMG treats higher central charges as normalized moments of the twist spectrum.
For modular data with quantum dimensions `d_a`, twists `θ_a`, and

```math
D^2 = \sum_a d_a^2,
```

the `n`-th higher central charge is

```math
\xi_n(C) = {1 \over D^2}\sum_a d_a^2 \theta_a^n.
```

This `D^2` normalization is the only HCC normalization used by the public
moment API.  It makes the coefficients

```math
\mu(a) = d_a^2 / D^2
```

a probability distribution on the simple objects, so `ξ_n` is the `n`-th
twist-spectrum moment and `ξ_0 = 1`.

## Exact API

```julia
using ACMG

data = fibonacci_modular_data()

higher_central_charge(data, 0)
higher_central_charge(data, 1)
higher_central_charges(data, 0:4)
higher_central_charge_period(data)
higher_central_charge_sequence(data)
```

`higher_central_charge_sequence(data)` returns one full period.  The period is
currently the twist conductor recorded on the modular data.

The structured generating-function helper stores the moment weights and
twists:

```julia
gf = higher_central_charge_generating_function(data)
gf(3) == higher_central_charge(data, 3)
```

Mathematically this represents

```math
\Xi_C(u) = \sum_{n \ge 0} \xi_n u^n
         = {1 \over D^2}\sum_a {d_a^2 \over 1 - \theta_a u}.
```

## Examples

For Semion,

```math
d = (1, 1), \qquad \theta = (1, i), \qquad D^2 = 2,
```

so

```math
\xi_n = {1 + i^n \over 2}
```

and the period is `4`.

For Fibonacci,

```math
d_1 = 1, \qquad d_\tau = \varphi, \qquad
\theta_1 = 1, \qquad \theta_\tau = \zeta_5^2,
```

so

```math
\xi_n = {1 + \varphi^2 \zeta_5^{2n} \over 1 + \varphi^2}
```

and the period is `5`.

For Ising,

```math
d = (1, \sqrt{2}, 1), \qquad
\theta = (1, \exp(\pi i/8), -1), \qquad D^2 = 4,
```

so

```math
\xi_n = {1 + 2\exp(n\pi i/8) + (-1)^n \over 4}
```

and the period is `16`.

## Finite Fields

The experimental finite-field API reduces the exact `D^2`-normalized value at
good primes:

```julia
higher_central_charge_modp(data, n, p; embedding = nothing)
higher_central_charge_sequence_modp(data, p; embedding = nothing)
```

For cyclotomic data, `embedding` can specify a prime-field root of unity as
`order => root`.  The Fibonacci regression at `p = 11` uses

```julia
fib = fibonacci_modular_data()
higher_central_charge_sequence_modp(fib, 11; embedding = 5 => 3)
# [1, 6, 7, 5, 9]
```

This corresponds to `ζ_5 ↦ 3 ∈ F_11`.  Then

```math
\varphi \mapsto 1 + 3 + 3^{-1} = 8, \qquad
\varphi^2 \mapsto 9, \qquad D^2 \mapsto 10, \qquad
\theta_\tau \mapsto 3^2 = 9,
```

and

```math
\xi_{n,11} = (1 + 9\cdot 9^n) / 10.
```

Bad primes are rejected explicitly.  In particular, primes dividing the twist
conductor, or primes where `D^2` is not invertible after reduction, are not
accepted.  ACMG does not coerce approximate complex numbers into finite
fields.

Prime-field reductions are implemented now.  If the requested root of unity
does not lie in `F_p`, ACMG reports that an extension field or larger
embedding is needed.

## Local Factors

The local-factor helpers are the first building block for future HCC zeta
experiments:

```julia
hcc_local_factor(data, n, p; embedding = nothing)
hcc_local_factors(data, p; embedding = nothing)
```

When `ξ_{n,p} ∈ F_p`, the local factor is

```math
P_{p,n}(T) = 1 - \xi_{n,p}T.
```

For the Fibonacci `p = 11` example, the one-period local factors are

```math
1 - T,\quad
1 - 6T,\quad
1 - 7T,\quad
1 - 5T,\quad
1 - 9T.
```

Full global HCC zeta functions are future work.  The current API provides the
exact moments, finite-field reductions, and local arithmetic factors needed for
those experiments.
