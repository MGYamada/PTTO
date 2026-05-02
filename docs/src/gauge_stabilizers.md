# Gauge Stabilizers

ACMG.jl v0.9.1 adds experimental infrastructure for computing gauge
stabilizers of already available solution data.  For a gauge action of a
group `G` on a solution space `X`, the automorphism group of a point is its
stabilizer:

```math
\operatorname{Aut}(x) = \operatorname{Stab}_G(x)
  = \{ g \in G \mid g \cdot x = x \}.
```

This page describes stabilizer computation infrastructure only.  ACMG.jl does
not implement quotient-stack classification, zeta functions, Lefschetz trace
sums, higher-central-charge trace sums, or étale cohomology.

## Why Stabilizers Matter

Stabilizers are the denominator needed by future stacky weights such as

```math
\frac{1}{|\operatorname{Aut}(x)|}.
```

They are also the local automorphism data that would be needed by future
Lefschetz-style weighted expressions of the form

```math
\sum_{[x]} \frac{\operatorname{trace}(x)}{|\operatorname{Aut}(x)|}.
```

The v0.9.1 API stops at computing or describing `Aut(x)`.  It intentionally
does not compute these sums.

## API

```julia
problem = StabilizerProblem(solution, gauge_group)
result = stabilizer(problem)

automorphisms(result)
stabilizer_order(result)
is_trivial_stabilizer(result)
stabilizer_metadata(result)
```

For a small explicitly enumerable finite group, pass the group elements and an
exact action:

```julia
using ACMG

p = 5
G = [ACMG.FpElem(i, p) for i in 1:(p - 1)]
action = (x, g) -> g * x

nonzero = stabilizer(ACMG.FpElem(2, p), G; action = action)
zero = stabilizer(ACMG.FpElem(0, p), G; action = action)

stabilizer_order(nonzero) # 1
stabilizer_order(zero)    # p - 1
```

The implementation uses exact equality and finite-field arithmetic.  It does
not convert exact or finite-field values to floating point.

## F/R Data

For multiplicity-free finite-field `FRData`, ACMG can construct the finite
toric gauge group internally:

```julia
fr = semion_fr_data_mod_p(17)
G = ACMG.finite_field_gauge_group(fr)
problem = StabilizerProblem(fr, G)

eqs = stabilizer_equations(problem)
```

`stabilizer_equations(problem)` returns character equations for the active
F/R symbol coordinates.  These equations are deterministic algebraic data for
future finite-field solving or counting code.

The finite-field group currently uses ACMG's full channel-scalar convention:
every channel with `N_ab^c = 1` contributes a parameter, including unit
channels.  Consequently the group includes the usual ineffective kernel.  The
group and equation metadata expose
`gauge_convention = :full_channel_scalar`,
`includes_unit_channels = true`, and
`includes_ineffective_kernel = true` so downstream counting code can
distinguish this from unit-normalized or effective toric quotients.

Brute-force enumeration is available for small finite gauge groups:

```julia
result = stabilizer(problem; max_enumeration = 100_000)
stabilizer_order(result)
```

If the finite gauge group is too large, `stabilizer` throws an explicit error
instead of returning a misleading partial result.  In that case,
`stabilizer_equations(problem)` is the intended fallback.

## Current Scope

Implemented in v0.9.1:

- `StabilizerProblem`, `StabilizerEquations`, and `StabilizerResult`
- brute-force stabilizer enumeration for small finite gauge groups
- exact structural comparison for `FRData`
- finite-field toric gauge-group support for multiplicity-free `FRData`
- deterministic stabilizer equation generation for active F/R coordinates

Not implemented in v0.9.1:

- zeta functions
- Euler factors
- Lefschetz trace sums
- higher-central-charge trace sums
- stacky point counting as a full feature
- automatic classification of quotient-stack points
