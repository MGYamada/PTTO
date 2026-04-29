# F/R equation infrastructure

ACMG v0.8.1 uses TensorCategories.jl as the single computation substrate for
pentagon and hexagon equations.  Public helpers return the same equations as
the internal `PentagonEquations.jl` / `HexagonEquations.jl` paths.

```julia
using ACMG

rules = fibonacci_fusion_rules()

pent = pentagon_equations(rules)
hex = hexagon_equations(rules)

system = fr_equation_system(rules)
validate_fr_system(system)

fixed = gauge_fix(system; strategy = :safe)
fp_system = reduce_mod_p(fixed, 11)
```

`pentagon_equations(rules)` is TensorCategories-backed and uses the same
variable ordering as `get_pentagon_system`.

`hexagon_equations(rules)` constructs a TensorCategories-backed polynomial
system in one ring containing both F-symbol variables and R-symbol variables.
It is the public wrapper for `get_hexagon_fr_system`.

```julia
H, hex, nvars = get_hexagon_fr_system(rules.N, rules.rank)
```

The older `hexagon_equations(rules, F_values; context = ctx)` specialization
has been retired as a public API.  The Phase-4 solver still uses
`get_hexagon_system(Nijk, rank, F_values; context = ctx)` internally after
solving pentagon equations, but public equation generation no longer requires
concrete F-values.

The old `fsymbol_variables()` and `rsymbol_variables()` helpers have been
removed.  The authoritative variable ordering now lives in the returned
TensorCategories/Oscar polynomial rings and in `get_pentagon_system` /
`get_hexagon_fr_system` metadata.

# Scope and limitations

Only multiplicity-free fusion rules are supported in this layer.  If a fusion
coefficient is greater than `1`, the API throws an explicit error; higher
multiplicity needs matrix-valued F/R symbols and is left for a later release.

`gauge_fix(system; strategy = :safe)` is not a canonical gauge.  It only adds
normalizations that are visibly forced by unit channels and records what was
fixed in metadata, along with residual gauge information.

Finite-field support is currently reduction metadata for FR equation systems.
`reduce_mod_p` checks that `p` is prime.  TensorCategories polynomial systems
are kept in their internal representation; a general finite-field FR solver is
not implemented.

Cyclotomic reconstruction is experimental.  The v0.8 API provides metadata and
validation hooks such as `frobenius_metadata` and `cyclotomic_reconstruct`, but
the full reconstruction algorithm is intentionally left as future work.
