# F/R equation infrastructure

ACMG v0.8 adds a backend-neutral layer for generating polynomial equations
from multiplicity-free fusion rules.  It is intended as infrastructure for
passing F-symbol, R-symbol, pentagon, hexagon, gauge, and finite-field data to
future algebra backends.

```julia
using ACMG

rules = fibonacci_fusion_rules()

fvars = fsymbol_variables(rules)
rvars = rsymbol_variables(rules)

pent = pentagon_equations(rules, fvars)
hex = hexagon_equations(rules, fvars, rvars)

system = fr_equation_system(rules)
validate_fr_system(system)

fixed = gauge_fix(system; strategy = :safe)
fp_system = reduce_mod_p(fixed, 11)
```

The equation representation stores sparse polynomial expressions using
`EquationVariable`, `EquationTerm`, `EquationExpr`, and `PolynomialEquation`.
It does not require Symbolics.jl, Nemo.jl, Oscar.jl, or TensorCategories.jl at
the representation boundary, although older exact Phase-4 APIs in ACMG still
use Oscar/TensorCategories internally.

# Scope and limitations

Only multiplicity-free fusion rules are supported in this layer.  If a fusion
coefficient is greater than `1`, the API throws an explicit error; higher
multiplicity needs matrix-valued F/R symbols and is left for a later release.

`gauge_fix(system; strategy = :safe)` is not a canonical gauge.  It only adds
normalizations that are visibly forced by unit channels and records what was
fixed in metadata, along with residual gauge information.

Finite-field support is currently reduction infrastructure.  `reduce_mod_p`
checks that `p` is prime and reduces integer/rational coefficients in the
equation system to `F_p`; a general finite-field FR solver is not implemented.

Cyclotomic reconstruction is experimental.  The v0.8 API provides metadata and
validation hooks such as `frobenius_metadata` and `cyclotomic_reconstruct`, but
the full reconstruction algorithm is intentionally left as future work.
