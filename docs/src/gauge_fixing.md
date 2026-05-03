# Gauge Fixing

## What this page covers

The gauge layer treats a change of fusion-channel bases as an algebraic group
action on `FRData`.  In v0.9.3, ACMG represents the full general gauge group
of a fusion tensor, while active canonical gauge fixing remains implemented
only for the multiplicity-free toric case.

## Mathematical meaning

F-symbols and R-symbols are coordinates, not invariants by themselves.  If the
basis vector of `Hom(a ⊗ b, c)` is rescaled by `u[a,b,c]`, then ACMG uses the
TensorCategories convention

```text
F'^{abc}_{d;e,f} = u[a,b,e] u[e,c,d] / (u[b,c,f] u[a,f,d]) F^{abc}_{d;e,f}
R'^{ab}_c        = u[a,b,c] / u[b,a,c] R^{ab}_c
```

The pentagon and hexagon equations are preserved by this action.  Gauge fixing
chooses a representative in this orbit by making explicit constraints, solving
for gauge parameters, and applying the resulting `GaugeAction`.

ACMG's current toric convention is the full channel-scalar gauge group.  Unit
channels such as `(1,a,a)`, `(a,1,a)`, and `(1,1,1)` are included among the
parameters.  This is useful for stacky automorphism counts, but it means the
action normally carries an ineffective kernel.  Metadata produced by finite
gauge-group helpers records this as `:gauge_convention =>
:full_channel_scalar`, `:includes_unit_channels => true`, and
`:includes_ineffective_kernel => true`.  Unit-normalized and effective toric
gauge groups are separate conventions, not the default used here.

## GeneralGauge

For general fusion multiplicities, gauge transformations are changes of basis
in the fusion spaces

```text
V_ab^c = Hom(a ⊗ b, c),      dim V_ab^c = N_ab^c.
```

The gauge group is therefore

```text
G = ∏_{a,b,c} GL(N_ab^c),
```

omitting zero-dimensional channels.  Its algebraic dimension is

```text
dim G = Σ_{a,b,c} (N_ab^c)^2.
```

The multiplicity-free case is the special case `N_ab^c ∈ {0,1}`, where every
positive factor is `GL(1) = G_m`; this is exactly the toric gauge group used
by the existing Smith-normal-form preconditioning.

The new foundation types are `FusionSpaceIndex`, `GaugeFactor`,
`GeneralGaugeData`, and `GaugeTransformation`.  Use
`general_gauge_data(fusion)` to construct the product of general linear
groups, `gauge_group_dimension(gauge)` to compute the formula above, and
`identity_gauge_transformation(gauge)` plus
`validate_gauge_transformation(gauge, transformation)` to validate basis
changes.  Validation checks that every positive-dimensional channel has a
matrix, that no zero channel appears, that matrix sizes match `N_ab^c`, and
that each matrix is invertible.

In v0.9.3, the general action API is intentionally conservative:
`apply_gauge_to_F`, `apply_gauge_to_R`, and `apply_gauge_to_FR` support the
identity transformation for all `GeneralGaugeData`, and the scalar toric
action for `GL(1)` factors.  Nontrivial nonabelian actions for factors
`GL(n)` with `n > 1` throw `GeneralGaugeActionNotImplementedError` rather
than silently applying an incorrect formula.  Stabilizer and quotient-facing
hooks such as `gauge_stabilizer` and `gauge_orbit_dimension` are placeholders
for future finite-field and GIT workflows; they do not fake stabilizer
results.

## API overview

- Records: `GaugeAction`, `GaugeTransform`, `GaugeParameters`,
  `GaugeChoice`, `GaugeFixingResult`
- Actions: `identity_gauge`, `apply_gauge`, `compose_gauge`,
  `inverse_gauge`, `validate_gauge_action`
- Constraints: `GaugeConstraint`, `FixUnitConstraints`,
  `FixSelectedFSymbols`, `FixSelectedRSymbols`, `NormalizationConstraint`
- Normal forms: `gauge_degrees_of_freedom`, `build_gauge_constraints`,
  `solve_gauge_constraints`, `gauge_normal_form`, `validate_gauge_fixed`
- Compatibility wrappers: `canonical_gauge`, `gauge_equivalent`,
  `gauge_transform`, `gauge_fixing_plan`, `is_gauge_fixed`
- Equation-layer helpers: `gauge_variables`, `gauge_fix`

## Relation to finite fields

`GaugeAction` is parameterized by the scalar values themselves.  No complex
phase, conjugation, square-root, or approximate comparison is required for the
basic action.  Over `FRData{FpElem}`, ACMG validates gauge moves by exact
finite-field pentagon and hexagon evaluation through `verify_pentagon`,
`verify_hexagon`, and `verify_FRData`.

The low-level toric helpers `gauge_weight_matrix`, `smith_gauge_split`,
`stabilizer_size_mod_p`, and `stacky_weight_mod_p` expose the character matrix
of the gauge torus.  They are useful for future quotient and zeta-function
work, but their detailed output is still experimental.

## Toric preconditioning in `classify_mtcs`

For multiplicity-free fusion rules, every nonzero trivalent channel
`Hom(a ⊗ b, c)` is one-dimensional.  The scalar gauge group is therefore a
torus with one parameter for each channel with `N_ab^c = 1`; over a finite
field these parameters are restricted to `F_p^×`.

In v0.9.2, `classify_mtcs_at_conductor` and `classify_mtcs_auto` compute this
toric data immediately before exact F/R reconstruction.  When the fusion rule
is multiplicity-free, ACMG chooses a deterministic Smith-normal-form F-symbol
slice and asks the reconstruction solver to set those selected coordinates to
`1`.  If any fusion coefficient is greater than one, the step is skipped by
default and the previous reconstruction path is used.

Use `gauge_fixing = :none` or `toric_gauge_fixing = false` to disable active
preconditioning.  The `gauge_fixing` keyword accepts:

```julia
:auto
:toric
:general
:none
```

`:auto` uses toric gauge fixing exactly when all `GeneralGauge` factors are
`GL(1)`.  If multiplicities greater than one occur, ACMG records
`GeneralGaugeData` and does not attempt nonabelian gauge fixing.  `:toric`
requires the toric condition and throws `ToricGaugeFixingError` otherwise.
`:general` constructs and validates the general gauge metadata but does not
claim a canonical representative.  `:none` disables active fixing.

```julia
using ACMG

result = classify_mtcs_at_conductor(20;
    max_rank = 2,
    primes = [41, 61],
    gauge_fixing = :auto,
    verbose = false,
)
```

The selected `(F,R)` representative may differ from older output by a gauge
change.  It should still satisfy the same pentagon, hexagon, and modular-data
roundtrip checks.

At solve time, `gauge_fix(fr_equation_system(rules))` constructs a Smith
normal form slice for the F-symbol torus action and substitutes those selected
F coordinates by `1` before the finite-field Groebner solve.  R-symbols are
not gauge-fixed in this pre-solver step; they remain part of the braiding
branch selected by the hexagon equations.

`toric_gauge_data(fr)` packages the solved F/R character matrix with its Smith
normal form split.  `toric_gauge_normal_form(fr)` reports the pre-solver
F-slice, stabilizer size, and stacky weight from metadata; it does not apply a
new post-hoc gauge transform.  Braid representations and finite-field higher
central charge evaluation consume the solved FRData directly.

## Normal forms and validation

The normal-form flow is deliberately separated into stages:

This keeps representative selection separate from validation.  At the
equation layer, `gauge_fix(system; strategy=:toric_snf)` uses an F-only toric
Smith-normal-form slice, substitutes those F-symbol coordinates by `1`, and
records the free-index map needed to expand reduced solver points back to
TensorCategories order.  Use `strategy=:none` to disable equation-level gauge
fixing.

## Migration notes

Existing compatibility code using `GaugeTransform`, `gauge_transform`,
`canonical_gauge`, `gauge_fixing_plan`, or `is_gauge_fixed` remains supported.
New code should prefer `GaugeAction` and `apply_gauge` for FRData-centered
workflows.  `GaugeTransform` is a legacy multiplicity-free wrapper, while
`GaugeAction` records explicit Hom-basis keys `(a,b,c,μ)` and backend metadata.

## Common pitfalls

- Scalar gauge actions currently require multiplicity-free vector-backed
  `FRData`.
- Gauge-fixed representatives are implementation choices, not category
  invariants.
- `gauge_equivalent` is a deterministic F-slice based comparison.  It is
  reliable when the selected F-slice is complete, or when the residual gauge
  action is known to be trivial on R-symbols.
- `target=:F` and `target=:R` in `apply_gauge` are diagnostic tools; a full
  braided category gauge move should normally use `target=:FR`.
- Avoid approximate tests for exact or finite-field data.  Use exact equality
  or the validation helpers whenever possible.

## Stability notes

`GaugeAction`, `identity_gauge`, `apply_gauge`, `compose_gauge`,
`inverse_gauge`, constraint records, and normal-form orchestration are public
infrastructure.  Low-level quotient diagnostics and future solver-specific
gauge-fixing strategies remain experimental.
