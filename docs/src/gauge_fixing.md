# Gauge Fixing

## What this page covers

The gauge layer records channel-scalar gauge parameters and conservative
normalization choices for multiplicity-free F/R data.

## Minimal example

```julia
using ACMG

data = semion_fr_data_mod_p(17)
plan = gauge_fixing_plan(data)

@assert validate_frdata_for_gauge(data)
@assert plan isa Vector
```

## Mathematical meaning

Changing bases in one-dimensional fusion channels rescales F-symbols and
R-symbols without changing the underlying braided category.  ACMG uses this as
solver preprocessing, representative selection, and unit-channel
normalization.

## API overview

- Records: `GaugeTransform`, `GaugeParameters`, `GaugeChoice`,
  `GaugeFixingResult`
- Transforms: `gauge_parameters`, `gauge_transform`,
  `gauge_fixing_plan`, `is_gauge_fixed`
- Equation-layer helpers: `gauge_variables`, `gauge_fix`

## Common pitfalls

- Scalar gauge transforms are currently for multiplicity-free F/R data.
- Gauge-fixed representatives are implementation choices, not new category
  invariants.
- Low-level finite-field gauge weights are useful for diagnostics but are not
  the stable user-facing interface.

## Stability notes

Typed gauge records, `gauge_parameters`, `gauge_fixing_plan`,
`is_gauge_fixed`, and `gauge_fix` are stable v0.8.6 API.  Low-level normal-form
internals, finite-field gauge action helpers, and stacky weights are
experimental/internal compatibility exports.
