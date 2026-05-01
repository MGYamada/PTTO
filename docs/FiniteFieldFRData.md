# FiniteFieldFRData

The maintained manual page is `docs/src/finite_field_fr_data.md`.

Finite-field F/R data in ACMG follows:

```text
fusion rules -> pentagon/hexagon scheme -> gauge fixing -> F_p solution
  -> FRData{FpElem} -> finite-field braid representation
```

Use `solve_fr_mod_p(rules, p)` or the convenience constructors
`semion_fr_data_mod_p`, `fibonacci_fr_data_mod_p`, and `ising_fr_data_mod_p`.
