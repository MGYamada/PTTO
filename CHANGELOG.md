# Changelog

## Unreleased

### Changed
- Public API cleanup: internal `_`-prefixed helpers
  `_number_of_variables_in_hexagon_equations` and `_coerce_complex` are no
  longer exported from `ACMG`.
- Added public wrapper names
  `number_of_variables_in_hexagon_equations` and `coerce_complex` for users
  who need equivalent external access.

### Added
- Phase 2 (`BlockU.jl`) now supports a selectable block search backend via
  `search_mode` with a solver-oriented `:groebner` path and exhaustive
  `:exhaustive` path.
- Added Phase-2 solver controls:
  `max_units_for_groebner`, `groebner_allow_fallback`,
  and `precheck_unit_axiom`, threaded through
  `find_mtcs_at_prime`, `classify_mtcs_at_conductor`, and
  `classify_mtcs_auto`.
- Added Phase-2 algebraic utilities and helpers for equation building and
  filtering (`build_cayley_link_equations`, `build_verlinde_unit_equations`,
  `passes_unit_axiom`, `is_orthogonal_mod_p`), along with deterministic
  candidate ordering for solver extraction outputs.

### Deprecation / API organization
- This is a compatibility-oriented API organization change to avoid exposing
  internal-only names as part of the default public surface.
