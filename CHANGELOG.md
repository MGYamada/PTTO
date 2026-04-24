# Changelog

## Unreleased

### Changed
- Public API cleanup: internal `_`-prefixed helpers
  `_number_of_variables_in_hexagon_equations` and `_coerce_complex` are no
  longer exported from `ACMG`.
- Added public wrapper names
  `number_of_variables_in_hexagon_equations` and `coerce_complex` for users
  who need equivalent external access.

### Deprecation / API organization
- This is a compatibility-oriented API organization change to avoid exposing
  internal-only names as part of the default public surface.
