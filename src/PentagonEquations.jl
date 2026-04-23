"""
    PentagonEquations

Pentagon equation generation via `TensorCategories.pentagon_equations`.

The TensorCategories package builds the pentagon equations from a
fusion rule N_ij^k (as an `Array{Int,3}`) and a unit vector, using the
categorical convention

    f = (id_X ⊗ α_{Y,Z,W}) ∘ α_{X, Y⊗Z, W} ∘ (α_{X,Y,Z} ⊗ id_W)
    g = α_{X,Y,Z⊗W} ∘ α_{X⊗Y, Z, W}
    pentagon: f = g

The result is a list of polynomials in Oscar's `QQMPolyRingElem` whose
variables correspond to the independent entries of every associator block
F^{ijk}_o. The **variable-to-block correspondence** is defined by
TensorCategories' internal traversal order (nested loops over i, j, k, o
with a `pop!` from a shared stack). `assign_F_to_associator!` in
`HexagonEquations.jl` mirrors exactly this order so that numerical
solutions from the pentagon can be injected back into a category object.

Depends on: TensorCategories.jl (>= 0.1), Oscar.jl
"""

using Oscar
using TensorCategories


"""
    get_pentagon_system(Nijk::Array{Int,3}, r::Int) -> (R, eqs, n)

Return the Oscar polynomial ring `R`, the list of non-trivial pentagon
equations `eqs`, and the number of F-symbol variables `n`.

Arguments
- `Nijk`: fusion coefficients, `Nijk[i, j, k] = N_{ij}^k ∈ Z_{≥0}`.
          Index 1 is assumed to be the unit object.
- `r`:    rank (= `size(Nijk, 1)`). Passed explicitly for consistency with
          other Phase 4 routines.

Output
- `R`:    an Oscar `QQMPolyRing` in `n` variables x_1, …, x_n, each
          representing one entry of an associator block F^{ijk}_o.
- `eqs`:  `Vector{QQMPolyRingElem}` containing the pentagon relations
          after filtering out trivially-zero ones (unit axioms handled
          by TensorCategories produce some zero polynomials).
- `n`:    `nvars(R)`, the number of pentagon variables.

Notes
- TensorCategories' pentagon_equations uses 1-indexed unit at position 1.
  We pass `one_vec = [1, 0, 0, ..., 0]` to mark this.
- The **variable ordering is an internal implementation detail** of
  TensorCategories. Any code that consumes a solution from
  `solve_pentagon_*` must route the values through
  `assign_F_to_associator!` (in HexagonEquations.jl), which respects the
  same ordering.
"""
function get_pentagon_system(Nijk::Array{Int,3}, r::Int)
    size(Nijk) == (r, r, r) || error("Nijk must have shape ($r, $r, $r), got $(size(Nijk))")

    one_vec = zeros(Int, r)
    one_vec[1] = 1

    C, eqs_raw = pentagon_equations(Nijk, one_vec)
    eqs = filter(eq -> !(eq isa Integer) && !iszero(eq), eqs_raw)

    isempty(eqs) && error("Pentagon produced no non-trivial equations " *
                          "(rank=$r Nijk may be pointed/trivial)")

    R = parent(eqs[1])
    n = nvars(R)
    return R, eqs, n
end

