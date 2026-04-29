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
with a `pop!` from a shared stack). Exact solvers should preserve this
ordering when interpreting a point of the pentagon variety.

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
  TensorCategories. Exact solver code should preserve this ordering when
  converting between polynomial coordinates and associator blocks.
"""
function get_pentagon_system(Nijk::Array{Int,3}, r::Int)
    size(Nijk) == (r, r, r) || error("Nijk must have shape ($r, $r, $r), got $(size(Nijk))")

    one_vec = zeros(Int, r)
    one_vec[1] = 1

    local C, eqs_raw
    try
        C, eqs_raw = TensorCategories.pentagon_equations(Nijk, one_vec)
    catch err
        if err isa ArgumentError && occursin("empty collection", sprint(showerror, err))
            R0, _ = polynomial_ring(QQ, 0, :x)
            return R0, elem_type(R0)[], 0
        end
        rethrow()
    end
    eqs = filter(eq -> !(eq isa Integer) && !iszero(eq), eqs_raw)

    if isempty(eqs)
        R0, _ = polynomial_ring(QQ, 0, :x)
        return R0, elem_type(R0)[], 0
    end

    R = parent(eqs[1])
    n = nvars(R)
    return R, eqs, n
end

function _multiplicity_free_fusion(Nijk::Array{Int,3})
    return all(x -> x == 0 || x == 1, Nijk)
end

function _pentagon_variable_metadata(Nijk::Array{Int,3}, r::Int, nvars_total::Int)
    metadata = Vector{NamedTuple}(undef, nvars_total)
    var_stack = collect(1:nvars_total)
    for i in 1:r, j in 1:r, k in 1:r, o in 1:r
        (i == 1 || j == 1 || k == 1) && continue
        left = [a for a in 1:r if Nijk[i, j, a] != 0 && Nijk[a, k, o] != 0]
        right = [b for b in 1:r if Nijk[j, k, b] != 0 && Nijk[i, b, o] != 0]
        for col in eachindex(right), row in eachindex(left)
            var_idx = pop!(var_stack)
            metadata[var_idx] = (i = i, j = j, k = k, o = o,
                                 a = left[row], b = right[col],
                                 row = row, col = col)
        end
    end
    isempty(var_stack) || error("pentagon metadata length mismatch: $(length(var_stack)) unassigned variables")
    return metadata
end

function _gauge_channel_basis(Nijk::Array{Int,3}, r::Int)
    channels = Tuple{Int,Int,Int}[]
    for i in 1:r, j in 1:r, k in 1:r
        Nijk[i, j, k] == 0 && continue
        (i == 1 || j == 1) && continue
        push!(channels, (i, j, k))
    end
    return Dict(ch => idx for (idx, ch) in enumerate(channels))
end

function _gauge_weight(meta, channel_index::Dict{Tuple{Int,Int,Int}, Int})
    w = zeros(Int, length(channel_index))
    for (ch, sgn) in (((meta.i, meta.j, meta.a), 1),
                      ((meta.a, meta.k, meta.o), 1),
                      ((meta.j, meta.k, meta.b), -1),
                      ((meta.i, meta.b, meta.o), -1))
        idx = get(channel_index, ch, 0)
        idx != 0 && (w[idx] += sgn)
    end
    return w
end

function _rank_int_rows(rows::Vector{Vector{Int}}, ncols::Int)
    isempty(rows) && return 0
    return rank(matrix(QQ, length(rows), ncols, vcat(rows...)))
end

function _select_pentagon_gauge_fixed_indices(Nijk::Array{Int,3}, r::Int, nvars_total::Int)
    _multiplicity_free_fusion(Nijk) || return Int[]
    metadata = _pentagon_variable_metadata(Nijk, r, nvars_total)
    channel_index = _gauge_channel_basis(Nijk, r)
    isempty(channel_index) && return Int[]

    selected = Int[]
    rows = Vector{Int}[]
    current_rank = 0
    for var_idx in 1:nvars_total
        w = _gauge_weight(metadata[var_idx], channel_index)
        any(!iszero, w) || continue
        trial_rows = vcat(rows, [w])
        trial_rank = _rank_int_rows(trial_rows, length(channel_index))
        if trial_rank > current_rank
            push!(selected, var_idx)
            push!(rows, w)
            current_rank = trial_rank
        end
    end
    return selected
end
