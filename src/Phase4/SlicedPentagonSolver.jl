"""
    SlicedPentagonSolver

Wire Kitaev/Hodge gauge slice into the pentagon HC solver.

Given a fusion rule and a prior F-symbol solution (used as a base point
for linearization), builds:
  1. The KitaevComplex F-coordinate gauge analysis,
  2. A mapping from KitaevComplex FKeys to TensorCategories pentagon
     variables x_1, ..., x_n,
  3. Linear slice constraints expressed in pentagon variables,
  4. The augmented polynomial system passed to HomotopyContinuation.

The slice kills the continuous gauge orbit (dim 1 for Ising, dim 1 for
Fibonacci), reducing the HC mixed volume correspondingly.

This module depends on:
  - KitaevComplex: for slice_basis
  - PentagonEquations: for TensorCategories pentagon system
  - HexagonEquations: for the FKey ↔ x_i mapping (mirrors `assign_F_to_associator!`)
"""
module SlicedPentagonSolver

using LinearAlgebra
using Oscar
using TensorCategories
using ACMG: FusionRule

export build_fkey_to_xvar_map
export slice_constraints_as_polynomials
export get_sliced_pentagon_system
export solve_pentagon_with_slice

"""
    build_fkey_to_xvar_map(Nijk, r, one_vec) -> Dict{NTuple{6,Int}, Int}

Build a Dict mapping F-symbol key (a,b,c,d,e,f) (1-based) to pentagon
variable index i such that x_i corresponds to F^{abc}_{d;e,f}.

Follows the exact traversal used by `assign_F_to_associator!`:
  - Nested loops over (i, j, k, o) ∈ 1:m × 1:m × 1:m × 1:m
  - Skip blocks where i, j, or k is the unit
  - Block at (i,j,k,o) has size (r_block, t_block)
  - Entries filled in row-major order into a stack, then pop! from the end
    during assignment

This produces the correspondence:
  variable x_p  ↔  F-symbol entry at some (i,j,k,o) block at (a, b)

In TensorCategories' convention, for a block at (i,j,k,o):
  - rows = paths (i⊗j → e, e⊗k → o) indexed by e ∈ i⊗j with e⊗k → o
  - cols = paths (j⊗k → f, i⊗f → o) indexed by f ∈ j⊗k with i⊗f → o
  - entry (row_idx_of_e, col_idx_of_f) = F^{ijk}_{o; e, f}

The channel indices e and f are enumerated in their natural 1:m order
(only allowed values included).
"""
function build_fkey_to_xvar_map(Nijk::Array{Int,3}, r::Int, one_vec::Vector{Int})
    m = r
    # Simulate the same traversal. We'll assign each block's (row, col) entry
    # a pentagon-variable index, then apply the pop! convention.

    # Step 1: for each (i,j,k,o), enumerate the row-basis e's and col-basis f's.
    block_info = Tuple{Int,Int,Int,Int,Vector{Int},Vector{Int}}[]  # (i,j,k,o, rows_e, cols_f)
    for i in 1:m, j in 1:m, k in 1:m, o in 1:m
        sum(one_vec[[i, j, k]]) > 0 && continue
        # Rows: e with N^{ij}_e ≥ 1 AND N^{ek}_o ≥ 1
        rows_e = Int[]
        for e in 1:m
            if Nijk[i, j, e] ≥ 1 && Nijk[e, k, o] ≥ 1
                push!(rows_e, e)
            end
        end
        # Cols: f with N^{jk}_f ≥ 1 AND N^{if}_o ≥ 1
        cols_f = Int[]
        for f in 1:m
            if Nijk[j, k, f] ≥ 1 && Nijk[i, f, o] ≥ 1
                push!(cols_f, f)
            end
        end
        isempty(rows_e) && continue
        isempty(cols_f) && continue
        push!(block_info, (i, j, k, o, rows_e, cols_f))
    end

    # Step 2: Build a flat list of (fkey = (i,j,k,o,e,f)) in ROW-MAJOR order per block,
    # ordered by block traversal. This is the order entries were POPPED.
    flat_fkeys = NTuple{6,Int}[]
    for (i, j, k, o, rows_e, cols_f) in block_info
        for e in rows_e, f in cols_f
            # row-major: entries[(a-1)*t + b] = M[a, b]
            # push in row-major: iterate row a (= index of e), then col b (= index of f)
            push!(flat_fkeys, (i, j, k, o, e, f))
        end
    end

    # Step 3: The pop! convention.  `assign_F_to_associator!` does:
    #   y = copy(F_values)
    #   for each block in traversal order:
    #     for each entry in the block (in row-major order, 1..r*t):
    #       entry = pop!(y)
    #
    # So the FIRST entry assigned (first block, first row-major position)
    # equals y[end] = F_values[end].
    #
    # Hence  flat_fkeys[1]  ↔  x_{length(F_values)}
    #        flat_fkeys[2]  ↔  x_{length(F_values) - 1}
    #        ...
    #        flat_fkeys[n]  ↔  x_1
    n = length(flat_fkeys)
    fkey_to_xvar = Dict{NTuple{6,Int}, Int}()
    for (pos, fkey) in enumerate(flat_fkeys)
        fkey_to_xvar[fkey] = n - pos + 1
    end
    return fkey_to_xvar
end

"""
    slice_constraints_as_polynomials(ga::GaugeAnalysis,
                                     fkey_to_xvar::Dict,
                                     R::Oscar.MPolyRing) -> Vector

Given the result of `analyze_gauge(fcs)`, translate each column of
`ga.Delta_gauge_eff` (a vector in F-coordinate space) into a linear
polynomial in the pentagon variables x_1, …, x_n. Setting this
polynomial = 0 kills that gauge direction.

Returns `Vector{elem_type(R)}` of `ga.gauge_orbit_dim` linear polynomials.

The polynomial for gauge direction `v ∈ R^nF` is:
    Σ_k v[k] * (x_{pent_var_idx} - F_value_at_key)

where key = ga.fcs.vars[k] and pent_var_idx comes from fkey_to_xvar.
We subtract the base F-value at each position so that the polynomial
vanishes AT the base F-solution (which sits on the slice by construction).

Note: entries with no corresponding pentagon variable (e.g. fixed by
TensorCategories' unit-axiom filtering) are dropped from the sum —
they contribute only a constant, which cancels at the base point.
"""
function slice_constraints_as_polynomials(ga::KitaevComplex.GaugeAnalysis,
                                          fkey_to_xvar::Dict,
                                          R)
    fcs = ga.fcs
    xs = gens(R)
    polys = Any[]

    # Eff gauge has ga.gauge_orbit_dim columns spanning the effective gauge image.
    # Use SVD of Delta_gauge_eff to get a basis of im Δ_gauge_eff.
    Dg_eff = ga.Delta_gauge_eff
    F = svd(Dg_eff; full = true)
    rank_g = ga.gauge_orbit_dim
    # Left singular vectors for nonzero singular values span im Δ_gauge_eff
    gauge_basis = F.U[:, 1:rank_g]   # nF × rank_g

    for col in 1:rank_g
        v = gauge_basis[:, col]
        poly = zero(R)
        for (k, key) in enumerate(fcs.vars)
            coeff = v[k]
            abs(coeff) < 1e-12 && continue
            if haskey(fkey_to_xvar, key)
                pent_idx = fkey_to_xvar[key]
                Fval = KitaevComplex.F_value(fcs, key)
                # Constraint: v · (x - F_base) = 0, so x is allowed to move
                # only perpendicular to v from the base point.
                poly += coeff * (xs[pent_idx] - Fval)
            end
            # Entries not in fkey_to_xvar correspond to unit-axiom-fixed
            # positions (TensorCategories doesn't generate variables for them).
            # These contribute a constant v[k] * (F_val_const - F_val_const) = 0
            # under the assumption that at base point x = F, so we can safely skip.
        end
        push!(polys, poly)
    end
    return polys
end

"""
    get_sliced_pentagon_system(Nijk, r, base_F_func) -> (R, eqs_augmented, n, n_slice)

Assemble the augmented pentagon system with Kitaev slice constraints.

Arguments
- `Nijk`: fusion array
- `r`: rank
- `base_F_func`: function (a,b,c,d,e,f) -> Float64 giving a base F-solution
  (used as the linearization point for the Kitaev gauge analysis)

Returns
- `R`: polynomial ring
- `eqs_augmented`: pentagon equations + slice constraints
- `n`: number of pentagon variables
- `n_slice`: number of slice constraints added (= gauge_orbit_dim)
"""
function get_sliced_pentagon_system(Nijk::Array{Int,3}, r::Int,
                                     base_F_func::Function)
    # Pentagon system from TensorCategories
    R, eqs, n = PentagonEquations.get_pentagon_system(Nijk, r)

    # KitaevComplex analysis
    fr = FusionRule(Nijk)
    fcs = KitaevComplex.build_F_coord_space(fr, base_F_func)
    ga  = KitaevComplex.analyze_gauge(fcs)

    # Build unit vector: unit is object 1
    one_vec = zeros(Int, r)
    one_vec[1] = 1

    # F-key ↔ pentagon-variable mapping
    fkey_to_xvar = build_fkey_to_xvar_map(Nijk, r, one_vec)

    # Slice constraints
    slice_polys = slice_constraints_as_polynomials(ga, fkey_to_xvar, R)

    eqs_aug = vcat(eqs, slice_polys)
    return R, eqs_aug, n, length(slice_polys)
end

end # module SlicedPentagonSolver
