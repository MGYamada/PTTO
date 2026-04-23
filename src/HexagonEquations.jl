"""
    HexagonEquations

Hexagon equation generation with F-symbols fixed to numerical values.

The strategy avoids symbolic matrix inversion (which is unreliable on
polynomial coefficients) by treating the forward braiding `c` and the
reverse braiding `c'` as **independent** variable families, with three
kinds of equations:

  (A) Left hexagon, using forward braiding c and forward associator α:
      (c_{X,Y} ⊗ id_Z) ∘ α_{Y,X,Z} ∘ (id_Y ⊗ c_{X,Z})
          = α_{Y,Z,X} ∘ c_{X, Y⊗Z} ∘ α_{X,Y,Z}

  (B) Left-hexagon SHAPE in a second category with α' = α⁻¹ and
      reverse braiding c'. Equivalent to the standard right hexagon with
      original α and c⁻¹.

  (C) Inverse consistency: R^{ij}_k · S^{ji}_k = I, enforcing that c ∘ c'
      is the identity.

α⁻¹ is precomputed numerically block-by-block via `invert_associator_numeric`.

The resulting polynomial system has `2 · r_var_count` variables over
AcbField() (Arb complex field). Solved by `HexagonSolver` via
HomotopyContinuation.

Depends on: Oscar, TensorCategories.
"""

using Oscar
using TensorCategories


"""
    _coerce_complex(c, K) -> K element

Convert a Julia `Number` (typically ComplexF64) into an element of the ring K.
If K is an AcbField, goes through real/imag parts. If K is a polynomial ring
over AcbField, first goes to the scalar AcbField then lifts to K.
"""
function _coerce_complex(c::Number, K)
    scalar = K isa AcbField ? K : base_ring(K)
    # Build scalar(real) + i * scalar(imag) in the AcbField
    re = scalar(Float64(real(c)))
    im_part = scalar(Float64(imag(c)))
    acb_val = re + onei(scalar) * im_part
    # Lift into K (identity if K is already the scalar field)
    return K(acb_val)
end

"""
    invert_associator_numeric(F_values, mult, one) -> Vector{ComplexF64}

Compute the inverse F-symbol values. Given `F_values` that would produce the
associator via `assign_F_to_associator!`, return a vector `F_inv_values` such
that `assign_F_to_associator!(_, F_inv_values)` produces the inverse
associator (block-wise inverse over ComplexF64).

Implementation: simulate `assign_F_to_associator!` exactly to extract each
block in-order, invert, and rebuild a fresh vector that would regenerate
those inverted blocks under the same `assign` convention.
"""
function invert_associator_numeric(F_values::Vector{ComplexF64},
                                   mult::Array{Int,3}, one_vec::Vector{Int})
    dummy = TensorCategories.six_j_category(QQ, mult)
    dummy.one = one_vec
    m = dummy.simples

    # First pass: extract blocks in the same order `assign_F_to_associator!` uses.
    # Each block is built from `entries = [pop!(y) for _ in 1:(r*t)]` and then
    # `matrix(K, r, t, entries)` fills row-major.
    y_stack = copy(F_values)
    extracted_blocks = Vector{Tuple{Int,Int,Matrix{ComplexF64}}}()  # (r, t, M)
    for i in 1:m, j in 1:m, k in 1:m, o in 1:m
        sum(one_vec[[i, j, k]]) > 0 && continue
        (r, t) = size(dummy.ass[i, j, k, o])
        @assert r == t "F-block at ($i,$j,$k,$o) is $r×$t — not square"
        # Mirror assign exactly:
        entries = ComplexF64[pop!(y_stack) for _ in 1:(r * t)]
        # `matrix(K, r, t, entries)` in Oscar fills row-major, so:
        # M[a, b] = entries[(a-1)*t + b]
        M = Matrix{ComplexF64}(undef, r, t)
        for a in 1:r, b in 1:t
            M[a, b] = entries[(a - 1) * t + b]
        end
        push!(extracted_blocks, (r, t, M))
    end
    @assert isempty(y_stack) "F_values had leftover entries"

    # Invert each block
    inverted_blocks = [(r, t, inv(M)) for (r, t, M) in extracted_blocks]

    # Rebuild a vector that, when passed to assign_F_to_associator!, reproduces
    # the inverted blocks. We need to "un-pop" the pop! sequence.
    # For each block (r, t, Minv):
    #   entries_new = [Minv[a,b] in row-major order] = flat_rowmajor(Minv)
    #   These entries must appear in F_inv_values such that `pop!` consumes them
    #   in REVERSE, i.e., F_inv_values[end] -> entries_new[1], etc.
    # So for this block alone, we append `reverse(entries_new)` to the OUTPUT
    # END... but blocks are processed in-order and each block is placed at
    # increasingly deeper positions. Concretely:
    # At start, y_stack = F_values. First block pops r1*t1 from the END.
    # → Those entries are F_values[end-r1*t1+1 : end] (in original order).
    # In pop order: F_values[end], F_values[end-1], ..., F_values[end-r1*t1+1]
    # become entries_new[1], entries_new[2], ..., entries_new[r1*t1].
    # So F_values[end - (p-1)] = entries_new[p]  for p = 1..r1*t1
    # Therefore to ENCODE entries_new into a slice of F_inv_values at the end:
    # F_inv_values[end - (p-1)] = inv_entries[p]  (same relation for the inverted block)

    F_inv = Vector{ComplexF64}(undef, length(F_values))
    cursor_end = length(F_values)
    for (r, t, Minv) in inverted_blocks
        inv_entries = Vector{ComplexF64}(undef, r * t)
        for a in 1:r, b in 1:t
            inv_entries[(a - 1) * t + b] = Minv[a, b]
        end
        # Place into F_inv[cursor_end - r*t + 1 : cursor_end] such that
        # F_inv[cursor_end - (p-1)] = inv_entries[p]
        for p in 1:(r * t)
            F_inv[cursor_end - (p - 1)] = inv_entries[p]
        end
        cursor_end -= r * t
    end
    @assert cursor_end == 0
    return F_inv
end

"""
    assign_F_to_associator!(poly_C, F_values)

Write a ComplexF64 F-symbol solution into `poly_C.ass` following the **same
traversal order** used by `TensorCategories.pentagon_equations`, so that the
i-th component of `F_values` lands in the slot that was called `x_i` during
Level II. The `pop!` convention matches `pentagon_equations` (pops from the
end of a copy of the variable vector).
"""
function assign_F_to_associator!(poly_C, F_values::Vector{<:Number})
    K = base_ring(poly_C)
    m = poly_C.simples
    one_vec = poly_C.one
    y = copy(F_values)     # treat as a stack; pop! takes from the end
    for i in 1:m, j in 1:m, k in 1:m, o in 1:m
        sum(one_vec[[i, j, k]]) > 0 && continue
        (r, t) = size(poly_C.ass[i, j, k, o])
        entries = [_coerce_complex(pop!(y), K) for _ in 1:(r * t)]
        poly_C.ass[i, j, k, o] = matrix(K, r, t, entries)
    end
    @assert isempty(y) "F_values length ($(length(F_values))) does not match the expected number of F-symbol slots"
    return poly_C
end

"""
    _number_of_variables_in_hexagon_equations(poly_C) -> Int

Count independent R-symbol entries. For each (i, j), the braiding lives in
`C.braiding[i, j, k]` as a matrix whose size is determined by multiplicity.
We count the total number of matrix entries across all (i, j, k) with
N_{ij}^k > 0.

Multiplicity-free case: each nonzero fusion coefficient contributes 1.
"""
function _number_of_variables_in_hexagon_equations(poly_C)
    m = poly_C.simples
    mult = poly_C.tensor_product  # N_{ij}^k  (TensorCategories stores this as tensor_product)
    n_r = 0
    for i in 1:m, j in 1:m, k in 1:m
        N_ijk = mult[i, j, k]
        N_ijk == 0 && continue
        # Multiplicity-free: N_ijk is 0 or 1, contributes 1
        # General case: would contribute N_ijk^2 (full N×N block)
        n_r += N_ijk * N_ijk
    end
    return n_r
end

"""
    hexagon_equations(mult, one, F_values) -> (poly_C_fwd, poly_C_rev, eqs)

Build the hexagon equations with F-symbols fixed to numerical values. To
avoid symbolic matrix inversion (which fails on polynomial coefficients),
we treat the forward braiding `c` and the reverse braiding `c'` as
INDEPENDENT variable families, with three kinds of equations:

  (A) Left hexagon, using forward braiding c and forward associator α:
      (c_{X,Y} ⊗ id_Z) ∘ α_{Y,X,Z} ∘ (id_Y ⊗ c_{X,Z})
          = α_{Y,Z,X} ∘ c_{X, Y⊗Z} ∘ α_{X,Y,Z}

  (B) Right hexagon, rewritten with reverse braiding c' and inverse
      associator α⁻¹ (both as separate data):
      (id_Z ⊗ c'_{X,Y}) ∘ α'_{Z,X,Y} ∘ (c'_{X,Z} ⊗ id_Y)
          = α'_{X,Y,Z} ∘ c'_{X⊗Y, Z} ∘ α'_{Y,Z,X}
      where α' = α⁻¹ is precomputed NUMERICALLY (no symbolic inv).

  (C) Inverse consistency: c_{X,Y} ∘ c'_{X,Y} = id  (equivalently c ⋅ c' = 1
      for each multiplicity-free R-symbol component).

Returns the two categories (forward and reverse) and the combined equation
list. The polynomial ring has 2 * r_var_count variables (first half = r,
second half = s = reverse braiding).
"""
function hexagon_equations(mult::Array{Int,3}, one_vec::Vector{Int},
                           F_values::Vector{ComplexF64})
    m = size(mult, 1)

    # Sanity check F_values length
    _dummy = TensorCategories.six_j_category(QQ, mult)
    _dummy.one = one_vec
    expected_n_F = TensorCategories._number_of_variables_in_pentagon_equations(_dummy)
    @assert length(F_values) == expected_n_F "F_values length $(length(F_values)) != expected $expected_n_F"

    # Precompute F⁻¹ values numerically (block-wise)
    Finv_values = invert_associator_numeric(F_values, mult, one_vec)

    # Count R variables (using the dummy category for consistent structure)
    _dummy_cc = TensorCategories.six_j_category(AcbField(), mult)
    _dummy_cc.one = one_vec
    r_var_count = _number_of_variables_in_hexagon_equations(_dummy_cc)

    # Build a single polynomial ring with 2 * r_var_count variables:
    # first half -> forward braiding (r), second half -> reverse braiding (s)
    R_ring, xs = polynomial_ring(AcbField(), 2 * r_var_count)
    r_vars = xs[1:r_var_count]
    s_vars = xs[r_var_count + 1 : 2 * r_var_count]

    # Two categories sharing the same R_ring:
    # poly_C_fwd: associator = F, braiding = r
    # poly_C_rev: associator = F⁻¹, braiding = s
    poly_C_fwd = TensorCategories.six_j_category(R_ring, mult)
    poly_C_fwd.one = one_vec
    assign_F_to_associator!(poly_C_fwd, F_values)

    poly_C_rev = TensorCategories.six_j_category(R_ring, mult)
    poly_C_rev.one = one_vec
    assign_F_to_associator!(poly_C_rev, Finv_values)

    # Fill braidings
    function _fill_braiding!(poly_C, vars)
        m_ = poly_C.simples
        braid_arr = Array{MatElem, 3}(undef, m_, m_, m_)
        y = copy(vars)  # stack
        for i in 1:m_, j in 1:m_, k in 1:m_
            N_ijk = mult[i, j, k]
            if N_ijk == 0
                braid_arr[i, j, k] = zero_matrix(R_ring, 0, 0)
                continue
            end
            entries = [pop!(y) for _ in 1:(N_ijk * N_ijk)]
            braid_arr[i, j, k] = matrix(R_ring, N_ijk, N_ijk, entries)
        end
        @assert isempty(y)
        TensorCategories.set_braiding!(poly_C, braid_arr)
    end
    _fill_braiding!(poly_C_fwd, r_vars)
    _fill_braiding!(poly_C_rev, s_vars)

    # Write down the equations categorically — NO symbolic inv is used.
    eqs = elem_type(R_ring)[]
    Ss_fwd = simples(poly_C_fwd)
    Ss_rev = simples(poly_C_rev)

    # (A) Left hexagon in poly_C_fwd
    for X in Ss_fwd, Y in Ss_fwd, Z in Ss_fwd
        lhs = (braiding(X, Y) ⊗ id(Z)) ∘ associator(Y, X, Z) ∘ (id(Y) ⊗ braiding(X, Z))
        rhs = associator(Y, Z, X) ∘ braiding(X, Y ⊗ Z) ∘ associator(X, Y, Z)
        append!(eqs, collect(matrix(lhs - rhs))[:])
    end

    # (B) Left hexagon SHAPE in poly_C_rev (which has α' = α⁻¹ and c' = reverse braiding).
    # The RIGHT hexagon of a braided category, when you treat c' and α' as new forward
    # data, becomes a left hexagon in the "reverse" category. Concretely:
    #   (c'_{X,Y} ⊗ id_Z) ∘ α'_{Y,X,Z} ∘ (id_Y ⊗ c'_{X,Z})
    #       = α'_{Y,Z,X} ∘ c'_{X, Y⊗Z} ∘ α'_{X,Y,Z}
    # This is equivalent to the standard right hexagon with original α and c⁻¹.
    for X in Ss_rev, Y in Ss_rev, Z in Ss_rev
        lhs = (braiding(X, Y) ⊗ id(Z)) ∘ associator(Y, X, Z) ∘ (id(Y) ⊗ braiding(X, Z))
        rhs = associator(Y, Z, X) ∘ braiding(X, Y ⊗ Z) ∘ associator(X, Y, Z)
        append!(eqs, collect(matrix(lhs - rhs))[:])
    end

    # (C) Inverse consistency: c_{X,Y} ∘ c'_{X,Y} = id_{X⊗Y}
    # Using the forward-category objects, braiding(X,Y) goes X⊗Y -> Y⊗X.
    # For c ∘ c' to make sense as X⊗Y -> X⊗Y, we compose:
    #   (c'_{Y,X} from poly_C_rev) ∘ (c_{X,Y} from poly_C_fwd) : X⊗Y -> Y⊗X -> X⊗Y
    # but that requires moving between the two categories, which share R_ring
    # but are different SixJCategory instances. We instead impose the equation
    # DIRECTLY at the level of R-symbols and S-symbols:
    # For each (i, j, k) with N_{ij}^k > 0:
    #   R^{ij}_k · S^{ji}_k = 1   (multiplicity-free case; treating each as a scalar)
    # Or in general as matrix equation: R^{ij}_k · S^{ji}_k = I
    # Rationale: c_{X_i, X_j}: X_i ⊗ X_j -> X_j ⊗ X_i has component X_k equal to
    # R^{ij}_k; the reverse braiding c'_{X_j, X_i}: X_j ⊗ X_i -> X_i ⊗ X_j has
    # component X_k equal to S^{ji}_k. Their composite on the X_k component
    # must be identity.
    #
    # We need to index r_vars and s_vars consistently. Recall _fill_braiding!
    # assigns variables by popping; we reproduce that indexing here.
    function _var_index(i::Int, j::Int, k::Int, a::Int, b::Int,
                        var_pool_size::Int, m_::Int)
        # Compute the position of entry (a, b) of braid[i,j,k] inside the vars
        # used by _fill_braiding!. Since pop! consumes from the end, the LAST
        # assigned block is (m_, m_, m_) and within it entries[1] = pop = vars[end].
        # General formula: count how many (i', j', k', a', b') come AFTER (i,j,k,a,b)
        # in the assignment order (= pop order).
        # Total positions BEFORE = count of earlier (i',j',k') blocks plus earlier
        # (a',b') within current block.
        # Easier: enumerate and memoize.
        error("use the precomputed lookup `_rs_index` instead")
    end

    # Build lookup: (i, j, k) -> (start_pos, size) for the R variable indices.
    rs_block_info = Dict{Tuple{Int,Int,Int}, Tuple{Int,Int}}()
    let pos = r_var_count
        for i in 1:m, j in 1:m, k in 1:m
            N_ijk = mult[i, j, k]
            N_ijk == 0 && continue
            sq = N_ijk * N_ijk
            # In pop order, the LAST block iterated is (m, m, m). It receives the
            # LAST sq variables (positions pos-sq+1 .. pos).
            # We walk in normal order and need the positions that pop! would assign.
            # Since pop! consumes from the end, the FIRST block in our loop
            # receives the LAST sq variables of the pool. But our loop order is
            # i in 1:m, j in 1:m, k in 1:m (nested), which means the LAST block
            # visited is (m, m, m), which pops FIRST from the remaining stack,
            # i.e., it gets the highest position.
            # Wait — `pop!(y)` where y = copy(vars) with vars=xs[1:r_var_count].
            # y starts as [x1, x2, ..., x_{r_var_count}]. pop! returns x_{r_var_count} first.
            # Loop order visits (1,1,1), (1,1,2), ..., (m,m,m). So (1,1,1) block
            # gets x_{r_var_count}, x_{r_var_count - 1}, ... (first sq entries in pop order).
            # So block (1,1,1) occupies positions r_var_count - sq + 1 .. r_var_count (but in reverse).
            # Hmm — this is getting complex. Let me just do it carefully via simulation.
        end
    end

    # Simulate _fill_braiding! indexing exactly to build lookup tables.
    # For forward (r_vars) positions within r_vars (1..r_var_count):
    r_block_positions = Dict{Tuple{Int,Int,Int}, Vector{Int}}()
    let y_positions = collect(1:r_var_count)  # positions in r_vars
        for i in 1:m, j in 1:m, k in 1:m
            N_ijk = mult[i, j, k]
            N_ijk == 0 && continue
            # Each block pops N_ijk*N_ijk positions from the end of y_positions.
            # entries[p] = pop! = y_positions[end - (p-1)] (original), then list shrinks.
            popped = Int[]
            for _ in 1:(N_ijk * N_ijk)
                push!(popped, pop!(y_positions))
            end
            # popped[p] is the r_vars-index that becomes entries[p], which fills
            # matrix[a,b] with (a-1)*N_ijk + b = p, i.e., row-major.
            # We store the list so we can map (a,b) -> global index later.
            r_block_positions[(i, j, k)] = popped
        end
        @assert isempty(y_positions)
    end

    # Same structure for s_vars (indices within s_vars are 1..r_var_count;
    # global index in xs is r_var_count + local_index).
    s_block_positions = Dict{Tuple{Int,Int,Int}, Vector{Int}}()
    let y_positions = collect(1:r_var_count)
        for i in 1:m, j in 1:m, k in 1:m
            N_ijk = mult[i, j, k]
            N_ijk == 0 && continue
            popped = Int[]
            for _ in 1:(N_ijk * N_ijk)
                push!(popped, pop!(y_positions))
            end
            s_block_positions[(i, j, k)] = popped
        end
        @assert isempty(y_positions)
    end

    # Now add the inverse-consistency equations R^{ij}_k · S^{ji}_k = I.
    # For multiplicity-free (N_ijk ∈ {0, 1}) this is simply r * s - 1 = 0.
    # In general: matrix R^{ij}_k (N_ijk × N_ijk) times S^{ji}_k (= N_jik × N_jik
    # = N_ijk × N_ijk since N_ijk = N_jik for symmetric tensor product, but we
    # should be careful).
    for i in 1:m, j in 1:m, k in 1:m
        N_ijk = mult[i, j, k]
        N_ijk == 0 && continue

        # Forward block r-matrix R^{ij}_k: lives in r_vars at positions r_block_positions[(i,j,k)]
        r_pos = r_block_positions[(i, j, k)]
        R_block = [r_vars[r_pos[(a - 1) * N_ijk + b]] for a in 1:N_ijk, b in 1:N_ijk]

        # Reverse block s-matrix S^{ji}_k: lives in s_vars at positions s_block_positions[(j,i,k)]
        # (note the swap of i and j)
        if !haskey(s_block_positions, (j, i, k))
            continue  # fusion rule is not symmetric for this (i,j,k)? skip defensively
        end
        s_pos = s_block_positions[(j, i, k)]
        N_jik = mult[j, i, k]
        @assert N_jik == N_ijk "Fusion multiplicity not symmetric: N^$k_{$i,$j}=$N_ijk vs N^$k_{$j,$i}=$N_jik"
        S_block = [s_vars[s_pos[(a - 1) * N_jik + b]] for a in 1:N_jik, b in 1:N_jik]

        # Product R_block * S_block = I
        for a in 1:N_ijk, c in 1:N_ijk
            expr = sum(R_block[a, b] * S_block[b, c] for b in 1:N_ijk)
            target = (a == c) ? one(R_ring) : zero(R_ring)
            push!(eqs, expr - target)
        end
    end

    return poly_C_fwd, poly_C_rev, filter(e -> !iszero(e), unique(eqs))
end

"""
    get_hexagon_system(Nijk, r, F_values) -> (R_ring, eqs, n_vars)

Wrapper analogous to `get_pentagon_system`. Returns polynomial ring, equation
list, and number of variables (= 2 * r_var_count, since both forward and
reverse braidings are variables).
"""
function get_hexagon_system(Nijk::Array{Int,3}, r::Int, F_values::Vector{<:Number})
    one_vec = zeros(Int, r)
    one_vec[1] = 1
    # F_values may arrive as Vector{Complex{Float64}}, make sure it's ComplexF64
    F_values_cf64 = ComplexF64.(F_values)
    poly_C_fwd, poly_C_rev, eqs = hexagon_equations(Nijk, one_vec, F_values_cf64)
    eqs_filt = filter(eq -> !(eq isa Integer) && !iszero(eq), eqs)
    @assert !isempty(eqs_filt) "Hexagon produced no equations — system may be trivially satisfied"
    R_ring = parent(eqs_filt[1])
    n_vars = nvars(R_ring)
    return R_ring, eqs_filt, n_vars
end

