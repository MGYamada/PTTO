# ============================================================
#  Step 1: Build block-diagonal atomic data from a stratum
# ============================================================

"""
    build_block_diagonal(catalog::Vector{AtomicIrrep}, stratum::Stratum)
        -> (S_atomic, T_atomic, K, blocks)

Given a stratum (catalog index → multiplicity) and an atomic catalog,
construct the block-diagonal atomic modular data as matrices over a
common cyclotomic field K = Q(ζ_N).

Returns:
- S_atomic: block-diagonal S matrix over K (size r × r where r = stratum.total_dim)
- T_atomic: diagonal T as a Vector (length r) over K
- K: the shared cyclotomic field
- blocks: Vector of (start_idx, end_idx, catalog_idx, copy_num) tuples
  describing the block structure

Currently assumes all atomic irreps in the catalog use the SAME cyclotomic
field K. This is true when `build_atomic_catalog(N)` is used consistently.
"""
function build_block_diagonal(catalog::Vector{AtomicIrrep}, stratum::Stratum)
    # Find the common cyclotomic field
    K = nothing
    for (idx, _) in stratum.multiplicities
        if K === nothing
            K = catalog[idx].K
        else
            K === catalog[idx].K || error("Atomic irreps use different cyclotomic fields; not supported yet")
        end
    end
    K === nothing && error("Empty stratum")

    r = stratum.total_dim

    # Build block-diagonal S as r × r zero matrix, then fill blocks
    S_atomic = zero_matrix(K, r, r)
    T_atomic = Vector{elem_type(K)}(undef, r)
    blocks = Tuple{Int, Int, Int, Int}[]

    pos = 1
    # Sort for determinism
    for (idx, m) in sort(collect(stratum.multiplicities); by = first)
        atom = catalog[idx]
        d = atom.dim
        for copy_num in 1:m
            # Place atom.S at positions (pos:pos+d-1, pos:pos+d-1)
            for i in 1:d
                for j in 1:d
                    S_atomic[pos + i - 1, pos + j - 1] = atom.S[i, j]
                end
            end
            for i in 1:d
                T_atomic[pos + i - 1] = atom.T[i, i]
            end
            push!(blocks, (pos, pos + d - 1, idx, copy_num))
            pos += d
        end
    end

    return S_atomic, T_atomic, K, blocks
end


# ============================================================
#  Step 3: T-eigenspace decomposition
# ============================================================

"""
    t_eigenspace_decomposition(T_Fp::Vector{Int}, p::Int)
        -> Dict{Int, Vector{Int}}

Given T as a vector of F_p eigenvalues (on the basis), return a Dict
mapping each distinct eigenvalue to the list of indices (1-based) where
it appears.
"""
function t_eigenspace_decomposition(T_Fp::Vector{Int}, p::Int)
    groups = Dict{Int, Vector{Int}}()
    for (i, val) in enumerate(T_Fp)
        if haskey(groups, val)
            push!(groups[val], i)
        else
            groups[val] = [i]
        end
    end
    return groups
end

"""
    parameter_dim(eigenspaces::Dict{Int, Vector{Int}}) -> Int

Compute parameter_dim = Σ_θ C(n_θ, 2) where n_θ = |eigenspaces[θ]|.
This is the continuous dimension of the block-U moduli.
"""
function parameter_dim(eigenspaces::Dict{Int, Vector{Int}})
    total = 0
    for (_, indices) in eigenspaces
        n = length(indices)
        total += div(n * (n - 1), 2)
    end
    return total
end

"""
    apply_o2_block(S::Matrix{Int}, idx_pair::Tuple{Int, Int},
                   u::Int, v::Int, det_sign::Int, p::Int) -> Matrix{Int}

Apply a 2×2 O(2) transformation to S in the 2-dimensional subspace
indexed by (i, j) = idx_pair, leaving all other rows/columns unchanged.

The transformation is U^T · S · U where U is the 2×2 block:
- det_sign = +1: U = [[u, -v], [v, u]]   (rotation)
- det_sign = -1: U = [[u,  v], [v, -u]]  (reflection)

In both cases (u, v) satisfies u² + v² ≡ 1 (mod p).
"""
function apply_o2_block(S::Matrix{Int}, idx_pair::Tuple{Int, Int},
                        u::Int, v::Int, det_sign::Int, p::Int)
    n = size(S, 1)
    U = Matrix{Int}(I_n(n))
    (i, j) = idx_pair
    if det_sign == +1
        U[i, i] = u
        U[i, j] = mod(-v, p)
        U[j, i] = v
        U[j, j] = u
    elseif det_sign == -1
        U[i, i] = u
        U[i, j] = v
        U[j, i] = v
        U[j, j] = mod(-u, p)
    else
        error("det_sign must be ±1")
    end
    UT = transpose_mod(U, p)
    return matmul_mod(matmul_mod(UT, S, p), U, p)
end

"""
    apply_block_U(S::Matrix{Int}, indices::Vector{Int},
                  U_block::Matrix{Int}, p::Int) -> Matrix{Int}

Apply an n×n block transformation U_block to S in the n-dimensional
subspace indexed by `indices` (length must equal size(U_block, 1)).
Other positions are left unchanged. Returns U^T · S · U (mod p), where
U is the full r×r matrix that is identity off the indexed block and
U_block on it.
"""
function apply_block_U(S::Matrix{Int}, indices::Vector{Int},
                       U_block::Matrix{Int}, p::Int)
    n_block = length(indices)
    size(U_block) == (n_block, n_block) ||
        error("U_block size $(size(U_block)) doesn't match |indices|=$n_block")

    r = size(S, 1)
    U_full = I_n(r)
    for ii in 1:n_block
        for jj in 1:n_block
            U_full[indices[ii], indices[jj]] = U_block[ii, jj]
        end
    end
    UT = transpose_mod(U_full, p)
    return matmul_mod(matmul_mod(UT, S, p), U_full, p)
end

function _apply_block_U_product(S_atomic::Matrix{Int},
                                block_choices,
                                p::Int)
    S_prime = copy(S_atomic)
    for choice in block_choices
        S_prime = apply_block_U(S_prime, choice.indices, choice.U, p)
    end
    return S_prime
end

function _orthogonal_block_product_fixes_S(S::Matrix{Int}, degenerate, p::Int)
    r = size(S, 1)
    for (_, indices) in degenerate
        idx_set = Set(indices)
        lambda = mod(S[indices[1], indices[1]], p)
        for i in indices, j in indices
            expected = i == j ? lambda : 0
            mod(S[i, j] - expected, p) == 0 || return false
        end
        for i in indices, j in 1:r
            j in idx_set && continue
            mod(S[i, j], p) == 0 || return false
            mod(S[j, i], p) == 0 || return false
        end
    end
    return true
end
