"""
Phase 2: Block-U parametrisation and MTC reconstruction.

This module implements the core of the Phase 2 pipeline:

1. `build_block_diagonal`: given a stratum {m_λ} and atomic catalog,
   build the block-diagonal atomic (S, T) on V = ⊕_λ V_λ^{m_λ}.

2. `reduce_to_Fp`: reduce an Oscar Q(ζ_N)-matrix to F_p (requires N | p-1
   and a chosen primitive N-th root of unity in F_p).

3. `t_eigenspace_decomposition`: given T (diagonal over F_p), group
   indices by eigenvalue. Returns a Dict{Int, Vector{Int}} mapping
   each T-eigenvalue to the list of basis indices.

4. `sweep_O2`: for a 2-dimensional degenerate eigenspace with indices
   (i, j), enumerate the O(2)(F_p)-rational circle points (u, v) with
   u² + v² = 1, and for each compute S' = U^T S U.

5. `verlinde_check`: compute fusion coefficients N_{ij}^k from an F_p
   modular datum (S', T'), check all are small non-negative integers.

The top-level driver is `find_mtcs_at_prime(catalog, stratum, p)`,
which runs the whole pipeline at a single prime.

Follow-up functions (in a separate CRT module, TODO) will combine
results across multiple primes.
"""

using Oscar

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

# ============================================================
#  Step 4: O(2) sweep for 2-dimensional degenerate eigenspace
# ============================================================

"""
    o2_circle_points(p::Int) -> Vector{Tuple{Int, Int}}

Enumerate all (u, v) ∈ F_p × F_p with u² + v² ≡ 1 (mod p).
For p ≡ 1 (mod 4), there are 2(p-1) such pairs; else 2(p+1).

For small p (p < 500 or so) a brute enumeration is instant.
"""
function o2_circle_points(p::Int)
    pts = Tuple{Int, Int}[]
    for u in 0:(p-1)
        u2 = mod(u * u, p)
        target = mod(1 - u2, p)
        # Is `target` a QR? Find sqrt
        for v in 0:(p-1)
            if mod(v * v, p) == target
                push!(pts, (u, v))
            end
        end
    end
    return pts
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

# Local helpers (small utilities, will move to FpArith if reused broadly)
function I_n(n::Int)
    M = zeros(Int, n, n)
    for i in 1:n
        M[i, i] = 1
    end
    return M
end

function transpose_mod(A::Matrix{Int}, p::Int)
    n, m = size(A)
    B = zeros(Int, m, n)
    for i in 1:m
        for j in 1:n
            B[i, j] = A[j, i]
        end
    end
    return B
end

# ============================================================
#  General O(n) parametrisation (Cayley) — n ≥ 2
# ============================================================

"""
    cayley_so_n(params::Vector{Int}, n::Int, p::Int)
        -> Union{Nothing, Matrix{Int}}

Construct a matrix in SO(n)(F_p) from `params`, the C(n, 2) coefficients
of the upper-triangular part of an antisymmetric matrix A.

`params` is indexed by (i, j) with i < j in lex order:
    params[1] = A[1, 2]
    params[2] = A[1, 3]
    ...
    params[n-1] = A[1, n]
    params[n] = A[2, 3]
    ...
    params[C(n,2)] = A[n-1, n]

Returns `(I - A)(I + A)^{-1}` mod p.
Returns `nothing` if (I + A) is singular mod p (Cayley exceptional locus).
"""
function cayley_so_n(params::Vector{Int}, n::Int, p::Int)
    length(params) == div(n * (n - 1), 2) ||
        error("expected $(div(n*(n-1), 2)) params for n=$n; got $(length(params))")

    # Build antisymmetric A
    A = zeros(Int, n, n)
    k = 1
    for i in 1:n
        for j in (i+1):n
            A[i, j] = mod(params[k], p)
            A[j, i] = mod(-params[k], p)
            k += 1
        end
    end

    I_mat = I_n(n)
    IpA = [mod(I_mat[i, j] + A[i, j], p) for i in 1:n, j in 1:n]
    IpA_inv = inverse_mod_p(IpA, p)
    IpA_inv === nothing && return nothing
    ImA = [mod(I_mat[i, j] - A[i, j], p) for i in 1:n, j in 1:n]
    return matmul_mod(ImA, IpA_inv, p)
end

"""
    inverse_mod_p(M::Matrix{Int}, p::Int) -> Union{Nothing, Matrix{Int}}

Compute the inverse of an n×n matrix mod p via Gauss-Jordan elimination.
Returns `nothing` if M is singular mod p.
"""
function inverse_mod_p(M::Matrix{Int}, p::Int)
    n = size(M, 1)
    n == size(M, 2) || error("M must be square")

    # Augment [M | I]
    aug = zeros(Int, n, 2n)
    for i in 1:n
        for j in 1:n
            aug[i, j] = mod(M[i, j], p)
        end
        aug[i, n + i] = 1
    end

    # Gauss-Jordan
    for col in 1:n
        # Find pivot
        pivot_row = 0
        for r in col:n
            if aug[r, col] != 0
                pivot_row = r
                break
            end
        end
        pivot_row == 0 && return nothing  # singular

        if pivot_row != col
            # Swap rows
            for j in 1:(2n)
                aug[col, j], aug[pivot_row, j] = aug[pivot_row, j], aug[col, j]
            end
        end

        # Scale pivot row
        pivot = aug[col, col]
        pivot_inv = invmod(pivot, p)
        for j in 1:(2n)
            aug[col, j] = mod(aug[col, j] * pivot_inv, p)
        end

        # Eliminate other rows
        for r in 1:n
            r == col && continue
            factor = aug[r, col]
            factor == 0 && continue
            for j in 1:(2n)
                aug[r, j] = mod(aug[r, j] - factor * aug[col, j], p)
            end
        end
    end

    return aug[:, (n+1):(2n)]
end

"""
    enumerate_so_n_Fp(n::Int, p::Int) -> Vector{Matrix{Int}}

Enumerate all SO(n)(F_p) matrices via Cayley parametrisation. There are
~p^{C(n,2)} such matrices (minus the singular Cayley locus). Returns
a Vector of n×n Int matrices.

WARNING: count grows as p^{n(n-1)/2}. For n=2 this is ~p (small);
for n=3 it is ~p^3 (~10^5 for p=73); for n≥4 likely infeasible by
brute force.
"""
function enumerate_so_n_Fp(n::Int, p::Int)
    n >= 1 || error("n must be ≥ 1")
    n == 1 && return [reshape([1], 1, 1)]  # SO(1) = {1}
    dim_params = div(n * (n - 1), 2)
    results = Matrix{Int}[]
    for params_tuple in Iterators.product(ntuple(_ -> 0:(p-1), dim_params)...)
        params = collect(params_tuple)
        U = cayley_so_n(params, n, p)
        U === nothing && continue
        push!(results, U)
    end
    return results
end

"""
    enumerate_o_n_Fp(n::Int, p::Int) -> Vector{Matrix{Int}}

Enumerate all O(n)(F_p) matrices.
O(n) = SO(n) ⊔ (reflection · SO(n)). We use a single fixed reflection
(diag(-1, 1, 1, ..., 1)) and concatenate the two cosets.

Note: the SO(n) part comes from Cayley, but the Cayley map misses some
SO(n) elements (those where (I + A) is singular). For most uses (e.g. MTC
sweep), this small loss is acceptable; the missing locus is
codimension ≥ 1.
"""
function enumerate_o_n_Fp(n::Int, p::Int)
    so_part = enumerate_so_n_Fp(n, p)
    # Reflection R = diag(-1, 1, 1, ..., 1)
    R = I_n(n)
    R[1, 1] = p - 1  # = -1 mod p
    refl_part = [matmul_mod(R, U, p) for U in so_part]
    return vcat(so_part, refl_part)
end


"""
    enumerate_o_n_Fp_groebner(n::Int, p::Int) -> Vector{Matrix{Int}}

Algebraic-search entry point for O(n)(F_p) block candidates.

This computes/caches a Gröbner basis for the O(n) orthogonality ideal and
extracts F_p-rational points from that basis. It does not fall back to the
Cayley/reflection exhaustive sweep; use `search_mode = :exhaustive`
explicitly for that backend.
"""
function enumerate_o_n_Fp_groebner(n::Int, p::Int)
    gb_data = _ensure_orthogonality_groebner_cache!(n, p)
    gb_data !== nothing || return Matrix{Int}[]

    gb_points = _extract_orthogonality_points(gb_data, n, p)
    if gb_points !== nothing && !isempty(gb_points)
        return gb_points
    end
    return Matrix{Int}[]
end

const _ORTHOGONALITY_GB_CACHE = Dict{Tuple{Int, Int}, Any}()
const _VERLINDE_GB_CACHE = Dict{Tuple{Int, Int, UInt, Tuple{Vararg{Int}}}, Any}()

"""
    build_orthogonality_equations(vars, n::Int, p::Int)
        -> Vector

Build polynomial equations for `U^T U = I` over `F_p`, where `vars`
encodes an `n×n` matrix in row-major order.
"""
function build_orthogonality_equations(vars, n::Int, p::Int)
    eqs = Any[]
    idx(i, j) = (i - 1) * n + j
    for i in 1:n
        for j in i:n
            acc = zero(vars[1])
            for k in 1:n
                acc += vars[idx(k, i)] * vars[idx(k, j)]
            end
            rhs = (i == j) ? one(vars[1]) : zero(vars[1])
            push!(eqs, acc - rhs)
        end
    end
    return eqs
end

"""
    _ensure_orthogonality_groebner_cache!(n::Int, p::Int)

Attempt to compute and cache a Gröbner basis for the orthogonality ideal
of `O(n)` over `F_p`. This is the first building block for a future
F4/FGLM-style solver backend.
"""
function _ensure_orthogonality_groebner_cache!(n::Int, p::Int)
    key = (n, p)
    haskey(_ORTHOGONALITY_GB_CACHE, key) && return _ORTHOGONALITY_GB_CACHE[key]

    gb_data = nothing
    try
        F = GF(p)
        R, vars = polynomial_ring(F, n * n, :u)
        eqs = build_orthogonality_equations(vars, n, p)
        I = ideal(R, eqs)
        G = groebner_basis(I)
        gb_data = (ring = R, vars = vars, equations = eqs, ideal = I, gb = G)
    catch err
        @warn "Groebner preprocessing failed for O($n)(F_$p); returning no algebraic O(n) points" exception = (err, catch_backtrace())
    end

    _ORTHOGONALITY_GB_CACHE[key] = gb_data
    return gb_data
end

function _point_entry(point, var)
    # Common container styles returned by algebra systems:
    # 1) associative map keyed by variable object
    # 2) keyed by Symbol/String of the variable name
    # 3) custom object supporting getindex(point, var)
    try
        return point[var]
    catch
    end

    vname = string(var)
    try
        return point[vname]
    catch
    end
    try
        return point[Symbol(vname)]
    catch
    end

    return nothing
end

function _coerce_point_to_matrix(point, vars, n::Int, p::Int)
    vals = Int[]
    for v in vars
        entry = _point_entry(point, v)
        entry === nothing && return nothing
        local intval
        try
            intval = Int(entry)
        catch
            return nothing
        end
        push!(vals, mod(intval, p))
    end

    M = zeros(Int, n, n)
    t = 1
    for i in 1:n
        for j in 1:n
            M[i, j] = vals[t]
            t += 1
        end
    end
    return M
end

function _dedupe_matrices(mats::Vector{Matrix{Int}})
    seen = Set{String}()
    out = Matrix{Int}[]
    for M in mats
        key = join(vec(M), ",")
        if !(key in seen)
            push!(seen, key)
            push!(out, M)
        end
    end
    return out
end

function _sort_matrices_lex(mats::Vector{Matrix{Int}})
    return sort(mats; by = M -> join(vec(M), ","))
end

"""
    is_orthogonal_mod_p(U::Matrix{Int}, p::Int) -> Bool

Check `U^T U ≡ I (mod p)`.
"""
function is_orthogonal_mod_p(U::Matrix{Int}, p::Int)
    n, m = size(U)
    n == m || return false
    for i in 1:n
        for j in 1:n
            acc = 0
            for k in 1:n
                acc = mod(acc + U[k, i] * U[k, j], p)
            end
            if i == j
                acc == 1 || return false
            else
                acc == 0 || return false
            end
        end
    end
    return true
end

function _filter_orthogonal_mats(mats::Vector{Matrix{Int}}, p::Int)
    return [M for M in mats if is_orthogonal_mod_p(M, p)]
end

"""
    _extract_orthogonality_points(gb_data, n::Int, p::Int)
        -> Union{Nothing, Vector{Matrix{Int}}}

Best-effort extraction of F_p-rational O(n) points from a Gröbner ideal.
If Oscar point extraction APIs are unavailable, returns `nothing`.
"""
function _extract_orthogonality_points(gb_data, n::Int, p::Int)
    I = gb_data.ideal
    vars = gb_data.vars

    # Try several Oscar entry points (version-dependent).
    candidates = _collect_points_from_ideal(I)
    isempty(candidates) && return _extract_orthogonality_points_via_groebner(gb_data, n, p)

    mats = Matrix{Int}[]
    for pts in candidates
        if pts isa AbstractVector
            for pt in pts
                M = _coerce_point_to_matrix(pt, vars, n, p)
                M === nothing || push!(mats, M)
            end
        end
    end

    isempty(mats) && return _extract_orthogonality_points_via_groebner(gb_data, n, p)
    mats = _dedupe_matrices(mats)
    mats = _filter_orthogonal_mats(mats, p)
    isempty(mats) && return _extract_orthogonality_points_via_groebner(gb_data, n, p)
    return _sort_matrices_lex(mats)
end

function _fp_candidate_values_for_var(polys, varidx::Int,
                                      assigned::Dict{Int, Any}, F, p::Int)
    candidates = nothing
    constrained = false
    for f in polys
        coeffs = _fp_partial_univariate_coeffs(f, varidx, assigned, F)
        coeffs === nothing && continue
        any(pow -> pow > 0, keys(coeffs)) || continue
        roots = _fp_roots_from_coeffs(coeffs, F, p)
        roots === nothing && continue
        constrained = true
        isempty(roots) && return Any[]
        candidates = candidates === nothing ? roots : [x for x in candidates if any(==(x), roots)]
        isempty(candidates) && return Any[]
    end
    constrained && return candidates
    return Any[F(a) for a in 0:(p - 1)]
end

function _extract_orthogonality_points_via_groebner(gb_data, n::Int, p::Int)
    polys = gb_data.gb === nothing ? gb_data.equations : collect(gb_data.gb)
    F = base_ring(gb_data.ring)
    nvars_total = length(gb_data.vars)
    mats = Matrix{Int}[]

    function descend(varidx::Int, assigned::Dict{Int, Any})
        if varidx == 0
            vals = [assigned[i] for i in 1:nvars_total]
            all(iszero(_eval_fp_poly(f, vals)) for f in gb_data.equations) || return
            M = zeros(Int, n, n)
            t = 1
            for i in 1:n, j in 1:n
                M[i, j] = _fp_elem_to_int(vals[t], p)
                t += 1
            end
            push!(mats, M)
            return
        end
        for v in _fp_candidate_values_for_var(polys, varidx, assigned, F, p)
            next_assigned = copy(assigned)
            next_assigned[varidx] = v
            descend(varidx - 1, next_assigned)
        end
    end

    descend(nvars_total, Dict{Int, Any}())
    isempty(mats) && return nothing
    mats = _dedupe_matrices(mats)
    mats = _filter_orthogonal_mats(mats, p)
    isempty(mats) && return nothing
    return _sort_matrices_lex(mats)
end

function _collect_points_from_ideal(I)
    candidates = Any[]
    for fname in (:rational_points, :variety, :points)
        if isdefined(Oscar, fname)
            f = getproperty(Oscar, fname)
            try
                pts = f(I)
                pts === nothing || push!(candidates, pts)
            catch
            end
        end
    end
    return candidates
end

function _coerce_point_to_block_matrix(point, varsU, n_block::Int, p::Int)
    vals = Int[]
    for v in varsU
        entry = _point_entry(point, v)
        entry === nothing && return nothing
        local intval
        try
            intval = Int(entry)
        catch
            return nothing
        end
        push!(vals, mod(intval, p))
    end

    M = zeros(Int, n_block, n_block)
    t = 1
    for i in 1:n_block
        for j in 1:n_block
            M[i, j] = vals[t]
            t += 1
        end
    end
    return M
end

function _extract_U_blocks_from_verlinde_system(sys, p::Int)
    candidates = Any[]
    if haskey(sys, :gb) && sys.gb !== nothing
        append!(candidates, _collect_points_from_ideal(sys.gb))
    end
    append!(candidates, _collect_points_from_ideal(sys.ideal))
    isempty(candidates) && return nothing

    nvars_u = length(sys.varsU)
    n_block = isqrt(nvars_u)
    n_block * n_block == nvars_u || return nothing
    mats = Matrix{Int}[]
    for pts in candidates
        if pts isa AbstractVector
            for pt in pts
                M = _coerce_point_to_block_matrix(pt, sys.varsU, n_block, p)
                M === nothing || push!(mats, M)
            end
        end
    end
    isempty(mats) && return nothing
    mats = _dedupe_matrices(mats)
    mats = _filter_orthogonal_mats(mats, p)
    isempty(mats) && return nothing
    return _sort_matrices_lex(mats)
end

function _block_var_index(n_block::Int, i::Int, j::Int)
    return (i - 1) * n_block + j
end

function _cayley_param_index(n::Int, i::Int, j::Int)
    i < j || error("require i < j")
    k = 0
    for a in 1:(i-1)
        k += n - a
    end
    k += j - i
    return k
end

function _A_entry_from_params(varsA, n::Int, i::Int, j::Int)
    if i == j
        return zero(varsA[1])
    elseif i < j
        return varsA[_cayley_param_index(n, i, j)]
    else
        return -varsA[_cayley_param_index(n, j, i)]
    end
end

"""
    build_cayley_link_equations(varsA, varsU, n_block::Int)
        -> Vector

Build equations encoding `(I + A)U = I - A`, where `A` is the antisymmetric
matrix parameterized by `varsA` (Cayley parameters) and `U` is an `n_block×n_block`
matrix represented by `varsU` in row-major order.
"""
function build_cayley_link_equations(varsA, varsU, n_block::Int)
    length(varsA) == div(n_block * (n_block - 1), 2) || error("varsA length mismatch")
    length(varsU) == n_block * n_block || error("varsU length mismatch")

    eqs = Any[]
    U(i, j) = varsU[_block_var_index(n_block, i, j)]
    for i in 1:n_block
        for j in 1:n_block
            lhs = zero(varsU[1])
            for k in 1:n_block
                Ipk = (i == k) ? one(varsU[1]) : zero(varsU[1])
                lhs += (Ipk + _A_entry_from_params(varsA, n_block, i, k)) * U(k, j)
            end
            rhs = (i == j ? one(varsU[1]) : zero(varsU[1])) - _A_entry_from_params(varsA, n_block, i, j)
            push!(eqs, lhs - rhs)
        end
    end
    return eqs
end

function _u_full_entry(varsU, indices_deg::Vector{Int}, i::Int, j::Int, p::Int)
    bi = findfirst(==(i), indices_deg)
    bj = findfirst(==(j), indices_deg)
    if bi !== nothing && bj !== nothing
        n_block = length(indices_deg)
        return varsU[_block_var_index(n_block, bi, bj)]
    end
    return i == j ? one(varsU[1]) : zero(varsU[1])
end

"""
    build_verlinde_unit_equations(varsU, varsW, S_atomic, indices_deg, p, unit_idx)
        -> Vector

Build polynomial equations for a fixed candidate unit `unit_idx`:
- orthogonality on the active block U
- invertibility witnesses `w_m * S′[unit_idx,m] - 1 = 0`
- unit axiom `N_{unit_idx,i}^k = δ_{i,k}`

where `S′ = U^T S U` and `U` acts only on `indices_deg`.
"""
function build_verlinde_unit_equations(varsU, varsW,
                                       S_atomic::Matrix{Int},
                                       indices_deg::Vector{Int},
                                       p::Int, unit_idx::Int)
    r = size(S_atomic, 1)
    n_block = length(indices_deg)
    size(S_atomic, 2) == r || error("S_atomic must be square")
    length(varsU) == n_block * n_block || error("varsU length mismatch")
    length(varsW) == r || error("varsW length mismatch")

    eqs = Any[]

    # O(n_block): U^T U = I
    append!(eqs, build_orthogonality_equations(varsU, n_block, p))

    # Build symbolic S' = U^T S U
    Sprime = Matrix{Any}(undef, r, r)
    for i in 1:r
        for j in 1:r
            acc = zero(varsU[1])
            for a in 1:r
                Uai = _u_full_entry(varsU, indices_deg, a, i, p)
                # Split second contraction for better readability:
                for b in 1:r
                    Ubj = _u_full_entry(varsU, indices_deg, b, j, p)
                    acc += Uai * S_atomic[a, b] * Ubj
                end
            end
            Sprime[i, j] = acc
        end
    end

    # Witness inverse entries for S'[u,m]
    for m in 1:r
        push!(eqs, varsW[m] * Sprime[unit_idx, m] - one(varsU[1]))
    end

    # Unit axiom only: N_{u,i}^k = δ_{i,k}
    for i in 1:r
        for k in 1:r
            acc = zero(varsU[1])
            for m in 1:r
                acc += Sprime[unit_idx, m] * Sprime[i, m] * Sprime[k, m] * varsW[m]
            end
            push!(eqs, acc - ((i == k) ? one(varsU[1]) : zero(varsU[1])))
        end
    end

    return eqs
end

"""
    _build_verlinde_groebner_system(S_atomic, indices_deg, p, unit_idx)

Create Gröbner-ready polynomial data for fixed `unit_idx` over variables
`U_block` and inverse witnesses `w_m`.
"""
function _build_verlinde_groebner_system(S_atomic::Matrix{Int},
                                         indices_deg::Vector{Int},
                                         p::Int, unit_idx::Int)
    n_block = length(indices_deg)
    r = size(S_atomic, 1)
    F = GF(p)
    nA = div(n_block * (n_block - 1), 2)
    R, vars = polynomial_ring(F, nA + n_block * n_block + r, :z)
    varsA = vars[1:nA]
    varsU = vars[(nA + 1):(nA + n_block * n_block)]
    varsW = vars[(nA + n_block * n_block + 1):end]
    eqs = build_verlinde_unit_equations(varsU, varsW, S_atomic, indices_deg, p, unit_idx)
    append!(eqs, build_cayley_link_equations(varsA, varsU, n_block))
    I = ideal(R, eqs)
    G = nothing
    try
        G = groebner_basis(I)
    catch
    end
    return (ring = R, vars = vars, varsA = varsA, varsU = varsU, varsW = varsW,
            equations = eqs, ideal = I, gb = G)
end

function _warmup_verlinde_groebner!(S_atomic::Matrix{Int}, indices_deg::Vector{Int},
                                    p::Int, unit_idx::Int)
    key = (p, unit_idx, hash(S_atomic), Tuple(indices_deg))
    haskey(_VERLINDE_GB_CACHE, key) && return _VERLINDE_GB_CACHE[key]
    sys = _build_verlinde_groebner_system(S_atomic, indices_deg, p, unit_idx)
    _VERLINDE_GB_CACHE[key] = sys
    return sys
end

"""
    _extract_U_blocks_via_verlinde_groebner(S_atomic, indices_deg, p)
        -> Union{Nothing, Vector{Matrix{Int}}}

Try fixed-unit Gröbner systems for all possible unit indices and collect
all recovered `U_block` candidates.
"""
function _extract_U_blocks_via_verlinde_groebner(S_atomic::Matrix{Int},
                                                 indices_deg::Vector{Int},
                                                 p::Int;
                                                 max_units::Int = typemax(Int))
    r = size(S_atomic, 1)
    all_blocks = Matrix{Int}[]
    for unit_idx in 1:min(r, max_units)
        local sys
        try
            sys = _warmup_verlinde_groebner!(S_atomic, indices_deg, p, unit_idx)
        catch
            continue
        end
        blocks = _extract_U_blocks_from_verlinde_system(sys, p)
        if blocks !== nothing && !isempty(blocks)
            append!(all_blocks, blocks)
        end
    end

    isempty(all_blocks) && return nothing
    return _sort_matrices_lex(_dedupe_matrices(all_blocks))
end

"""
    enumerate_block_candidates(n_block::Int, p::Int, search_mode::Symbol)
        -> Vector{Matrix{Int}}

Choose the Phase-2 block-U search backend.

Supported modes:
- `:exhaustive`: Cayley + reflection sweep.
- `:groebner`: Gröbner basis plus finite-field point extraction.
"""
function enumerate_block_candidates(n_block::Int, p::Int, search_mode::Symbol)
    validate_search_mode(search_mode)
    if search_mode == :exhaustive
        return enumerate_o_n_Fp(n_block, p)
    elseif search_mode == :groebner
        return enumerate_o_n_Fp_groebner(n_block, p)
    end
end

"""
    validate_search_mode(search_mode::Symbol)

Validate Phase-2 search backend selector.
"""
function validate_search_mode(search_mode::Symbol)
    search_mode in (:exhaustive, :groebner) ||
        error("Unknown search_mode=$(search_mode). Expected :exhaustive or :groebner")
    return nothing
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

# ============================================================
#  Step 5: Verlinde integrality check
# ============================================================

"""
    signed_Fp(x::Int, p::Int) -> Int

Lift x ∈ [0, p) to the symmetric interval [-⌊p/2⌋, ⌊p/2⌋].
"""
function signed_Fp(x::Int, p::Int)
    return x <= p ÷ 2 ? x : x - p
end

"""
    verlinde_find_unit(S_Fp::Matrix{Int}, p::Int; threshold::Int = 5)
        -> Union{Nothing, Tuple{Int, Array{Int,3}}}

Given an F_p modular datum S, try each basis index u as the candidate
"unit object". For unit u to be valid:
1. S[u, m] ≠ 0 for all m (so 1/S[u,m] is defined).
2. N_{u, i}^{k} = δ_{i,k} (unit fusion is identity).
3. All N_{i,j}^k small non-negative integers (|n| ≤ threshold).

If found, returns (u, N) where N is the r×r×r tensor. Otherwise nothing.

The Verlinde formula (assuming S^2 = I, self-dual MTC):
    N_{ij}^k = Σ_m (S[i,m] · S[j,m] · S[k,m]) / S[u,m]
"""
function verlinde_find_unit(S_Fp::Matrix{Int}, p::Int; threshold::Int = 5)
    r = size(S_Fp, 1)

    for u in 1:r
        S_u_row = [S_Fp[u, m] for m in 1:r]
        any(x -> x == 0, S_u_row) && continue

        try
            S_u_inv = [invmod(x, p) for x in S_u_row]

            # First check unit axiom: N[u, i, k] = δ_{i,k}
            is_unit = true
            for i in 1:r
                for k in 1:r
                    val = 0
                    for m in 1:r
                        term = mod(S_Fp[u, m] * S_Fp[i, m], p)
                        term = mod(term * S_Fp[k, m], p)
                        term = mod(term * S_u_inv[m], p)
                        val = mod(val + term, p)
                    end
                    expected = (i == k) ? 1 : 0
                    if val != expected
                        is_unit = false
                        break
                    end
                end
                is_unit || break
            end
            is_unit || continue

            # Compute full tensor, check non-negative integer
            N = zeros(Int, r, r, r)
            ok = true
            for i in 1:r
                for j in 1:r
                    for k in 1:r
                        val = 0
                        for m in 1:r
                            term = mod(S_Fp[i, m] * S_Fp[j, m], p)
                            term = mod(term * S_Fp[k, m], p)
                            term = mod(term * S_u_inv[m], p)
                            val = mod(val + term, p)
                        end
                        s = signed_Fp(val, p)
                        if s < 0 || abs(s) > threshold
                            ok = false
                            break
                        end
                        N[i, j, k] = s
                    end
                    ok || break
                end
                ok || break
            end
            ok && return (u, N)
        catch
            continue
        end
    end
    return nothing
end

"""
    passes_unit_axiom(S_Fp::Matrix{Int}, p::Int, u::Int) -> Bool

Fast check of unit axiom only:
`N_{u,i}^k = δ_{i,k}` using the Verlinde formula over `F_p`.
"""
function passes_unit_axiom(S_Fp::Matrix{Int}, p::Int, u::Int)
    r = size(S_Fp, 1)
    S_u_row = [S_Fp[u, m] for m in 1:r]
    any(x -> x == 0, S_u_row) && return false

    local S_u_inv
    try
        S_u_inv = [invmod(x, p) for x in S_u_row]
    catch
        return false
    end

    for i in 1:r
        for k in 1:r
            val = 0
            for m in 1:r
                term = mod(S_Fp[u, m] * S_Fp[i, m], p)
                term = mod(term * S_Fp[k, m], p)
                term = mod(term * S_u_inv[m], p)
                val = mod(val + term, p)
            end
            expected = (i == k) ? 1 : 0
            val == expected || return false
        end
    end
    return true
end

"""
    solve_cayley_unit_filtered_blocks(S_atomic, indices_deg, p; max_units)
        -> Vector{Matrix{Int}}

Enumerate Cayley-parametrized O(n)(F_p) blocks, but retain only those whose
transformed `S'` satisfies the unit axiom for at least one candidate unit
index in `1:min(r, max_units)`.

This is a deterministic algebraic filter stage used by `:groebner` mode
when direct Gröbner point extraction is unavailable.
"""
function solve_cayley_unit_filtered_blocks(S_atomic::Matrix{Int},
                                           indices_deg::Vector{Int},
                                           p::Int;
                                           max_units::Int = typemax(Int))
    n_block = length(indices_deg)
    d = div(n_block * (n_block - 1), 2)
    if n_block == 1
        return [reshape([1], 1, 1)]
    end

    so_blocks = Matrix{Int}[]
    for params_tuple in Iterators.product(ntuple(_ -> 0:(p-1), d)...)
        U = cayley_so_n(collect(params_tuple), n_block, p)
        U === nothing || push!(so_blocks, U)
    end

    # Include reflection coset as in enumerate_o_n_Fp
    R = I_n(n_block)
    R[1, 1] = p - 1
    all_blocks = vcat(so_blocks, [matmul_mod(R, U, p) for U in so_blocks])

    r = size(S_atomic, 1)
    unit_max = min(r, max_units)
    kept = Matrix{Int}[]
    for U_block in all_blocks
        S_prime = apply_block_U(S_atomic, indices_deg, U_block, p)
        ok = false
        for u in 1:unit_max
            if passes_unit_axiom(S_prime, p, u)
                ok = true
                break
            end
        end
        ok && push!(kept, U_block)
    end
    return _sort_matrices_lex(_dedupe_matrices(kept))
end

# ============================================================
#  Step 6: Top-level single-prime driver
# ============================================================

"""
    MTCCandidate

Result of a successful block-U sweep at a single prime.

Fields:
- p: the prime used
- U_params: parameters of the block-U used (e.g. (u, v, det) for O(2))
- S_Fp: the S matrix mod p after block-U transformation
- T_Fp: the T eigenvalues mod p (unchanged by block-U)
- unit_index: which basis index is the unit object
- N: fusion tensor N[i][j][k] (signed integers)
- d: quantum dimensions d_i = S[unit, i] / S[unit, unit]
- D2: total quantum dimension squared
"""
struct MTCCandidate
    p::Int
    U_params::Any
    S_Fp::Matrix{Int}
    T_Fp::Vector{Int}
    unit_index::Int
    N::Array{Int, 3}
    d::Vector{Int}
    D2::Int
end

function Base.show(io::IO, c::MTCCandidate)
    d_signed = [signed_Fp(x, c.p) for x in c.d]
    # U_params can be a Matrix (general O(n)) or a Tuple (legacy O(2))
    params_str = if isa(c.U_params, AbstractMatrix)
        n = size(c.U_params, 1)
        "U=$(n)×$(n) block"
    elseif isa(c.U_params, AbstractVector) && all(x -> x isa NamedTuple && haskey(x, :U), c.U_params)
        dims = [size(x.U, 1) for x in c.U_params]
        "U blocks=$(join(dims, ","))"
    else
        "params=$(c.U_params)"
    end
    print(io, "MTCCandidate(p=$(c.p), unit=$(c.unit_index), ",
          "d=$d_signed, D²=$(signed_Fp(c.D2, c.p)), ",
          "$params_str)")
end

function _block_U_candidates_for_eigenspace(S_atomic::Matrix{Int},
                                            indices_deg::Vector{Int},
                                            p::Int;
                                            search_mode::Symbol,
                                            max_units_for_groebner::Int,
                                            groebner_allow_fallback::Bool,
                                            use_verlinde_warmup::Bool)
    n_block = length(indices_deg)
    U_blocks = Matrix{Int}[]

    # O(2)(F_p) is tiny and the explicit circle parametrisation is exact.
    # Avoid Gröbner warmup here; for p≈100 it is much slower than checking
    # the ~2p orthogonal blocks directly.
    n_block == 2 && return enumerate_o_n_Fp(n_block, p)

    if search_mode == :groebner
        if use_verlinde_warmup
            try
                gb_u_blocks = _extract_U_blocks_via_verlinde_groebner(S_atomic, indices_deg, p;
                                                                       max_units = max_units_for_groebner)
                if gb_u_blocks !== nothing && !isempty(gb_u_blocks)
                    U_blocks = gb_u_blocks
                end
            catch err
                @warn "Verlinde Gröbner warmup failed; trying orthogonality Gröbner extraction" exception = (err, catch_backtrace())
            end
        end

        if isempty(U_blocks)
            U_blocks = enumerate_o_n_Fp_groebner(n_block, p)
        end
    end

    # Explicit exhaustive fallback is opt-in for :groebner mode.
    if isempty(U_blocks)
        if search_mode == :groebner && groebner_allow_fallback
            U_blocks = enumerate_o_n_Fp(n_block, p)
        elseif search_mode != :groebner
            U_blocks = enumerate_block_candidates(n_block, p, search_mode)
        end
    end

    return U_blocks
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

function _candidate_from_fixed_S(S_atomic::Matrix{Int}, T_atomic::Vector{Int},
                                 p::Int, U_params;
                                 verlinde_threshold::Int)
    r = size(S_atomic, 1)
    result = verlinde_find_unit(S_atomic, p; threshold = verlinde_threshold)
    result === nothing && return MTCCandidate[]
    (u, N_tensor) = result
    Suu_inv = invmod(S_atomic[u, u], p)
    d = [mod(S_atomic[u, i] * Suu_inv, p) for i in 1:r]
    D2 = sum(mod(d[i] * d[i], p) for i in 1:r) % p
    return [MTCCandidate(p, U_params, copy(S_atomic), copy(T_atomic),
                         u, N_tensor, d, D2)]
end

"""
    find_mtcs_at_prime(catalog::Vector{AtomicIrrep}, stratum::Stratum,
                       p::Int; verlinde_threshold::Int = 3,
                       max_block_dim::Int = 3,
                       search_mode::Symbol = :groebner,
                       max_units_for_groebner::Int = typemax(Int),
                       groebner_allow_fallback::Bool = false,
                       precheck_unit_axiom::Bool = true)
        -> Vector{MTCCandidate}

Top-level Phase 2 driver at a single prime.

Steps:
1. Build block-diagonal atomic (S, T) from stratum + catalog.
2. Reduce to F_p (requires N | p-1 where N = catalog[].N).
3. Decompose T into eigenspaces; compute parameter_dim.
4. For each degenerate eigenspace of dimension n_θ ≥ 2, solve for
   block-U candidates using the selected backend, then take the Cartesian
   product across all degenerate eigenspaces.
5. For each block-U candidate, apply transformation to S, then check
   Verlinde integrality.
6. Return all MTC candidates found.

`max_block_dim` is a safety guard against very large block systems: if any
degenerate eigenspace has n_θ > max_block_dim, an error is raised.

`search_mode` selects the block-U candidate backend:
- `:groebner` (default): Gröbner basis plus finite-field point extraction
- `:exhaustive`: Cayley+reflection sweep

In `:groebner` mode, the driver performs best-effort Gröbner system
construction for fixed-unit variants and attempts direct U-block point
extraction. If that produces no candidates, it next tries Gröbner point
extraction for the orthogonality block itself.

`max_units_for_groebner` can cap how many unit indices are used for
fixed-unit Gröbner extraction (default: all).

`groebner_allow_fallback` controls whether `:groebner` mode is allowed to
fall back to the explicit `:exhaustive` Cayley/reflection backend when
solver extraction returns no candidates. It is disabled by default.

`precheck_unit_axiom` enables a fast unit-axiom prefilter before full
`verlinde_find_unit` tensor construction.
"""
function find_mtcs_at_prime(catalog::Vector{AtomicIrrep}, stratum::Stratum,
                            p::Int; verlinde_threshold::Int = 3,
                            max_block_dim::Int = 3,
                            search_mode::Symbol = :groebner,
                            max_units_for_groebner::Int = typemax(Int),
                            groebner_allow_fallback::Bool = false,
                            precheck_unit_axiom::Bool = true)
    validate_search_mode(search_mode)
    max_units_for_groebner >= 1 || error("max_units_for_groebner must be ≥ 1")
    # Common N
    N = catalog[first(keys(stratum.multiplicities))].N

    # Build block-diagonal atomic data
    S_K, T_K, K, blocks = build_block_diagonal(catalog, stratum)

    # Reduce to F_p
    zeta_N_Fp = find_zeta_in_Fp(N, p)
    S_atomic = [cyclotomic_to_Fp(S_K[i, j], zeta_N_Fp, p) for i in 1:size(S_K, 1), j in 1:size(S_K, 2)]
    T_atomic = [cyclotomic_to_Fp(t, zeta_N_Fp, p) for t in T_K]

    # T-eigenspace decomposition
    eigenspaces = t_eigenspace_decomposition(T_atomic, p)
    p_dim = parameter_dim(eigenspaces)

    degenerate = sort([(theta, indices) for (theta, indices) in eigenspaces if length(indices) >= 2];
                      by = x -> x[1])

    # Case 1: no degeneracy → atomic S is already an MTC candidate (modulo signs)
    if isempty(degenerate)
        return _candidate_from_fixed_S(S_atomic, T_atomic, p, :atomic;
                                       verlinde_threshold = verlinde_threshold)
    end

    if _orthogonal_block_product_fixes_S(S_atomic, degenerate, p)
        return _candidate_from_fixed_S(S_atomic, T_atomic, p, :block_invariant;
                                       verlinde_threshold = verlinde_threshold)
    end

    # Case 2+: one or more degenerate eigenspaces.
    block_candidate_sets = Vector{Vector{NamedTuple}}()
    use_verlinde_warmup = length(degenerate) == 1
    for (theta_deg, indices_deg) in degenerate
        n_block = length(indices_deg)
        n_block <= max_block_dim || error(
            "Degenerate eigenspace dim n=$n_block exceeds max_block_dim=$max_block_dim. " *
            "Increase max_block_dim only if the selected block-U backend can handle it.")

        U_blocks = _block_U_candidates_for_eigenspace(S_atomic, indices_deg, p;
                                                      search_mode = search_mode,
                                                      max_units_for_groebner = max_units_for_groebner,
                                                      groebner_allow_fallback = groebner_allow_fallback,
                                                      use_verlinde_warmup = use_verlinde_warmup)
        isempty(U_blocks) && return MTCCandidate[]
        push!(block_candidate_sets,
              [(theta = theta_deg, indices = indices_deg, U = U_block) for U_block in U_blocks])
    end

    candidates = MTCCandidate[]
    r = size(S_atomic, 1)
    for block_tuple in Iterators.product(block_candidate_sets...)
        block_choices = collect(block_tuple)
        S_prime = _apply_block_U_product(S_atomic, block_choices, p)
        if precheck_unit_axiom
            any_unit = any(u -> passes_unit_axiom(S_prime, p, u), 1:r)
            any_unit || continue
        end
        result = verlinde_find_unit(S_prime, p; threshold = verlinde_threshold)
        if result !== nothing
            (u_idx, N_tensor) = result
            Suu_inv = invmod(S_prime[u_idx, u_idx], p)
            d = [mod(S_prime[u_idx, i] * Suu_inv, p) for i in 1:r]
            D2 = sum(mod(d[m] * d[m], p) for m in 1:r) % p
            U_params = length(block_choices) == 1 ? block_choices[1].U : block_choices
            push!(candidates, MTCCandidate(
                p, U_params,
                copy(S_prime), copy(T_atomic),
                u_idx, N_tensor, d, D2))
        end
    end
    return candidates
end
