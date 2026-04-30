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
