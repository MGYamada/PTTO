"""
Phase 3: Multi-prime CRT reconstruction.

Given MTC candidates found at multiple primes via `find_mtcs_at_prime`
(Phase 2), this module:

1. Matches candidates across primes using fusion-tensor invariants
   (the integer fusion tensor N is prime-independent, so candidates
   with the same N_{ij}^k belong to the "same MTC")

2. Reconstructs S'-matrix entries as integers in Z via chinese
   remainder theorem (CRT)

3. Attempts to recognize each entry as an element of a small algebraic
   extension (e.g. Z[√d] for small d) given a basis guess

4. Cross-validates at primes not used in reconstruction

Usage:

    # At each prime, run find_mtcs_at_prime
    results_by_prime = Dict(p => find_mtcs_at_prime(...) for p in primes)

    # Group across primes by matching fusion tensors
    groups = group_mtcs_by_fusion(results_by_prime)

    # For each group, attempt CRT reconstruction
    for group in groups
        recon = reconstruct_mtc(group, used_primes)
        validate_at(recon, fresh_primes)
    end
"""

# ============================================================
#  Step 1: CRT primitives
# ============================================================

"""
    crt2(r1::Int, m1::Int, r2::Int, m2::Int) -> (Int, Int)

Solve x ≡ r1 (mod m1), x ≡ r2 (mod m2) via extended Euclidean.
Returns (x, m1*m2) with 0 ≤ x < m1*m2.

Requires gcd(m1, m2) = 1 (coprime moduli).
"""
function crt2(r1::Integer, m1::Integer, r2::Integer, m2::Integer)
    # x = r1 + m1 * k, find k s.t. r1 + m1*k ≡ r2 (mod m2)
    # m1 * k ≡ r2 - r1 (mod m2)
    # k ≡ (r2 - r1) * m1^{-1} (mod m2)
    m1_inv = invmod(m1 % m2, m2)
    k = mod((r2 - r1) * m1_inv, m2)
    M = m1 * m2
    return (mod(r1 + m1 * k, M), M)
end

"""
    acmg_crt(residues::Vector{<:Integer}, moduli::Vector{<:Integer})
        -> (Int, Int)

CRT across multiple moduli. Pairwise fold via `crt2`.
Returns (x, M) where M = prod(moduli) and x is the unique residue
in [0, M).

Named `acmg_crt` to avoid name collision with `Oscar.crt` (also exported
by AbstractAlgebra, Nemo, Hecke, Singular).
"""
function acmg_crt(residues::Vector{<:Integer}, moduli::Vector{<:Integer})
    length(residues) == length(moduli) || error("length mismatch")
    length(residues) >= 1 || error("empty input")

    x = BigInt(residues[1])
    M = BigInt(moduli[1])
    for i in 2:length(residues)
        (x, M) = crt2(x, M, BigInt(residues[i]), BigInt(moduli[i]))
    end
    return (x, M)
end

"""
    rational_reconstruct(c::Integer, M::Integer; bound::Integer = isqrt(M ÷ 2))
        -> Union{Nothing, Tuple{BigInt, BigInt}}

Given c mod M, attempt to find (a, b) with a/b ≡ c (mod M) and
|a|, b ≤ bound (default ⌊√(M/2)⌋).

Returns (a, b) or nothing. Uses the extended Euclidean method
(Wang–Monagan form).
"""
function rational_reconstruct(c::Integer, M::Integer;
                              bound::Union{Integer,Nothing} = nothing)
    if bound === nothing
        bound = isqrt(M ÷ 2)
    end
    # Extended Euclidean on (M, c)
    r0, r1 = BigInt(M), BigInt(mod(c, M))
    s0, s1 = BigInt(0), BigInt(1)
    while r1 != 0
        if r1 <= bound && abs(s1) <= bound
            return (BigInt(r1), BigInt(s1))
        end
        q = r0 ÷ r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    end
    return nothing
end

# ============================================================
#  Step 2: Fusion-tensor matching across primes
# ============================================================

"""
    fusion_signature(c::MTCCandidate) -> Array{Int, 3}

The fusion tensor of an MTC candidate — prime-independent integer values
that serve as a canonical invariant for matching across primes.
"""
fusion_signature(c::MTCCandidate) = c.N

function _label_signature(N::Array{Int, 3}, i::Int)
    r = size(N, 1)
    s1 = 0
    s2 = 0
    s3 = 0
    q1 = 0
    q2 = 0
    q3 = 0
    d1 = 0
    d2 = 0
    d3 = 0
    @inbounds for j in 1:r, k in 1:r
        a = N[i, j, k]
        b = N[j, i, k]
        c = N[j, k, i]
        s1 += a
        s2 += b
        s3 += c
        q1 += a * a
        q2 += b * b
        q3 += c * c
    end
    @inbounds for j in 1:r
        d1 += N[i, j, j]
        d2 += N[j, i, j]
        d3 += N[j, j, i]
    end
    return (s1, s2, s3, q1, q2, q3, d1, d2, d3)
end

function _perm_lex_less(N::Array{Int, 3}, perm::Vector{Int}, best::Vector{Int})
    idx = 1
    r = size(N, 1)
    @inbounds for k in 1:r, j in 1:r, i in 1:r
        v = N[perm[i], perm[j], perm[k]]
        b = best[idx]
        if v < b
            return true
        elseif v > b
            return false
        end
        idx += 1
    end
    return false
end

function _flatten_permuted(N::Array{Int, 3}, perm::Vector{Int})
    r = size(N, 1)
    out = Vector{Int}(undef, r * r * r)
    idx = 1
    @inbounds for k in 1:r, j in 1:r, i in 1:r
        out[idx] = N[perm[i], perm[j], perm[k]]
        idx += 1
    end
    return out
end

function _for_each_permutation!(f::F, v::Vector{Int}) where {F}
    function rec(start::Int)
        if start > length(v)
            f(v)
            return
        end
        for i in start:length(v)
            v[start], v[i] = v[i], v[start]
            rec(start + 1)
            v[start], v[i] = v[i], v[start]
        end
    end
    rec(1)
    return nothing
end

"""
    canonical_rule(rule::AbstractArray{<:Integer, 3}) -> String

Return a deterministic, order-insensitive key for a fusion rule tensor.
All simple-object labels are treated as unlabeled (full permutation
canonicalization), integer types are normalized, and a blockwise
lexicographic minimization is used as payload:

- labels are first partitioned by permutation-invariant signatures;
- only labels inside the same signature block are permuted;
- we minimize lexicographically over the product of per-block permutations.

This reduces search from `r!` to `∏_b |b|!` without changing key stability.
"""
function canonical_rule(rule::AbstractArray{<:Integer, 3})
    size(rule, 1) == size(rule, 2) == size(rule, 3) ||
        error("canonical_rule expects a rank-r cubic tensor")
    N_int = Array{Int, 3}(rule)
    r = size(N_int, 1)

    labels = collect(1:r)
    sigs = [_label_signature(N_int, i) for i in labels]
    order = sortperm(labels; by = i -> sigs[i])
    sorted_labels = labels[order]
    sorted_sigs = sigs[order]

    blocks = Vector{Vector{Int}}()
    i = 1
    while i <= r
        j = i
        while j < r && sorted_sigs[j + 1] == sorted_sigs[i]
            j += 1
        end
        push!(blocks, sorted_labels[i:j])
        i = j + 1
    end

    best = nothing
    current = Vector{Int}(undef, r)

    function dfs_block(bi::Int, offset::Int)
        if bi > length(blocks)
            if best === nothing || _perm_lex_less(N_int, current, best)
                best = _flatten_permuted(N_int, current)
            end
            return
        end

        block = copy(blocks[bi])
        _for_each_permutation!(block) do p
            @inbounds for t in eachindex(p)
                current[offset + t - 1] = p[t]
            end
            dfs_block(bi + 1, offset + length(p))
        end
    end
    dfs_block(1, 1)

    payload = join(best, ",")
    return "r=$(r)|$payload"
end

"""
    group_mtcs_by_fusion(results_by_prime::Dict{Int, Vector{MTCCandidate}})
        -> Vector{Dict{Int, MTCCandidate}}

Group MTC candidates from multiple primes into groups, where each group
corresponds to "the same MTC" (as identified by matching fusion tensors).

Returns a Vector of Dicts, each mapping prime → MTCCandidate.

Note: two candidates match only if their full N-tensor and unit_index
agree. For Galois-conjugate MTC pairs (which have identical fusion
tensors but differ on S-matrix entries involving √d irrationals), both
conjugates will appear in the same group here — downstream consumers
must pick a consistent Galois sector across primes (see
`group_mtcs_galois_aware`).
"""
function group_mtcs_by_fusion(results_by_prime::Dict{Int, Vector{MTCCandidate}};
                              debug_stable_key::Bool = false)
    groups_by_key = Dict{String, Dict{Int, MTCCandidate}}()
    key_to_raw = Dict{String, Set{String}}()
    for (p, candidates) in results_by_prime
        for c in candidates
            key = canonical_rule(c.N)
            group = get!(groups_by_key, key, Dict{Int, MTCCandidate}())
            haskey(group, p) || (group[p] = c)
            if debug_stable_key
                raw = "unit=$(c.unit_index)|" * join(vec(c.N), ",")
                push!(get!(key_to_raw, key, Set{String}()), raw)
            end
        end
    end
    if debug_stable_key
        for key in sort!(collect(keys(key_to_raw)))
            println("[stable_key] key=$key")
            for raw in sort!(collect(key_to_raw[key]))
                println("  raw_rule=$raw")
            end
        end
    end
    groups = collect(values(groups_by_key))
    return groups
end

"""
    group_mtcs_galois_aware(results_by_prime::Dict{Int, Vector{MTCCandidate}},
                             anchor_prime::Int;
                             quadratic_d::Int = 3,
                             reconstruction_bound::Int = 50,
                             sqrtd_fn = (d, p) -> compute_sqrt3_cyclotomic_mod_p(p))
        -> Vector{Dict{Int, MTCCandidate}}

Group MTC candidates across primes such that each group has a CONSISTENT
Galois sector — i.e., all (S_Fp, T_Fp) values at all primes lift to the
same element of Z[ζ_N] (not just up to Galois conjugation).

Method: for each candidate at `anchor_prime`, seed a new group; then at
each other prime, pick the candidate whose S_Fp (together with T_Fp)
lifts consistently — detected by attempting a 2-prime reconstruction of
`scale_factor · √d · S_Fp` in Z[√d].

The `sqrtd_fn(d, p)` callback supplies a specific cyclotomic branch of
`√d mod p` — it MUST match the `d` / branch convention used by
downstream `reconstruct_S_matrix` or the final CRT reconstruction will
disagree with the Galois alignment performed here. The caller should
pass the same `sqrtd_fn` as `classify_from_group` / `reconstruct_S_matrix`.
"""
function group_mtcs_galois_aware(results_by_prime::Dict{Int, Vector{MTCCandidate}},
                                  anchor_prime::Int;
                                  quadratic_d::Int = 3,
                                  reconstruction_bound::Int = 50,
                                  sqrtd_fn = (d, p) -> compute_sqrt3_cyclotomic_mod_p(p),
                                  branch_sign_getter = nothing,
                                  branch_sign_setter = nothing,
                                  kwargs...)
    isempty(kwargs) || error("unsupported keyword arguments: $(collect(keys(kwargs)))")
    haskey(results_by_prime, anchor_prime) || error("anchor_prime $anchor_prime not in results")
    anchor_cands = results_by_prime[anchor_prime]

    # For each anchor candidate, seed a group. For each other prime's
    # candidate, accept it if `2·√quadratic_d · S_Fp` from {anchor, p} can be
    # reconstructed as Z[√quadratic_d] with small bounded coefficients. This
    # is the same test `verify_reconstruction` does, just truncated to
    # two primes.

    groups = Vector{Dict{Int, MTCCandidate}}()

    for anchor_c in anchor_cands
        anchor_key = canonical_rule(anchor_c.N)
        group = Dict{Int, MTCCandidate}(anchor_prime => anchor_c)

        for (p, cands) in results_by_prime
            p == anchor_prime && continue

            best_match = nothing
            best_sign = nothing
            best_score = typemax(Int)
            for c in cands
                # Different fusion tensor → skip
                # unit_index may legitimately differ across primes for
                # equivalent candidates; use fusion-tensor equality as the
                # compatibility guard here.
                canonical_rule(c.N) == anchor_key || continue

                # Check: can we reconstruct 2·√d·S as Z[√d] from just {anchor, p}?
                trial_signs = [nothing]
                if branch_sign_getter !== nothing && branch_sign_setter !== nothing
                    cur_sign = branch_sign_getter(p)
                    trial_signs = cur_sign == 1 ? [1, -1] : [-1, 1]
                end

                for sgn in trial_signs
                    if branch_sign_setter !== nothing && sgn !== nothing
                        branch_sign_setter(p, sgn)
                    end
                    try
                        s_anchor = sqrtd_fn(quadratic_d, anchor_prime)
                        s_p = sqrtd_fn(quadratic_d, p)
                        two_s_anchor = mod(2 * s_anchor, anchor_prime)
                        two_s_p = mod(2 * s_p, p)
                        nrow = size(anchor_c.S_Fp, 1)
                        matrix_by_prime = Dict{Int, Matrix{Int}}()
                        matrix_by_prime[anchor_prime] = [
                            mod(two_s_anchor * anchor_c.S_Fp[i, j], anchor_prime)
                            for i in 1:nrow, j in 1:nrow]
                        matrix_by_prime[p] = [
                            mod(two_s_p * c.S_Fp[i, j], p)
                            for i in 1:nrow, j in 1:nrow]
                        recon = reconstruct_matrix_in_Z_sqrt_d(
                            matrix_by_prime, quadratic_d;
                            bound = reconstruction_bound,
                            sqrtd_fn = sqrtd_fn)
                        # Success: score by coefficient size to prefer the
                        # most "integral-small" reconstruction.
                        score = 0
                        for entry in recon
                            score += abs(entry[1]) + abs(entry[2])
                        end
                        if score < best_score ||
                           (score == best_score && best_sign !== nothing && sgn !== best_sign)
                            best_match = c
                            best_sign = sgn
                            best_score = score
                        elseif score == best_score && best_sign === nothing
                            best_match = c
                            best_sign = sgn
                            best_score = score
                        end
                    catch
                        # Not compatible with this candidate/sign; try next
                        continue
                    end
                end
                if best_match !== nothing
                    break
                end
            end

            if best_match !== nothing
                group[p] = best_match
                if branch_sign_setter !== nothing && best_sign !== nothing
                    branch_sign_setter(p, best_sign)
                end
            end
        end

        push!(groups, group)
    end

    return groups
end

"""
    build_sqrtd_selector(quadratic_d::Int, primes::Vector{Int}, anchor_prime::Int; verbose::Bool = false)
        -> NamedTuple

Build a per-`d` square-root selector that is consistent across primes.

- `d ∈ {2,3,5}`: use cyclotomic branch with per-prime sign cache
  (`mode = :cyclotomic_signed`).
- otherwise: use anchored mode (`mode = :anchored`) with a sign cache
  relative to `anchor_prime`. The returned `sqrtd_fn` is a transform
  layer on top of raw Tonelli-style roots, so callers can align branches
  by updating `branch_sign_setter`.
"""
function build_sqrtd_selector(quadratic_d::Int, primes::Vector{Int}, anchor_prime::Int;
                              verbose::Bool = false)
    if quadratic_d in (2, 3, 5)
        sign_by_prime = Dict{Int, Int}(anchor_prime => 1)
        cyclotomic_base = function (p::Int)
            if quadratic_d == 2
                return compute_sqrt2_cyclotomic_mod_p(p)
            elseif quadratic_d == 3
                return compute_sqrt3_cyclotomic_mod_p(p)
            else
                return compute_sqrt5_cyclotomic_mod_p(p)
            end
        end
        branch_sign_getter = function (p::Int)
            return get(sign_by_prime, p, 1)
        end
        branch_sign_setter = function (p::Int, sgn::Int)
            sign_by_prime[p] = sgn >= 0 ? 1 : -1
            return nothing
        end
        sqrtd_fn = function (d::Int, p::Int)
            d == quadratic_d || error("cyclotomic selector was built for d=$quadratic_d, got d=$d")
            base = cyclotomic_base(p)
            sgn = branch_sign_getter(p)
            return sgn == 1 ? base : mod(-base, p)
        end
        verbose && println("  √d selector: d=$quadratic_d, branch=cyclotomic-signed (anchor prime=$anchor_prime)")
        return (sqrtd_fn = sqrtd_fn,
                mode = :cyclotomic_signed,
                anchor_prime = anchor_prime,
                branch_sign_getter = branch_sign_getter,
                branch_sign_setter = branch_sign_setter)
    end

    sign_by_prime = Dict{Int, Int}(anchor_prime => 1)
    raw_cache = Dict{Int, Int}()
    raw_root = function (p::Int)
        if !haskey(raw_cache, p)
            s = compute_sqrt_d_mod_p(quadratic_d, p)
            s === nothing && error("$quadratic_d is not a QR mod $p")
            raw_cache[p] = s
        end
        return raw_cache[p]
    end
    branch_sign_getter = function (p::Int)
        return get(sign_by_prime, p, 1)
    end
    branch_sign_setter = function (p::Int, sgn::Int)
        sign_by_prime[p] = sgn >= 0 ? 1 : -1
        return nothing
    end
    sqrtd_fn = function (d::Int, p::Int)
        d == quadratic_d || error("anchored selector was built for d=$quadratic_d, got d=$d")
        s = raw_root(p)
        sgn = branch_sign_getter(p)
        return sgn == 1 ? s : mod(-s, p)
    end

    verbose && println("  √d selector: d=$quadratic_d, branch=anchored (anchor prime=$anchor_prime)")
    return (sqrtd_fn = sqrtd_fn,
            mode = :anchored,
            anchor_prime = anchor_prime,
            branch_sign_getter = branch_sign_getter,
            branch_sign_setter = branch_sign_setter)
end

# ============================================================
#  Step 3: Algebraic reconstruction of matrix entries
# ============================================================

"""
    reconstruct_rational(values_by_prime::Dict{Int, Int})
        -> Union{Nothing, Tuple{BigInt, BigInt}}

Try to reconstruct a rational a/b consistent with `values_by_prime[p] ≡ a/b (mod p)`
for all given primes.

Returns (a, b) or nothing.
"""
function reconstruct_rational(values_by_prime::Dict{Int, Int})
    primes = collect(keys(values_by_prime))
    residues = [values_by_prime[p] for p in primes]
    (c, M) = acmg_crt(residues, primes)
    return rational_reconstruct(c, M)
end

"""
    reconstruct_in_Z_sqrt_d(values_by_prime::Dict{Int, Int},
                             d::Int,
                             sqrtd_by_prime::Dict{Int, Int})
        -> Union{Nothing, Tuple{BigInt, BigInt}}

Try to reconstruct x ∈ Z[√d] in the form x = a + b·√d where a, b ∈ Z,
given:
- values_by_prime[p] = x mod p
- sqrtd_by_prime[p] = √d mod p (chosen branch, consistent across primes)

Returns (a, b) or nothing. Searches (a, b) with small absolute value.

Algorithm: search small (a, b) by brute force up to some bound,
checking consistency at all primes.
"""
function reconstruct_in_Z_sqrt_d(values_by_prime::Dict{Int, Int},
                                  d::Int,
                                  sqrtd_by_prime::Dict{Int, Int};
                                  bound::Int = 5)
    primes = collect(keys(values_by_prime))
    for a in -bound:bound
        for b in -bound:bound
            all_ok = true
            for p in primes
                s3 = sqrtd_by_prime[p]
                expected = mod(a + b * s3, p)
                if expected != values_by_prime[p]
                    all_ok = false
                    break
                end
            end
            if all_ok
                return (BigInt(a), BigInt(b))
            end
        end
    end
    return nothing
end

"""
    compute_sqrt_d_mod_p(d::Int, p::Int) -> Union{Nothing, Int}

Compute one √d mod p (a specific branch chosen by returning the smaller
representative). Returns nothing if d is not a QR mod p.

Note: this chooses a branch by "smaller representative", which is NOT
Galois-consistent across primes. For multi-prime CRT reconstruction,
use `compute_sqrt3_cyclotomic_mod_p(p, zeta24)` or similar functions
that use a cyclotomic formula for consistent branch selection.
"""
function compute_sqrt_d_mod_p(d::Int, p::Int)
    d_mod = mod(d, p)
    if d_mod == 0
        return 0
    end
    # Tonelli-Shanks would be proper; for small p brute is fine
    for x in 1:(p - 1)
        if mod(x * x, p) == d_mod
            return min(x, p - x)  # smaller branch
        end
    end
    return nothing
end

"""
    compute_sqrt3_cyclotomic_mod_p(p::Int) -> Int

Compute √3 mod p using the cyclotomic identity √3 = ζ₂₄² + ζ₂₄⁻².
This gives a branch consistent with the choice of ζ₂₄ returned by
`find_zeta_in_Fp(24, p)`, so it is Galois-consistent across primes
(all primes see √3 as coming from the same Galois-orbit element).

Requires 24 | p - 1.
"""
function compute_sqrt3_cyclotomic_mod_p(p::Int)
    zeta24 = find_zeta_in_Fp(24, p)
    # √3 = ζ₂₄² + ζ₂₄⁻² = ζ₂₄² + ζ₂₄²²
    return mod(powermod(zeta24, 2, p) + powermod(zeta24, 22, p), p)
end

"""
    compute_sqrt2_cyclotomic_mod_p(p::Int) -> Int

Compute √2 mod p using the cyclotomic identity √2 = ζ₈ + ζ₈⁻¹ = ζ₂₄³ + ζ₂₄²¹.
Galois-consistent across primes with 24 | p - 1.
"""
function compute_sqrt2_cyclotomic_mod_p(p::Int)
    zeta24 = find_zeta_in_Fp(24, p)
    return mod(powermod(zeta24, 3, p) + powermod(zeta24, 21, p), p)
end

"""
    compute_sqrt5_cyclotomic_mod_p(p::Int) -> Int

Compute √5 mod p using the quadratic Gauss-sum identity
`√5 = 1 + 2·(ζ_5 + ζ_5^4)`. Galois-consistent across primes satisfying
`5 | p - 1` (which is required for ζ_5 ∈ F_p).
"""
function compute_sqrt5_cyclotomic_mod_p(p::Int)
    zeta5 = find_zeta_in_Fp(5, p)
    # √5 = 1 + 2·(ζ_5 + ζ_5^4)
    return mod(1 + 2 * (zeta5 + powermod(zeta5, 4, p)), p)
end

"""
    reconstruct_matrix_in_Z_sqrt_d(matrix_by_prime::Dict{Int, Matrix{Int}},
                                    d::Int; bound::Int = 5,
                                    sqrtd_fn = compute_sqrt_d_mod_p)
        -> Matrix{Tuple{Int, Int}}

Reconstruct a matrix of Z[√d] entries from F_p reductions across primes.
Each entry is returned as a pair (a, b) meaning a + b·√d.

`sqrtd_fn(d, p)` should return a specific branch of √d mod p. The same
branch must be chosen consistently across all input primes for CRT to
succeed. Default `compute_sqrt_d_mod_p` chooses the smaller representative
(which may NOT be Galois-consistent). For d = 3 with primes dividing
24, pass `(d, p) -> compute_sqrt3_cyclotomic_mod_p(p)` for consistency.

Requires d to be a QR at all input primes.
"""
function reconstruct_matrix_in_Z_sqrt_d(matrix_by_prime::Dict{Int, Matrix{Int}},
                                         d::Int; bound::Int = 5,
                                         sqrtd_fn = compute_sqrt_d_mod_p)
    primes = collect(keys(matrix_by_prime))
    # Precompute √d mod p for each prime
    sqrtd_by_prime = Dict{Int, Int}()
    for p in primes
        s = sqrtd_fn(d, p)
        s === nothing && error("$d is not a QR mod $p")
        sqrtd_by_prime[p] = s
    end

    # Get matrix dimensions from any prime
    first_mat = first(values(matrix_by_prime))
    nrow, ncol = size(first_mat)

    result = Matrix{Tuple{Int, Int}}(undef, nrow, ncol)
    for i in 1:nrow
        for j in 1:ncol
            values_ij = Dict(p => matrix_by_prime[p][i, j] for p in primes)
            recon = reconstruct_in_Z_sqrt_d(values_ij, d, sqrtd_by_prime; bound = bound)
            if recon === nothing
                error("Failed to reconstruct entry ($i, $j) in Z[√$d] with bound=$bound. " *
                      "This may indicate: (a) inconsistent sqrt branch across primes " *
                      "(try passing a cyclotomic-consistent sqrtd_fn), " *
                      "(b) entry is not in Z[√$d] (try larger bound or different d), " *
                      "(c) the MTC candidate matching was wrong.")
            end
            result[i, j] = (Int(recon[1]), Int(recon[2]))
        end
    end
    return result
end

# ============================================================
#  Step 4: End-to-end reconstruction of an MTC across primes
# ============================================================

"""
    reconstruct_S_matrix(group::Dict{Int, MTCCandidate}, D::Int;
                         bound::Int = 5)
        -> Matrix{Tuple{Int, Int}}

Given a group of matching MTCCandidates indexed by prime, reconstruct
the (2D · S')-matrix as a matrix of Z[√(D² * something)] entries.

For D² = 12 (SU(2)_4), we use d = 3 and reconstruct 2√3 · S' in Z[√3].
This is the standard normalization used in the Python pipeline.

The factor-of-2D scaling is chosen so S = (1/(2D)) · M has M with integer
(or Z[√d]) entries.

Actually more cleanly: for D² = 12 we scale by 2√3 (not 2D = 4√3), because
the Kac-Peterson S has an extra 1/2 factor in its natural form. The result
M has entries in Z[√3].

This function takes the parameter `quadratic_d`: the d such that the scaling
is √(quadratic_d). For SU(2)_4 (D² = 12), quadratic_d = 3.
"""
function reconstruct_S_matrix(group::Dict{Int, MTCCandidate};
                              quadratic_d::Int = 3,
                              bound::Int = 5,
                              sqrtd_fn = compute_sqrt_d_mod_p,
                              kwargs...)
    isempty(kwargs) || error("unsupported keyword arguments: $(collect(keys(kwargs)))")
    primes = collect(keys(group))

    # Build matrix_by_prime: S_scaled[p][i, j] = (2 · √quadratic_d · S'[i, j]) mod p
    matrix_by_prime = Dict{Int, Matrix{Int}}()
    for p in primes
        s = sqrtd_fn(quadratic_d, p)
        s === nothing && error("$quadratic_d is not a QR mod $p")
        two_s = mod(2 * s, p)
        S_p = group[p].S_Fp
        M_p = [mod(two_s * S_p[i, j], p) for i in 1:size(S_p, 1), j in 1:size(S_p, 2)]
        matrix_by_prime[p] = M_p
    end

    return reconstruct_matrix_in_Z_sqrt_d(matrix_by_prime, quadratic_d; bound = bound, sqrtd_fn = sqrtd_fn)
end

"""
    verify_reconstruction(recon::Matrix{Tuple{Int, Int}},
                          candidate::MTCCandidate,
                          d::Int;
                          scale::Int = 2)
        -> Bool

Verify that the reconstructed matrix (with entries (a, b) meaning a + b·√d)
reproduces `candidate.S_Fp` when reduced mod candidate.p.

The assumption is that the reconstructed matrix represents `scale · √d · S`,
so we check (scale · √d · S_Fp) ≡ recon (mod p) pointwise.
"""
function verify_reconstruction(recon::Matrix{Tuple{Int, Int}},
                                candidate::MTCCandidate,
                                d::Int;
                                scale::Int = 2,
                                sqrtd_fn = compute_sqrt_d_mod_p)
    p = candidate.p
    s = sqrtd_fn(d, p)
    s === nothing && return false
    scale_s = mod(scale * s, p)

    nrow, ncol = size(recon)
    for i in 1:nrow
        for j in 1:ncol
            (a, b) = recon[i, j]
            expected = mod(a + b * s, p)
            actual = mod(scale_s * candidate.S_Fp[i, j], p)
            if expected != actual
                return false
            end
        end
    end
    return true
end

"""
    describe_matrix(M::Matrix{Tuple{Int, Int}}, d::Int) -> String

Human-readable rendering of a matrix whose entries are (a, b) = a + b√d.
"""
function describe_matrix(M::Matrix{Tuple{Int, Int}}, d::Int)
    nrow, ncol = size(M)
    sqrtd_sym = "√$d"

    # Format each entry
    fmt(pair) = begin
        (a, b) = pair
        if a == 0 && b == 0
            "0"
        elseif b == 0
            "$a"
        elseif a == 0
            if b == 1
                sqrtd_sym
            elseif b == -1
                "-" * sqrtd_sym
            else
                "$b" * sqrtd_sym
            end
        else
            if b == 1
                "$a+" * sqrtd_sym
            elseif b == -1
                "$a-" * sqrtd_sym
            elseif b > 0
                "$a+$b" * sqrtd_sym
            else
                "$a$b" * sqrtd_sym
            end
        end
    end

    lines = String[]
    for i in 1:nrow
        row_strs = [fmt(M[i, j]) for j in 1:ncol]
        # pad each entry to same width
        w = maximum(length, row_strs)
        padded = [lpad(s, w) for s in row_strs]
        push!(lines, "  " * join(padded, "  "))
    end
    return join(lines, "\n")
end
