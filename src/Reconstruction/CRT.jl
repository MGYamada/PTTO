"""
Phase 3: Multi-prime CRT reconstruction.

Given MTC candidates found at multiple primes via `find_mtcs_at_prime`
(Phase 2), this module:

1. Matches candidates across primes using fusion-tensor invariants
   (the integer fusion tensor N is prime-independent, so candidates
   with the same N_{ij}^k belong to the "same MTC")

2. Provides CRT and rational-reconstruction primitives shared by the
   cyclotomic modular-data lift.

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

Returns a Vector of Dicts, each mapping prime → MTCCandidate.  Multiple
finite-field candidates can have the same canonical fusion tensor, for example
from diagonal sign choices preserving `T`.  These variants are kept as
separate cross-prime groups and downstream cyclotomic reconstruction decides
which combinations are compatible.

Note: two candidates match by the canonicalized fusion tensor. Downstream
cyclotomic reconstruction and fresh-prime checks decide whether the finite
field modular data are compatible.
"""
function group_mtcs_by_fusion(results_by_prime::Dict{Int, Vector{MTCCandidate}};
                              debug_stable_key::Bool = false)
    buckets_by_key = Dict{String, Dict{Int, Vector{MTCCandidate}}}()
    key_to_raw = Dict{String, Set{String}}()
    for (p, candidates) in results_by_prime
        for c in candidates
            key = canonical_rule(c.N)
            prime_buckets = get!(buckets_by_key, key, Dict{Int, Vector{MTCCandidate}}())
            bucket = get!(prime_buckets, p, MTCCandidate[])
            sig = (c.unit_index, Tuple(vec(c.N)), Tuple(vec(c.S_Fp)), Tuple(c.T_Fp))
            if !any(existing -> (existing.unit_index, Tuple(vec(existing.N)),
                                 Tuple(vec(existing.S_Fp)), Tuple(existing.T_Fp)) == sig,
                    bucket)
                push!(bucket, c)
            end
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
    groups = Dict{Int, MTCCandidate}[]
    for key in sort!(collect(keys(buckets_by_key)))
        by_prime = buckets_by_key[key]
        primes = sort!(collect(keys(by_prime)))
        function emit(pi::Int, current::Dict{Int, MTCCandidate})
            if pi > length(primes)
                push!(groups, copy(current))
                return
            end
            p = primes[pi]
            for c in by_prime[p]
                current[p] = c
                emit(pi + 1, current)
            end
            delete!(current, p)
        end
        emit(1, Dict{Int, MTCCandidate}())
    end
    return groups
end

# ============================================================
#  Step 3: Rational reconstruction helper
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
