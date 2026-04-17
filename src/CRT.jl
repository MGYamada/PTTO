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
    crt(residues::Vector{<:Integer}, moduli::Vector{<:Integer})
        -> (Int, Int)

CRT across multiple moduli. Pairwise fold via `crt2`.
Returns (x, M) where M = prod(moduli) and x is the unique residue
in [0, M).
"""
function crt(residues::Vector{<:Integer}, moduli::Vector{<:Integer})
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

"""
    group_mtcs_by_fusion(results_by_prime::Dict{Int, Vector{MTCCandidate}})
        -> Vector{Dict{Int, MTCCandidate}}

Group MTC candidates from multiple primes into groups, where each group
corresponds to "the same MTC" (as identified by matching fusion tensors).

Returns a Vector of Dicts, each mapping prime → MTCCandidate. A group
should contain one candidate per prime for successful CRT.

Note: two candidates match only if their full N-tensor and unit_index
agree. This requires the F_p results at each prime to produce the unit
at the same basis index, which depends on candidate ordering.
"""
function group_mtcs_by_fusion(results_by_prime::Dict{Int, Vector{MTCCandidate}})
    groups = Vector{Dict{Int, MTCCandidate}}()
    # For each prime, for each candidate, find matching group or create new
    for (p, candidates) in results_by_prime
        for c in candidates
            matched = false
            for g in groups
                # Pick any representative and compare N + unit
                rep = first(values(g))
                if rep.N == c.N && rep.unit_index == c.unit_index
                    # Same MTC — add this prime's candidate
                    # (if this prime already in group, skip — first wins)
                    if !haskey(g, p)
                        g[p] = c
                    end
                    matched = true
                    break
                end
            end
            if !matched
                new_group = Dict{Int, MTCCandidate}(p => c)
                push!(groups, new_group)
            end
        end
    end
    return groups
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
    (c, M) = crt(residues, primes)
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
    reconstruct_matrix_in_Z_sqrt_d(matrix_by_prime::Dict{Int, Matrix{Int}},
                                    d::Int; bound::Int = 5)
        -> Matrix{Tuple{Int, Int}}

Reconstruct a matrix of Z[√d] entries from F_p reductions across primes.
Each entry is returned as a pair (a, b) meaning a + b·√d.

Requires d to be a QR at all input primes.
"""
function reconstruct_matrix_in_Z_sqrt_d(matrix_by_prime::Dict{Int, Matrix{Int}},
                                         d::Int; bound::Int = 5)
    primes = collect(keys(matrix_by_prime))
    # Precompute √d mod p for each prime
    sqrtd_by_prime = Dict{Int, Int}()
    for p in primes
        s = compute_sqrt_d_mod_p(d, p)
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
                error("Failed to reconstruct entry ($i, $j) in Z[√$d] with bound=$bound")
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

This function takes the parameter `scale_d`: the d such that the scaling
is √(scale_d). For SU(2)_4 (D² = 12), scale_d = 3.
"""
function reconstruct_S_matrix(group::Dict{Int, MTCCandidate};
                              scale_d::Int = 3,
                              bound::Int = 5)
    primes = collect(keys(group))

    # Build matrix_by_prime: S_scaled[p][i, j] = (2 · √scale_d · S'[i, j]) mod p
    matrix_by_prime = Dict{Int, Matrix{Int}}()
    for p in primes
        s = compute_sqrt_d_mod_p(scale_d, p)
        s === nothing && error("$scale_d is not a QR mod $p")
        two_s = mod(2 * s, p)
        S_p = group[p].S_Fp
        M_p = [mod(two_s * S_p[i, j], p) for i in 1:size(S_p, 1), j in 1:size(S_p, 2)]
        matrix_by_prime[p] = M_p
    end

    return reconstruct_matrix_in_Z_sqrt_d(matrix_by_prime, scale_d; bound = bound)
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
                                scale::Int = 2)
    p = candidate.p
    s = compute_sqrt_d_mod_p(d, p)
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
