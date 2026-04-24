"""
Core data types for modular data and fusion rules.

Design:
- ModularDatumFp: (S, T) over F_p together with metadata (α, D², conductor)
- FusionRule: ring-agnostic, N_ij^k ∈ Z (lifted from F_p, verified non-negative)
- Separation of concerns: ModularDatumFp is a point in M_{r,N} ⊗ F_p;
  FusionRule is a point in F_r (fusion ring scheme)
"""

# ----- Modular datum over F_p -----

"""
    ModularDatumFp

A modular datum realized over F_p.

Fields:
- rank:       number of simple objects r
- p:          prime
- N:          conductor (N | p-1, so ζ_N ∈ F_p)
- S:          r × r symmetric matrix over F_p, with S² = C (charge conj)
- T:          diagonal vector of length r, each entry an N-th root of unity in F_p
- α:          central charge phase, α ∈ μ_{8N}(F_p), satisfying (ST)³ = α · S²
- D_squared:  D² = Σ d_i² ∈ F_p (sum of squared quantum dimensions)
- D:          D ∈ F_p if D² is a square (else 0 as sentinel)
- C:          charge conjugation permutation (as Vector{Int}, 1-indexed)
"""
struct ModularDatumFp
    rank::Int
    p::Int
    N::Int
    S::Matrix{Int}
    T::Vector{Int}
    α::Int
    D_squared::Int
    D::Int
    C::Vector{Int}
end

# ----- Ring-agnostic fusion rule -----

"""
    FusionRule

A fusion rule N_ij^k with non-negative integer coefficients.

Fields:
- rank:         number of simple objects r
- N:           r × r × r array, N[i, j, k] = N_ij^k ∈ Z_{≥0}
- dual:        dual map i → i*, i.e., i + 1 → dual[i+1] (1-indexed)

Invariants (enforced by constructor via `validate`):
- N[1, j, k] = δ_{jk}                  (unit)
- N[i, 1, k] = δ_{ik}                  (unit, other side)
- N[i, j, k] = N[j, i, k]              (commutativity, for MTC fusion)
- N[i, dual[i], 1] = 1, N[i, j, 1] = 0 for j ≠ dual[i]  (duality)
- Associativity: Σ_m N[i,j,m] N[m,k,ℓ] = Σ_m N[j,k,m] N[i,m,ℓ]
"""
struct FusionRule
    rank::Int
    N::Array{Int, 3}
    dual::Vector{Int}
end

"""
    FusionRule(N::Array{Int, 3}) -> FusionRule

Construct a fusion rule from N_ij^k array, deriving dual automatically
and validating axioms. Throws if axioms fail.
"""
function FusionRule(N::Array{Int, 3})
    r = size(N, 1)
    size(N) == (r, r, r) || error("N must be a cube")
    # Derive dual: i* is the unique j with N[i+1, j, 1] = 1 (0-indexed: N_i^0 has support at i*)
    # Here we use 1-indexed: object 1 is unit, so N[i, j, 1] = 1 iff j = dual[i]
    dual = zeros(Int, r)
    for i in 1:r
        found = 0
        for j in 1:r
            if N[i, j, 1] == 1
                found == 0 || error("object $i has ambiguous dual")
                found = j
            elseif N[i, j, 1] != 0
                error("N[$i, $j, 1] = $(N[i,j,1]) not in {0,1}")
            end
        end
        found == 0 && error("object $i has no dual")
        dual[i] = found
    end
    fr = FusionRule(r, N, dual)
    validate(fr)
    return fr
end

"""
    validate(fr::FusionRule)

Verify all axioms. Throws on failure.
"""
function validate(fr::FusionRule)
    r = fr.rank
    N = fr.N
    # Unit axiom
    for j in 1:r, k in 1:r
        expected = (j == k) ? 1 : 0
        N[1, j, k] == expected || error("unit axiom fails at (1, $j, $k): $(N[1,j,k]) ≠ $expected")
        N[j, 1, k] == expected || error("unit axiom fails at ($j, 1, $k): $(N[j,1,k]) ≠ $expected")
    end
    # Commutativity
    for i in 1:r, j in 1:r, k in 1:r
        N[i, j, k] == N[j, i, k] || error("commutativity fails at ($i,$j,$k)")
    end
    # Non-negativity
    all(N .>= 0) || error("fusion coefficients must be non-negative")
    # Associativity
    for i in 1:r, j in 1:r, k in 1:r, ℓ in 1:r
        lhs = sum(N[i, j, m] * N[m, k, ℓ] for m in 1:r)
        rhs = sum(N[j, k, m] * N[i, m, ℓ] for m in 1:r)
        lhs == rhs || error("associativity fails at ($i,$j,$k,$ℓ): $lhs ≠ $rhs")
    end
    return true
end

# ----- Convenience: retrieve fusion matrix N_i -----

"""
    fusion_matrix(fr::FusionRule, i::Int) -> Matrix{Int}

Return the matrix (N_i)_{jk} := N_ij^k, the action of i on the fusion ring.
"""
fusion_matrix(fr::FusionRule, i::Int) = fr.N[i, :, :]

# ----- Equality on fusion rules (up to basis relabeling fixing unit) -----

"""
    fusion_isomorphic(fr1, fr2) -> Bool

Test whether two fusion rules are isomorphic as based rings with unit fixed.
Brute-force over permutations of {2, ..., r} (unit at position 1 fixed).
"""
function fusion_isomorphic(fr1::FusionRule, fr2::FusionRule)
    fr1.rank == fr2.rank || return false
    r = fr1.rank
    r == 1 && return true
    # Try all permutations σ of 2:r (unit 1 fixed)
    for σ_rest in _permutations(2:r)
        σ = [1; σ_rest...]
        match = true
        for i in 1:r, j in 1:r, k in 1:r
            if fr1.N[σ[i], σ[j], σ[k]] != fr2.N[i, j, k]
                match = false
                break
            end
        end
        match && return true
    end
    return false
end

# Simple permutation generator (for small ranks only — up to rank 8 or so)
function _permutations(v)
    v = collect(v)
    length(v) <= 1 && return [v]
    result = Vector{Vector{eltype(v)}}()
    for i in 1:length(v)
        rest = vcat(v[1:i-1], v[i+1:end])
        for p in _permutations(rest)
            push!(result, [v[i]; p])
        end
    end
    return result
end

# ===== Integrated legacy modules =====

"""
Finite field arithmetic over F_p.

Uses Julia's native Int/BigInt arithmetic with `mod(·, p)` wrapping
rather than a wrapper struct, for minimal overhead in tight loops.

All functions assume p is prime and 2 < p < typemax(Int).
"""

# ----- Euler criterion -----
"""
    is_square(a, p) -> Bool

Test if a is a quadratic residue mod p (p odd prime).
Uses Euler's criterion: a^((p-1)/2) ≡ ±1 (mod p).
"""
function is_square(a::Integer, p::Integer)
    a = mod(a, p)
    a == 0 && return true
    return powermod(a, (p - 1) ÷ 2, p) == 1
end

# ----- Tonelli-Shanks: square root mod p -----
"""
    sqrt_mod(a, p) -> Int

Return some r with r^2 ≡ a (mod p), or throw if a is NQR.
Returns the "smaller" of the two roots (r ≤ p/2) for determinism.
"""
function sqrt_mod(a::Integer, p::Integer)
    a = mod(a, p)
    a == 0 && return 0
    is_square(a, p) || error("$a is not a quadratic residue mod $p")

    # p ≡ 3 (mod 4) easy case
    if p % 4 == 3
        r = powermod(a, (p + 1) ÷ 4, p)
        return min(r, p - r)
    end

    # Tonelli-Shanks for p ≡ 1 (mod 4)
    # Decompose p - 1 = q * 2^s with q odd
    q, s = p - 1, 0
    while q % 2 == 0
        q ÷= 2
        s += 1
    end

    # Find a non-residue z
    z = 2
    while is_square(z, p)
        z += 1
    end

    m = s
    c = powermod(z, q, p)
    t = powermod(a, q, p)
    r = powermod(a, (q + 1) ÷ 2, p)

    while true
        t == 1 && return min(r, p - r)
        # find least i, 0 < i < m, such that t^(2^i) = 1
        i, temp = 0, t
        while temp != 1
            temp = (temp * temp) % p
            i += 1
            i >= m && error("Tonelli-Shanks failed (shouldn't happen)")
        end
        b = powermod(c, 1 << (m - i - 1), p)
        m = i
        c = (b * b) % p
        t = (t * c) % p
        r = (r * b) % p
    end
end

# ----- Primitive root of F_p^* -----
"""
    primitive_root(p) -> Int

Return the smallest positive primitive root modulo p.
"""
function primitive_root(p::Integer)
    p == 2 && return 1
    order_target = p - 1
    # Factor p - 1 to check orders (qualify as Primes.factor to avoid clash with Oscar)
    factors = [f for (f, _) in Primes.factor(order_target)]
    for g in 2:(p-1)
        is_primitive = true
        for q in factors
            if powermod(g, order_target ÷ q, p) == 1
                is_primitive = false
                break
            end
        end
        is_primitive && return Int(g)
    end
    error("no primitive root found mod $p (shouldn't happen for prime p)")
end

# ----- n-th root of unity -----
"""
    root_of_unity(n, p) -> Int

Return a primitive n-th root of unity in F_p.
Requires n | (p - 1).
"""
function root_of_unity(n::Integer, p::Integer)
    n == 1 && return 1
    (p - 1) % n == 0 || error("ζ_$n does not exist in F_$p (need $n | $(p-1))")
    g = primitive_root(p)
    return powermod(g, (p - 1) ÷ n, p)
end

"""
    roots_of_unity(n, p) -> Vector{Int}

Return all n-th roots of unity in F_p, i.e. [ζ_n^0, ζ_n^1, ..., ζ_n^{n-1}].
"""
function roots_of_unity(n::Integer, p::Integer)
    ζ = root_of_unity(n, p)
    result = Vector{Int}(undef, n)
    result[1] = 1
    for k in 2:n
        result[k] = (result[k-1] * ζ) % p
    end
    return result
end

"""
    primitive_roots_of_unity(n, p) -> Vector{Int}

Return all *primitive* n-th roots of unity in F_p.
These are ζ_n^k for gcd(k, n) = 1.
"""
function primitive_roots_of_unity(n::Integer, p::Integer)
    all_roots = roots_of_unity(n, p)
    return [all_roots[k+1] for k in 1:(n-1) if gcd(k, n) == 1]
end

# ----- Matrix helpers over F_p -----
"""
    matmul_mod(A, B, p) -> Matrix{Int}

Matrix multiplication mod p.
"""
function matmul_mod(A::AbstractMatrix{<:Integer}, B::AbstractMatrix{<:Integer}, p::Integer)
    m, n = size(A)
    n2, k = size(B)
    n == n2 || error("matrix dimensions mismatch: $(size(A)) * $(size(B))")
    C = zeros(Int, m, k)
    for i in 1:m, j in 1:k
        s = 0
        for ℓ in 1:n
            s = (s + Int(A[i, ℓ]) * Int(B[ℓ, j])) % p
        end
        C[i, j] = s
    end
    return C
end

"""
    matpow_mod(A, n, p) -> Matrix{Int}

Compute A^n mod p for square matrix A, using fast exponentiation.
"""
function matpow_mod(A::AbstractMatrix{<:Integer}, n::Integer, p::Integer)
    n == 0 && return Matrix{Int}(I, size(A, 1), size(A, 1))
    n == 1 && return mod.(Int.(A), p)
    result = Matrix{Int}(I, size(A, 1), size(A, 1))
    base = mod.(Int.(A), p)
    while n > 0
        if n & 1 == 1
            result = matmul_mod(result, base, p)
        end
        n >>= 1
        if n > 0
            base = matmul_mod(base, base, p)
        end
    end
    return result
end

"""
    diagmul_right(A, d, p) -> Matrix{Int}

Compute A * diag(d) mod p. Scales column j of A by d[j].
"""
function diagmul_right(A::AbstractMatrix{<:Integer}, d::AbstractVector{<:Integer}, p::Integer)
    m, n = size(A)
    length(d) == n || error("diagonal length mismatch")
    C = zeros(Int, m, n)
    for i in 1:m, j in 1:n
        C[i, j] = (Int(A[i, j]) * Int(d[j])) % p
    end
    return C
end

"""
    diagmul_left(d, A, p) -> Matrix{Int}

Compute diag(d) * A mod p. Scales row i of A by d[i].
"""
function diagmul_left(d::AbstractVector{<:Integer}, A::AbstractMatrix{<:Integer}, p::Integer)
    m, n = size(A)
    length(d) == m || error("diagonal length mismatch")
    C = zeros(Int, m, n)
    for i in 1:m, j in 1:n
        C[i, j] = (Int(d[i]) * Int(A[i, j])) % p
    end
    return C
end

"""
    lift_symmetric(x, p) -> Int

Lift x ∈ F_p = {0,...,p-1} to the symmetric range [-p/2, p/2].
Useful for recognizing small integers (e.g., fusion coefficients).
"""
function lift_symmetric(x::Integer, p::Integer)
    x = mod(x, p)
    return x <= p ÷ 2 ? Int(x) : Int(x) - Int(p)
end

"""
Modular data validation and construction over F_p.

A modular datum (S, T) must satisfy:
- S is symmetric: S[i,j] = S[j,i]
- S² = C (charge conjugation permutation matrix), C² = I
- T is diagonal with T[0,0] = 1, entries are N-th roots of unity
- (ST)³ = α · S² for some α ∈ μ_{8N}(F_p) (central charge phase)
- S[i, 0] ≠ 0 for all i (so quantum dimensions are well-defined)

This module provides:
- `validate_modular_data`: check all axioms given (S, T, p, N)
- `build_modular_datum`: construct ModularDatumFp from (S, T, p, N) with validation
- `compute_alpha`: compute α from (ST)³ and S²
- `compute_charge_conjugation`: extract C from S²
"""

"""
    compute_charge_conjugation(S::Matrix{Int}, p::Int) -> Union{Vector{Int}, Nothing}

Given S mod p, compute S² and extract charge conjugation permutation C,
where S² = C (as permutation matrix). Returns C as a vector where C[i] is
the image of i under the permutation (1-indexed).

Returns nothing if S² is not a permutation matrix.
"""
function compute_charge_conjugation(S::AbstractMatrix{Int}, p::Int)
    r = size(S, 1)
    S2 = matmul_mod(S, S, p)
    C = zeros(Int, r)
    for i in 1:r
        found = 0
        for j in 1:r
            if S2[i, j] == 1
                found == 0 || return nothing  # multiple 1s in row
                found = j
            elseif S2[i, j] != 0
                return nothing  # non-0, non-1 entry
            end
        end
        found == 0 && return nothing  # no 1 in row
        C[i] = found
    end
    # Check it's an involution
    for i in 1:r
        C[C[i]] == i || return nothing
    end
    return C
end

"""
    compute_alpha(S::Matrix{Int}, T::Vector{Int}, p::Int) -> Union{Int, Nothing}

Compute α such that (ST)³ = α · S². Returns α ∈ F_p if (ST)³ is indeed a
scalar multiple of S² (entry-wise), otherwise nothing.
"""
function compute_alpha(S::AbstractMatrix{Int}, T::AbstractVector{Int}, p::Int)
    r = size(S, 1)
    ST = diagmul_right(S, T, p)
    ST3 = matpow_mod(ST, 3, p)
    S2 = matmul_mod(S, S, p)

    # Find α: (ST)³[i,j] = α * S²[i,j] for all i,j
    α = -1  # sentinel
    for i in 1:r, j in 1:r
        lhs = ST3[i, j]
        rhs = S2[i, j]
        if rhs == 0
            lhs == 0 || return nothing
        else
            α_candidate = (lhs * invmod(rhs, p)) % p
            if α == -1
                α = α_candidate
            elseif α != α_candidate
                return nothing  # inconsistent
            end
        end
    end
    return α == -1 ? 1 : α  # if S² ≡ 0 (shouldn't happen for modular data), default 1
end

"""
    validate_modular_data(S, T, p, N) -> NamedTuple

Perform all modular datum axiom checks. Returns a NamedTuple:
- valid::Bool
- reason::String (empty if valid, else describes first failure)
- α::Int (central charge phase, if computable)
- C::Vector{Int} (charge conjugation, if computable)
- D_squared::Int
- D::Int (0 if D² is not a square)
"""
function validate_modular_data(S::AbstractMatrix{Int}, T::AbstractVector{Int}, p::Int, N::Int)
    r = size(S, 1)
    # Shape checks
    size(S) == (r, r) || return (valid=false, reason="S not square", α=0, C=Int[], D_squared=0, D=0)
    length(T) == r || return (valid=false, reason="T wrong length", α=0, C=Int[], D_squared=0, D=0)
    # Symmetry
    for i in 1:r, j in (i+1):r
        S[i, j] == S[j, i] || return (valid=false, reason="S not symmetric at ($i,$j)", α=0, C=Int[], D_squared=0, D=0)
    end
    # T[1] = 1 (unit)
    T[1] == 1 || return (valid=false, reason="T[1] ≠ 1", α=0, C=Int[], D_squared=0, D=0)
    # T entries are N-th roots of unity
    for i in 1:r
        powermod(T[i], N, p) == 1 || return (valid=false, reason="T[$i] is not an N-th root of unity", α=0, C=Int[], D_squared=0, D=0)
    end
    # S[1, i] ≠ 0 (quantum dimensions well-defined)
    for i in 1:r
        S[1, i] != 0 || return (valid=false, reason="S[1,$i] = 0", α=0, C=Int[], D_squared=0, D=0)
    end
    # S² = C (permutation matrix, involution)
    C = compute_charge_conjugation(S, p)
    C === nothing && return (valid=false, reason="S² is not a permutation matrix", α=0, C=Int[], D_squared=0, D=0)
    # (ST)³ = α · S²
    α = compute_alpha(S, T, p)
    α === nothing && return (valid=false, reason="(ST)³ not proportional to S²", α=0, C=C, D_squared=0, D=0)
    # Compute D² from quantum dimensions
    # d_i = S[i, 1] / S[1, 1], D² = Σ d_i²
    S00_inv = invmod(S[1, 1], p)
    D_squared = 0
    for i in 1:r
        d_i = (Int(S[i, 1]) * S00_inv) % p
        D_squared = (D_squared + d_i * d_i) % p
    end
    # D is sqrt of D² if exists
    D = is_square(D_squared, p) ? sqrt_mod(D_squared, p) : 0
    return (valid=true, reason="", α=α, C=C, D_squared=D_squared, D=D)
end

"""
    build_modular_datum(S, T, p, N) -> Union{ModularDatumFp, Nothing}

Construct a ModularDatumFp from raw (S, T, p, N) with full validation.
Returns nothing if any axiom fails.
"""
function build_modular_datum(S::AbstractMatrix{Int}, T::AbstractVector{Int}, p::Int, N::Int)
    r = size(S, 1)
    S_mod = mod.(Int.(S), p)
    T_mod = mod.(Int.(T), p)
    result = validate_modular_data(S_mod, T_mod, p, N)
    result.valid || return nothing
    return ModularDatumFp(r, p, N, S_mod, T_mod, result.α, result.D_squared, result.D, result.C)
end

"""
Extract fusion rule from modular datum over F_p via Verlinde formula.

Key formula (BNRW convention):
    N_ij^k = Σ_m S_im * S_jm * S*_km / S_0m

where S*_km is the entry of S† (Hermitian conjugate). For a real S-matrix
in the standard MTC sense, S* = S̄, and over F_p this becomes the Galois
conjugate under Frob_p. For our p-fix enumeration with all objects self-dual
or when we work in the representation where S is symmetric and entries are
in μ_∞, S̄ is related to S via charge conjugation: S̄ = CS.

For simplicity in the prototype, we assume the conjugation is given by
the precomputed C (charge conjugation permutation): S†_km = S_{C[k], m}.
This is correct when all entries of S lie in a field where conjugation acts
as the charge conjugation permutation.
"""

"""
    verlinde_coefficient(md::ModularDatumFp, i, j, k) -> Int

Compute N_ij^k ∈ F_p via Verlinde formula.
Assumes 1-indexed i, j, k ∈ 1:r.
"""
function verlinde_coefficient(md::ModularDatumFp, i::Int, j::Int, k::Int)
    r = md.rank
    p = md.p
    S = md.S
    C = md.C  # charge conjugation: C[k] is the index of k-bar
    kbar = C[k]  # k* in 1-indexed
    # S^†_km in our convention = S_{k*, m}
    total = 0
    for m in 1:r
        S0m_inv = invmod(S[1, m], p)  # S_0m ≠ 0 for modular data
        term = (Int(S[i, m]) * Int(S[j, m])) % p
        term = (term * Int(S[kbar, m])) % p
        term = (term * S0m_inv) % p
        total = (total + term) % p
    end
    return total
end

"""
    extract_fusion_rule_Fp(md::ModularDatumFp) -> Array{Int, 3}

Compute all N_ij^k ∈ F_p. Returns the full rank × rank × rank array.
"""
function extract_fusion_rule_Fp(md::ModularDatumFp)
    r = md.rank
    N_Fp = zeros(Int, r, r, r)
    for i in 1:r, j in 1:r, k in 1:r
        N_Fp[i, j, k] = verlinde_coefficient(md, i, j, k)
    end
    return N_Fp
end

"""
    lift_fusion_to_Z(N_Fp::Array{Int, 3}, p::Int; bound::Int = 0) -> Union{Array{Int, 3}, Nothing}

Attempt to lift F_p fusion coefficients to non-negative integers.
Each coefficient x ∈ {0, ..., p-1} is accepted as itself if x ≤ bound,
otherwise as x - p if p - x ≤ bound (but then it's negative — rejected).

If `bound` is 0 (default), use p ÷ 4 as a conservative heuristic.

Returns the integer array if all coefficients lift to non-negative values
within the bound, otherwise nothing.
"""
function lift_fusion_to_Z(N_Fp::Array{Int, 3}, p::Int; bound::Int = 0)
    bound == 0 && (bound = p ÷ 4)
    r = size(N_Fp, 1)
    N_Z = zeros(Int, r, r, r)
    for i in 1:r, j in 1:r, k in 1:r
        x = N_Fp[i, j, k]
        if 0 <= x <= bound
            N_Z[i, j, k] = x
        else
            return nothing  # cannot lift to small non-negative integer
        end
    end
    return N_Z
end

"""
    extract_and_lift(md::ModularDatumFp; bound::Int = 0) -> Union{FusionRule, Nothing}

Full pipeline: compute Verlinde coefficients, lift to Z, wrap as FusionRule
with axiom validation. Returns nothing if any step fails.
"""
function extract_and_lift(md::ModularDatumFp; bound::Int = 0)
    N_Fp = extract_fusion_rule_Fp(md)
    N_Z = lift_fusion_to_Z(N_Fp, md.p; bound = bound)
    N_Z === nothing && return nothing
    try
        return FusionRule(N_Z)
    catch
        return nothing  # axiom violation
    end
end
