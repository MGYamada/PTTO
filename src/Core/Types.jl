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

Base.:(==)(a::FusionRule, b::FusionRule) =
    a.rank == b.rank && a.N == b.N && a.dual == b.dual
Base.hash(fr::FusionRule, h::UInt) = hash((fr.rank, fr.N, fr.dual), h)

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

"""
    fusion_automorphisms(fr_or_Nijk) -> Vector{Vector{Int}}

Return all based-ring automorphisms of a fusion rule, represented as
1-indexed permutations fixing the tensor unit at index `1`.
"""
function fusion_automorphisms(Nijk::Array{Int, 3})
    r = size(Nijk, 1)
    size(Nijk) == (r, r, r) || error("Nijk must be a cube")
    r == 1 && return [Int[1]]
    autos = Vector{Vector{Int}}()
    for rest in _permutations(2:r)
        perm = [1; rest...]
        ok = true
        @inbounds for i in 1:r, j in 1:r, k in 1:r
            if Nijk[perm[i], perm[j], perm[k]] != Nijk[i, j, k]
                ok = false
                break
            end
        end
        ok && push!(autos, perm)
    end
    return autos
end

fusion_automorphisms(fr::FusionRule) = fusion_automorphisms(fr.N)

"""
    is_fusion_automorphism(fr_or_Nijk, perm) -> Bool

Test whether `perm` is a unit-fixing automorphism of a fusion rule.
"""
function is_fusion_automorphism(Nijk::Array{Int, 3}, perm::AbstractVector{<:Integer})
    r = size(Nijk, 1)
    size(Nijk) == (r, r, r) || return false
    length(perm) == r || return false
    collect(sort(perm)) == collect(1:r) || return false
    perm[1] == 1 || return false
    @inbounds for i in 1:r, j in 1:r, k in 1:r
        Nijk[perm[i], perm[j], perm[k]] == Nijk[i, j, k] || return false
    end
    return true
end

is_fusion_automorphism(fr::FusionRule, perm::AbstractVector{<:Integer}) =
    is_fusion_automorphism(fr.N, perm)

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
