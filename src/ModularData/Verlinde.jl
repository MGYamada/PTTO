"""
Verlinde extraction and integer lifting for modular data over F_p.

This file converts validated finite-field modular data into fusion tensors
and ring-agnostic FusionRule values.
"""

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
