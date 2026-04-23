"""
    ModularDataLift

Bridge from ACMG v0.2 Phase 3 output (F_p + ℤ[√d]) to Phase 4 input (ℂ).

v0.2 Phase 3 produces two kinds of data per MTC candidate:
- `MTCCandidate.T_Fp :: Vector{Int}`  —  T-eigenvalues mod p
- `reconstruct_S_matrix(group) :: Matrix{Tuple{Int,Int}}`  —  entries
  `(a, b)` meaning `a + b·√d` in ℤ[√d]
- `MTCCandidate.N :: Array{Int,3}`  —  fusion coefficients (Galois-invariant)

Phase 4 (pentagon/hexagon solver) consumes:
- `Nijk :: Array{Int,3}`                 —  used as-is
- `(S, T) :: (Matrix{ComplexF64}, Vector{ComplexF64})` —  for Verify layer
  (not strictly needed for pentagon itself, but needed for round-trip check)

This module provides:
- `lift_T_Fp_to_complex`:  discrete log T_Fp → ζ_N powers → ComplexF64
- `lift_S_sqrtd_to_complex`: (a, b) in ℤ[√d] → a + b·√d as ComplexF64
- `lift_mtc_candidate`:    convenience wrapper bundling both

Depends on: Primes (already in v0.2 via ACMG); no new dependencies.
"""

using LinearAlgebra
using Primes

# Access v0.2 utilities. We must reach into the outer ACMG module.
# These functions are expected to be defined at ACMG top level:
#   find_zeta_in_Fp(N, p) :: Int
#   primitive_root(p) :: Int
# both live in ACMG.BlockU / ACMG.FpArith.


# ============================================================
#  Discrete log infrastructure
# ============================================================

"""
    DiscreteLogTable(N, p, zeta_Fp)

Precomputed table mapping ζ_N^k (mod p) → k for k = 0, …, N-1.

Building this table is O(N) work; looking up is O(1). Used to invert
`cyclotomic_to_Fp` when we only know the F_p image.

Fields:
- N, p, zeta_Fp: parameters used to build the table
- lookup: Dict{Int, Int} mapping F_p residue → exponent in [0, N)
"""
struct DiscreteLogTable
    N::Int
    p::Int
    zeta_Fp::Int
    lookup::Dict{Int, Int}
end

function DiscreteLogTable(N::Int, p::Int, zeta_Fp::Int)
    lookup = Dict{Int, Int}()
    acc = 1
    for k in 0:(N - 1)
        if haskey(lookup, acc)
            error("ζ_$N has order < $N in F_$p: ζ^$k = ζ^$(lookup[acc]) = $acc. " *
                  "zeta_Fp may not be a primitive N-th root of unity.")
        end
        lookup[acc] = k
        acc = mod(acc * zeta_Fp, p)
    end
    acc == 1 || error("ζ_$N does not have order $N in F_$p (loop closed at $acc ≠ 1)")
    return DiscreteLogTable(N, p, zeta_Fp, lookup)
end

"""
    discrete_log(tbl::DiscreteLogTable, x::Int) -> Int

Return k ∈ [0, N) such that ζ_N^k ≡ x (mod p), using the precomputed
lookup. Throws if `x` is not in the ⟨ζ_N⟩ cyclic subgroup.

Note: only monomials ζ_N^k are invertible this way. General elements
`a + b·ζ + c·ζ² + …` cannot be recovered without extra structure
(e.g., knowing we're in a cyclic subgroup, as is the case for T).
"""
function discrete_log(tbl::DiscreteLogTable, x::Int)
    haskey(tbl.lookup, mod(x, tbl.p)) ||
        error("F_$(tbl.p) value $x is not a power of ζ_$(tbl.N); " *
              "cannot discrete-log. " *
              "(This value has components outside the ⟨ζ_N⟩ subgroup.)")
    return tbl.lookup[mod(x, tbl.p)]
end

# ============================================================
#  T lift: F_p → ℂ
# ============================================================

"""
    lift_T_Fp_to_complex(T_Fp, N, p, zeta_Fp) -> Vector{ComplexF64}

Lift T-eigenvalues from F_p back to ℂ by discrete log.

Each entry of T is assumed to be an N-th root of unity (axiomatic for
modular data): T_Fp[i] ∈ ⟨ζ_N⟩ ⊂ F_p*. We recover the exponent k_i
such that T_Fp[i] = ζ_N^{k_i} mod p, then embed as
`exp(2π·i·k_i / N) ∈ ℂ`.

Arguments:
- `T_Fp::Vector{Int}`:   T residues mod p, 1-indexed
- `N::Int`:              conductor (order of ζ_N)
- `p::Int`:              prime (must satisfy N | p-1)
- `zeta_Fp::Int`:        a chosen primitive N-th root of unity in F_p.
                         Must match the one used when reducing T to F_p.

Returns: Vector{ComplexF64} of the same length, each entry on the unit
circle with argument 2π·k/N for some integer k ∈ [0, N).

Consistency requirement: The `zeta_Fp` passed here must be THE SAME
primitive root used when producing `T_Fp`. In v0.2's pipeline this is
`find_zeta_in_Fp(N, p)`, which uses `primitive_root(p)` as a fixed base.
"""
function lift_T_Fp_to_complex(T_Fp::Vector{Int}, N::Int, p::Int, zeta_Fp::Int)
    (p - 1) % N == 0 || error("N=$N does not divide p-1=$(p-1)")
    tbl = DiscreteLogTable(N, p, zeta_Fp)
    T_complex = Vector{ComplexF64}(undef, length(T_Fp))
    two_pi_i_over_N = (2π * im) / N
    for i in eachindex(T_Fp)
        k = discrete_log(tbl, T_Fp[i])
        T_complex[i] = exp(two_pi_i_over_N * k)
    end
    return T_complex
end

# ============================================================
#  S lift: ℤ[√d] → ℂ
# ============================================================

"""
    lift_S_sqrtd_to_complex(recon_S, d; scale=2) -> Matrix{ComplexF64}

Convert a matrix of ℤ[√d] entries (as `(a, b)` tuples meaning
`a + b·√d`) into a ℂ matrix, undoing any normalization scaling.

`reconstruct_S_matrix` in v0.2 returns `scale · √d · S'` for some
`scale` (default 2 in the SU(2)_4 / d=3 setting). To recover `S'` we
divide by `scale · √d`.

Arguments:
- `recon_S::Matrix{Tuple{Int,Int}}`:   each entry `(a, b)` = `a + b√d`
- `d::Int`:                             d under the square root
- `scale::Int=2`:                       normalization factor used at
                                        reconstruction time
                                        (v0.2 default, see
                                         `reconstruct_S_matrix` docstring)

Returns: `Matrix{ComplexF64}` representing `S'` as originally a
cyclotomic-valued matrix.

Note: The choice of branch for √d matches the `sqrtd_fn` used at
reconstruction time. Here we simply take the positive real branch
`sqrt(d)` for the ℂ embedding; Galois conjugates should already be
encoded in the sign of `b`.
"""
function lift_S_sqrtd_to_complex(recon_S::Matrix{Tuple{Int,Int}}, d::Int;
                                 scale::Int = 2)
    nrow, ncol = size(recon_S)
    sd = sqrt(Float64(d))
    denom = scale * sd
    result = Matrix{ComplexF64}(undef, nrow, ncol)
    for i in 1:nrow, j in 1:ncol
        (a, b) = recon_S[i, j]
        result[i, j] = ComplexF64((a + b * sd) / denom)
    end
    return result
end

# ============================================================
#  Convenience wrapper
# ============================================================

"""
    lift_mtc_candidate(candidate, recon_S; d, scale=2) ->
        (S::Matrix{ComplexF64}, T::Vector{ComplexF64}, Nijk::Array{Int,3})

Combine the two lifts into a single call that produces Phase 4-ready
modular data from one v0.2 MTCCandidate plus a reconstructed S.

Arguments:
- `candidate`:    an `ACMG.MTCCandidate` from Phase 2 at some prime `p`.
                  Must carry `.p`, `.T_Fp`, `.N` (fusion tensor).
- `recon_S`:      the matrix from `reconstruct_S_matrix(group; …)`.
- `d::Int`:       the d such that entries live in ℤ[√d] (keyword).
- `scale::Int=2`: scaling factor used at reconstruction (keyword).

Returns `(S, T, Nijk)` suitable for calling Phase 4 routines:
- `Nijk = candidate.N` (pass directly to `get_pentagon_system`)
- `(S, T)` for downstream verification

Requirements at the call site:
- `find_zeta_in_Fp` and `primitive_root` must be accessible.
  If called from within the ACMG package, use
  `ACMG.find_zeta_in_Fp(N, candidate.p)` and pass as `zeta_Fp`.
  The N (conductor) must be available — it's not stored on
  MTCCandidate directly, so it must be passed explicitly.
"""
function lift_mtc_candidate(candidate, recon_S::Matrix{Tuple{Int,Int}};
                            d::Int, N::Int, zeta_Fp::Int, scale::Int = 2)
    S_ℂ = lift_S_sqrtd_to_complex(recon_S, d; scale = scale)
    T_ℂ = lift_T_Fp_to_complex(candidate.T_Fp, N, candidate.p, zeta_Fp)
    Nijk = candidate.N
    return (S_ℂ, T_ℂ, Nijk)
end

