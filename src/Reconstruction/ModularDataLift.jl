"""
    ModularDataLift

Bridge from finite-field Phase 3 output to exact cyclotomic modular data.

Phase 3 produces finite-field candidates and reconstructs modular data
directly over `Q(ζ_N)`.

Phase 4 consumes:
- `Nijk :: Array{Int,3}`                 —  used as-is
- `(S, T)` over `Q(ζ_N)`                 —  for exact modular-data checks

This module provides:
- `lift_T_Fp_to_cyclotomic`:  discrete log T_Fp → powers of ζ_N

Depends on: Primes (already in v0.2 via ACMG); no new dependencies.
"""

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
#  T lift: F_p → Q(ζ_N)
# ============================================================

"""
    lift_T_Fp_to_cyclotomic(T_Fp, N, p, zeta_Fp) -> Vector

Lift T-eigenvalues from F_p back to the cyclotomic field `Q(ζ_N)` by
discrete log.

Each entry of T is assumed to be an N-th root of unity (axiomatic for
modular data): T_Fp[i] ∈ ⟨ζ_N⟩ ⊂ F_p*. We recover the exponent k_i
such that T_Fp[i] = ζ_N^{k_i} mod p, then return the exact field element
`ζ_N^k`.

Arguments:
- `T_Fp::Vector{Int}`:   T residues mod p, 1-indexed
- `N::Int`:              conductor (order of ζ_N)
- `p::Int`:              prime (must satisfy N | p-1)
- `zeta_Fp::Int`:        a chosen primitive N-th root of unity in F_p.
                         Must match the one used when reducing T to F_p.

Returns: a vector of elements in `Q(ζ_N)`.

Consistency requirement: The `zeta_Fp` passed here must be THE SAME
primitive root used when producing `T_Fp`. In v0.2's pipeline this is
`find_zeta_in_Fp(N, p)`, which uses `primitive_root(p)` as a fixed base.
"""
function lift_T_Fp_to_cyclotomic(T_Fp::Vector{Int}, N::Int, p::Int, zeta_Fp::Int)
    (p - 1) % N == 0 || error("N=$N does not divide p-1=$(p-1)")
    K, z = cyclotomic_field(N)
    tbl = DiscreteLogTable(N, p, zeta_Fp)
    T = Vector{elem_type(K)}(undef, length(T_Fp))
    for i in eachindex(T_Fp)
        k = discrete_log(tbl, T_Fp[i])
        T[i] = z^k
    end
    return T
end

lift_T_Fp_to_complex(args...; kwargs...) =
    error("lift_T_Fp_to_complex was removed; use lift_T_Fp_to_cyclotomic for exact Q(ζ_N) output")
