# mtc_types.jl
# Type system aligned with Ng-Rowell-Wang-Wen (arXiv:2203.14829).
# April 16, 2026 — Refactoring stage 1
#
# This file introduces the paper's vocabulary as Julia types. It deliberately
# duplicates no logic from mtc_classifier.jl; instead it offers clean
# constructors that wrap existing matrices, plus predicates that the old
# enumerate_modular_data path doesn't yet use. Once these are stable, the
# classifier can be rewritten to speak in these types.
#
# Paper reference: every docstring cites the specific result (Theorem 2.1,
# 3.4, 3.7, Remark 3.5/3.8, Proposition 3.9, Definition 3.1/3.2).

using LinearAlgebra

# ============================================================
#  Stage 1 types: objects produced by each step of the paper
# ============================================================

"""
    AtomicIrrep

One irreducible `SL₂(ℤ/p^a ℤ)`-representation from Appendix A of the paper,
or a tensor product of such over distinct prime-power levels (Chinese
Remainder Theorem, §3.5 of the paper).

Fields follow the paper's notation `d^{a,m}_{l,k}`:
  `dim`    – dimension d of the representation
  `level`  – l = ord(ρ(t)); always a positive integer
  `label`  – human-readable identifier such as "2_5^1" or "2_5 ⊗ 3_7"
  `S, T`   – matrices over the cyclotomic field `K = Q(ζ_N)`, where N is
             the conductor in which the classifier is working
  `parity` – +1 if ρ(s²) = +I (even, `C` is self-dual-compatible),
             -1 if ρ(s²) = -I (odd), 0 if neither (mixed — only possible
             for tensor products of mixed-parity factors; treat with care)
  `K`      – the cyclotomic field the entries live in
  `N`      – conductor (so that K = Q(ζ_N) and N is a multiple of level)

Remark 3.5 of the paper: if ρ is symmetric and irreducible, then ρ(s) is
either real-symmetric (even) or purely imaginary (odd). Our GAP-sourced
atomic irreps are always pure even or pure odd. Direct sums of atomics
with opposite parities are neither — those are `IrrepSumRep`s below.
"""
struct AtomicIrrep
    dim::Int
    level::Int
    label::String
    S                    # MatElem over K
    T                    # MatElem over K, diagonal
    parity::Int
    K                    # AbsSimpleNumField = Q(ζ_N)
    N::Int
end

"""
    IrrepSumRep

A direct sum of atomic irreps, as produced in Step 1 of the paper's
classification strategy (§3.4). The "type" (Definition 3.2) is the
sorted list of atomic dimensions — a Young diagram with `rank` boxes.

An `IrrepSumRep` is **symmetric** (Theorem 3.3) and **congruence** but
is generally not yet in the right basis to be modular data: it is
block-diagonal, whereas a genuine MD matrix has nonzero entries in the
unit row because every quantum dimension is nonzero.

Fields:
  `components` – vector of (atomic_index, multiplicity) pairs
                 referring to a catalog of AtomicIrrep
  `type`       – Young diagram: sorted list of atomic dims, non-increasing
  `S, T`       – block-diagonal matrices over K, size = rank × rank
  `rank`       – total rank = sum of component dims × multiplicities
  `K, N`       – as in AtomicIrrep
"""
struct IrrepSumRep
    components::Vector{Tuple{Int,Int}}   # (catalog index, multiplicity)
    type::Vector{Int}                    # Young diagram, sorted ≥
    S                                    # block-diagonal MatElem over K
    T                                    # block-diagonal MatElem over K
    rank::Int
    K
    N::Int
end

"""
    OrderedIrrepSum

An `IrrepSumRep` conjugated by a permutation matrix P so that the
diagonal entries of T, interpreted as topological spins via
`arg(T[i,i]) = 2π s_i`, are sorted in a canonical order.

Corresponds to equation (3.13) of the paper: ρ̃ = P ρ_isum P^T.

The permutation `perm` maps new indices to old: so T_new[i,i] =
T_old[perm[i], perm[i]].
"""
struct OrderedIrrepSum
    base::IrrepSumRep
    S                 # = P · S · P^T
    T                 # = P · T · P^T (diagonal, sorted)
    perm::Vector{Int}
end

"""
    pMDRep (pseudo-modular-data representation)

Definition 3.1 of the paper: ρ is a pMD representation if V·ρ·V is an MD
representation for some signed diagonal V. Equivalently, it is the image
of an OrderedIrrepSum under an orthogonal U that commutes with T
(Theorem 3.4). Structurally:

    ρ_pMD = U · ρ̃ · U^T

where U is *block-diagonal orthogonal* with block sizes equal to the
multiplicities of the eigenvalues of ρ̃(t). The blocks corresponding to
non-degenerate eigenvalues are just ±1 (absorbed into V later); the
interesting blocks are for degenerate eigenvalues, where U is a genuine
SO(k) element.
"""
struct pMDRep
    base::OrderedIrrepSum
    U                 # real orthogonal matrix (in K if rational-valued)
    S                 # = U · base.S · U^T
    T                 # = base.T (unchanged because U commutes with T)
    K
    N::Int
end

"""
    MDRep (modular-data representation)

The final object: a pair (S, T) that satisfies all conditions of
Theorem 2.1 of the paper. Constructed from a pMDRep by conjugating with
a signed diagonal matrix V (Remark 3.8), chosen so that the quantum
dimensions `d_i = S[1,i] / S[1,1]` are all positive real.

Fields include the derived invariants that downstream code (Pentagon /
Hexagon solvers) needs: fusion coefficients, central charge mod 8,
global quantum dimension squared.
"""
struct MDRep
    base::pMDRep
    V::Vector{Int}                  # ±1 entries, the signed diagonal
    S                               # = V · base.S · V (entry-wise)
    T                               # = base.T
    Nijk::Array{Int,3}              # fusion coefficients from Verlinde
    central_charge::Rational{Int}   # c mod 8 (numerator over 8)
    K
    N::Int
end

# ============================================================
#  Predicates — one per Theorem 2.1 condition, plus paper's §3
# ============================================================

"""
    compute_parity(ρ_S; tol=1e-10)

Evaluate ρ(s²) and classify:
   +1 if ρ(s²) = +I (representation is **even**, per Remark 3.5);
   -1 if ρ(s²) = -I (representation is **odd**);
    0 otherwise (mixed — expected only for direct sums with opposite-parity atomics).

The caller is expected to pass the numeric form of ρ(s) as a complex
matrix. This function does not attempt exact arithmetic because the
purpose of parity is classification, not certification — the exact
condition comes later via Theorem 3.4-based checks.
"""
function compute_parity(ρ_S::AbstractMatrix{<:Number}; tol::Real=1e-10)
    s² = ρ_S * ρ_S
    r = size(s², 1)
    diff_plus  = s² - Matrix(I, r, r)
    diff_minus = s² + Matrix(I, r, r)
    maximum(abs, diff_plus)  < tol && return +1
    maximum(abs, diff_minus) < tol && return -1
    return 0
end

"""
    has_dim_real(S, K, N; tol=1e-10) -> Bool

Theorem 2.1(3): the quantum dimensions `d_i = S[1,i]/S[1,1]` must be
**real** numbers. Equivalently, S[1,i]/S[1,1] = conj(S[1,i]/S[1,1]) in
Q(ζ_N).

This fails when e.g. the atomic 2_5^1 is used in its raw GAP form,
because the factor (s_5^1)^{-1} = 1/(ζ_5 - ζ_5^{-1}) has not been pulled
out. In paper's language, that atomic is *odd* and its S is purely
imaginary; normalization brings it to real form.

Returns true if all ratios are real, false otherwise.
"""
function has_dim_real_numeric(S_num::AbstractMatrix{<:Number}; tol::Real=1e-10)
    r = size(S_num, 1)
    S11 = S_num[1, 1]
    abs(S11) < tol && return false
    for i in 1:r
        d_i = S_num[1, i] / S11
        abs(imag(d_i)) > tol && return false
    end
    return true
end

"""
    fp_ratio_nonnegative(S, tol=1e-10) -> Bool

Theorem 2.1(4), last clause: there exists an index j (the
Frobenius-Perron row) such that S[i,j]/S[0,j] ∈ [1, +∞) for every i.

In a valid MD, this is the Perron-Frobenius eigenvector row of the
fusion matrices; picking S[:, 0] (the unit column) normally works since
S[i, 0] = d_i ≥ 1 for every object in a pseudo-unitary MTC.

Returns true if *some* column j gives all ratios ≥ 1. Uses float
arithmetic (this is a consistency check, not an algebraic certificate).
"""
function fp_ratio_nonnegative(S_num::AbstractMatrix{<:Number}; tol::Real=1e-10)
    r = size(S_num, 1)
    for j in 1:r
        s0j = S_num[1, j]
        abs(s0j) < tol && continue
        ok = true
        for i in 1:r
            ratio = S_num[i, j] / s0j
            if abs(imag(ratio)) > tol || real(ratio) < 1 - tol
                ok = false
                break
            end
        end
        ok && return true
    end
    return false
end

"""
    non_degenerate_block(S_num, T_num; tol=1e-10)

Return (indices, block) where `indices` is the list of positions i such
that T[i,i] is a non-degenerate (multiplicity-1) eigenvalue of T, and
`block = S[indices, indices]`.

Proposition 3.9 of the paper: for any (symmetric) SL₂(ℤ) representation
ρ̃ that is equivalent to an MD representation, the entries of ρ̃(s)_{ndeg}
are cyclotomic numbers in Q_{ord(ρ̃)}. This gives a computationally cheap
test: after extracting the non-degenerate block, check that every entry
lies in Q_N (not just in some larger field).

The caller is responsible for the Q_N membership test, which in Oscar
amounts to confirming that the minimal polynomial divides the N-th
cyclotomic polynomial.
"""
function non_degenerate_block(S_num::AbstractMatrix{<:Number},
                              T_num::AbstractMatrix{<:Number};
                              tol::Real=1e-10)
    r = size(T_num, 1)
    t_diag = [T_num[i, i] for i in 1:r]
    # Group equal eigenvalues
    multiplicities = zeros(Int, r)
    for i in 1:r, j in 1:r
        abs(t_diag[i] - t_diag[j]) < tol && (multiplicities[i] += 1)
    end
    indices = findall(m -> m == 1, multiplicities)
    block = S_num[indices, indices]
    return indices, block
end

# ============================================================
#  Stage 0: wrap existing atomic catalog entries
# ============================================================

"""
    atomic_from_catalog_entry(entry, K, N)

Given a NamedTuple `(S, T, rank, name)` from the existing
`inspect_atomic_irreps` / `enumerate_modular_data` pipeline, build an
`AtomicIrrep`. Computes parity numerically.

The `name` field from the existing pipeline is just "2d" or "2d⊗3d"
and doesn't carry the full `d^{a,m}_{l,k}` information; we keep it as-is
for now. A full label-aware catalog is a later refactoring step (see
implementation_plan.md step [B]).
"""
function atomic_from_catalog_entry(entry, K, N::Int;
                                   matrix_to_complex_fn, degree_fn)
    r = entry.rank
    zeta = exp(2π * im / N)
    deg = degree_fn(K)
    S_num = matrix_to_complex_fn(entry.S, zeta, deg)
    T_num = matrix_to_complex_fn(entry.T, zeta, deg)
    par = compute_parity(S_num)

    # Infer level from T (ord of the diagonal)
    level = _numeric_order(T_num, r, N)

    return AtomicIrrep(
        r,            # dim
        level,        # level
        entry.name,   # label (coarse for now)
        entry.S,      # S over K
        entry.T,      # T over K
        par,          # parity
        K,
        N,
    )
end

function _numeric_order(T_num::AbstractMatrix{<:Number}, r::Int, N::Int;
                       tol::Real=1e-8)
    # The smallest m ≥ 1 with T^m = I.  Since we're over Q(ζ_N), m | N.
    for m in 1:N
        ok = true
        for i in 1:r
            if abs(T_num[i, i]^m - 1) > tol
                ok = false
                break
            end
        end
        ok && return m
    end
    return N  # fallback
end

# ============================================================
#  Pretty printing
# ============================================================

function Base.show(io::IO, a::AtomicIrrep)
    par_str = a.parity == +1 ? "even" :
              a.parity == -1 ? "odd"  : "mixed"
    print(io, "AtomicIrrep(dim=$(a.dim), level=$(a.level), $par_str, label=$(a.label))")
end

function Base.show(io::IO, s::IrrepSumRep)
    type_str = "(" * join(s.type, ",") * ")"
    print(io, "IrrepSumRep(rank=$(s.rank), type=$type_str, N=$(s.N))")
end

function Base.show(io::IO, p::pMDRep)
    print(io, "pMDRep(rank=$(p.base.base.rank), N=$(p.N))")
end

function Base.show(io::IO, m::MDRep)
    print(io, "MDRep(rank=$(m.base.base.base.rank), c=$(m.central_charge)/8, N=$(m.N))")
end
