"""
SL(2, ℤ/N) irreducible representation catalog (atomic irreps).

This module provides:
- `AtomicIrrep`: data structure for a single `SL(2, ℤ/N)` irrep realized over
  `Q(ζ_N)`
- `build_atomic_catalog`: enumerate irreps needed by the conductor-first
  classification pipeline using GAP's SL2Reps package through Oscar.GAP

Dependency: Oscar.jl only. GAP is reached through Oscar's bundled GAP instance
(Oscar.GAP). The SL2Reps package must be installed in GAP (once per machine):

    julia> using Oscar
    julia> Oscar.GAP.Packages.install("SL2Reps")

The helper layer records the arithmetic data used by the pipeline:
- `matrix_conductor`: minimal cyclotomic conductor of matrix entries
- `gap_to_oscar_matrix`: conversion from GAP cyclotomics to Oscar matrices
- `compute_parity`: sign `ε` in `S² = ε·C`
- `exact_order`: order of the diagonal `T` matrix
"""

using Oscar
using LinearAlgebra: diag

# ============================================================
#  Atomic irrep data structure
# ============================================================

"""
    AtomicIrrep

A single irreducible representation of `SL(2, ℤ/N)`, realized as matrices
over the cyclotomic field `Q(ζ_N)`.

Atomic irreps are the representation-theoretic building blocks for strata.
The classification pipeline forms direct sums of these irreps, then searches
for modular-data bases inside the resulting representation space.

Fields:
- dim:    dimension of the irrep (= number of simple objects it provides)
- level:  order of the T matrix (divides N; the minimal conductor for this irrep)
- label:  human-readable identifier, e.g. "3d_16" = 3-dimensional, T of order 16
- S:      matrix representing ρ(𝔰) over K = Q(ζ_N)
- T:      matrix representing ρ(𝔱) over K = Q(ζ_N), diagonal
- parity: sign ε ∈ {+1, -1} such that S² = ε · C (C permutation matrix)
- K:      Oscar cyclotomic field Q(ζ_N)
- N:      the conductor of the embedding field (may be > level)
"""
struct AtomicIrrep
    dim::Int
    level::Int
    label::String
    S::Any   # Oscar matrix over K
    T::Any   # Oscar matrix over K (diagonal)
    parity::Int
    K::Any   # Oscar.CyclotomicField
    N::Int
end

function Base.show(io::IO, atom::AtomicIrrep)
    print(io, "AtomicIrrep($(atom.label), dim=$(atom.dim), level=$(atom.level), parity=$(atom.parity))")
end

# ============================================================
#  Helper functions
# ============================================================

"""
    all_divisors(N::Int) -> Vector{Int}

All positive divisors of N in increasing order.
"""
function all_divisors(N::Int)
    return [d for d in 1:N if N % d == 0]
end

"""
    matrix_conductor(M_gap, r::Int) -> Int

Compute the minimal conductor N_min such that all entries of the r×r GAP
matrix M lie in Q(ζ_{N_min}).

Uses GAP's `Conductor` (via Oscar.GAP) on each entry.
"""
function matrix_conductor(M_gap, r::Int)
    cond = 1
    for i in 1:r, j in 1:r
        entry = M_gap[i, j]
        c = Int(Oscar.GAP.Globals.Conductor(entry))
        cond = lcm(cond, c)
    end
    return cond
end

"""
    cyclotomic_to_oscar(x_gap, K, N::Int) -> elt of K

Convert a single GAP cyclotomic number (living in Q(ζ_M) for some M | N)
to an element of the Oscar cyclotomic field K = Q(ζ_N).

GAP represents cyclotomics as CoeffsCyc(x, N) = coefficients vector
[c_0, ..., c_{N-1}] such that x = Σ c_i · ζ_N^i.
"""
function cyclotomic_to_oscar(x_gap, K, N::Int)
    # Get coefficient vector of x with respect to ζ_N
    coeffs_gap = Oscar.GAP.Globals.CoeffsCyc(x_gap, N)
    z = gen(K)  # this is ζ_N in Oscar
    result = zero(K)
    for i in 1:N
        c = Rational{BigInt}(coeffs_gap[i])
        if c != 0
            result += c * z^(i - 1)
        end
    end
    return result
end

"""
    gap_to_oscar_matrix(M_gap, r::Int, K, N::Int) -> Oscar matrix

Convert an r×r GAP cyclotomic matrix to an Oscar matrix over K = Q(ζ_N).
"""
function gap_to_oscar_matrix(M_gap, r::Int, K, N::Int)
    M = zero_matrix(K, r, r)
    for i in 1:r, j in 1:r
        M[i, j] = cyclotomic_to_oscar(M_gap[i, j], K, N)
    end
    return M
end

"""
    compute_parity(S) -> Int

Compute ε ∈ {+1, -1} such that S² = ε · C for some permutation matrix C.
The parity distinguishes two classes of SL(2, ℤ/N) irreducibles.

Returns +1 if S² is close to a permutation matrix with +1 entries,
-1 if entries are -1. Throws if neither.
"""
function compute_parity(S)
    r = nrows(S)
    S2 = S * S
    K = base_ring(S)
    # Each row should have exactly one non-zero entry, value ±1
    first_sign = 0
    for i in 1:r
        non_zero_idx = 0
        non_zero_val = zero(K)
        for j in 1:r
            if !iszero(S2[i, j])
                if non_zero_idx != 0
                    error("S² row $i has multiple non-zero entries; not a (signed) permutation")
                end
                non_zero_idx = j
                non_zero_val = S2[i, j]
            end
        end
        non_zero_idx == 0 && error("S² row $i is all zero")
        sign = if non_zero_val == one(K)
            1
        elseif non_zero_val == -one(K)
            -1
        else
            error("S² entry ($i, $non_zero_idx) = $non_zero_val; not ±1")
        end
        first_sign == 0 && (first_sign = sign)
        sign == first_sign || error("S² has mixed signs")
    end
    return first_sign
end

"""
    exact_order(T, r::Int, N::Int) -> Int

Compute the order of the diagonal matrix T (= lcm of orders of its diagonal
entries as N-th roots of unity).
"""
function exact_order(T, r::Int, N::Int)
    K = base_ring(T)
    z = gen(K)
    ord = 1
    for i in 1:r
        θ = T[i, i]
        k = findfirst(j -> θ == z^(j - 1), 1:N)
        k === nothing && error("T[$i,$i] is not a power of ζ_$N")
        k = k - 1
        # Order of ζ_N^k in μ_N is N / gcd(k, N)
        d = gcd(k, N)
        entry_ord = N ÷ d
        # When k = 0, this gives ord = 1 (θ = 1), correct
        if d == N
            entry_ord = 1
        end
        ord = lcm(ord, entry_ord)
    end
    return ord
end

# ============================================================
#  Atomic catalog construction
# ============================================================

"""
    load_sl2reps_package()

Ensure the GAP SL2Reps package is loaded in Oscar.GAP. Throws informative
error if not available.
"""
function load_sl2reps_package()
    ok = Oscar.GAP.Globals.LoadPackage(Oscar.GAP.GapObj("SL2Reps"))
    # GAP returns `true` (a GAP boolean), Oscar.GAP converts to Julia Bool automatically
    ok == true || error(
        "Failed to load GAP SL2Reps package. Install once via:\n" *
        "    julia> using Oscar\n" *
        "    julia> Oscar.GAP.Packages.install(\"SL2Reps\")"
    )
    return nothing
end

"""
    build_atomic_catalog(N::Int; max_rank::Int = 20, verbose::Bool = true) -> Vector{AtomicIrrep}

Build a catalog of `SL(2, ℤ/N)` irreducible representations of dimension
`≤ max_rank`, realized as Oscar matrices over `Q(ζ_N)`.

Enumerates irreps at all levels (divisors of N), requiring that both
cond(S) | N and cond(T) | N so the irrep genuinely lives in Q(ζ_N).
Deduplicates by T-spectrum + dimension.
"""
function build_atomic_catalog(N::Int; max_rank::Int = 20, verbose::Bool = true)
    load_sl2reps_package()

    K, z = cyclotomic_field(N)

    catalog = AtomicIrrep[]
    seen = Set{UInt64}()

    for lev in all_divisors(N)
        reps = Oscar.GAP.evalstr("SL2IrrepsOfLevel($lev)")
        n_reps = Int(Oscar.GAP.Globals.Length(reps))
        for i in 1:n_reps
            rep = reps[i]
            r = Int(rep.degree)
            r > max_rank && continue

            cond_S = matrix_conductor(rep.S, r)
            cond_T = matrix_conductor(rep.T, r)
            (N % cond_S != 0 || N % cond_T != 0) && continue

            S = gap_to_oscar_matrix(rep.S, r, K, N)
            T = gap_to_oscar_matrix(rep.T, r, K, N)

            # Dedup signature: (dim, sorted exact T eigenvalues)
            t_entries = sort([string(T[k, k]) for k in 1:r])
            sig = hash((r, t_entries))
            sig in seen && continue
            push!(seen, sig)

            par = compute_parity(S)
            level = exact_order(T, r, N)
            label = "$(r)d_$level"
            push!(catalog, AtomicIrrep(r, level, label, S, T, par, K, N))
        end
    end

    if verbose
        println("Atomic catalog for N=$N: $(length(catalog)) irreps")
        for (i, a) in enumerate(catalog)
            println("  [$i] $a")
        end
    end
    return catalog
end
