"""
SL(2, ℤ/N) irreducible representation catalog (atomic irreps).

This module provides:
- `AtomicIrrep`: data structure for a single SL(2, ℤ/N) irrep realized over Q(ζ_N)
- `build_atomic_catalog`: enumerate all SL(2, ℤ/N) irreps of a given conductor N
  using GAP's SL2Reps package (accessed via Oscar.GAP)

Dependency: Oscar.jl only. GAP is reached through Oscar's bundled GAP instance
(Oscar.GAP). The SL2Reps package must be installed in GAP (once per machine):

    julia> using Oscar
    julia> Oscar.GAP.Packages.install("SL2Reps")

Based on v5's `build_atomic_catalog` (mtc_pipeline.jl). Helpers that were in
v5's unavailable `mtc_types.jl` are reimplemented here:
- matrix_conductor: min N such that entries of M lie in Q(ζ_N)
- gap_to_oscar_matrix: convert GAP cyclotomic matrix to Oscar matrix over Q(ζ_N)
- matrix_to_complex: numerical complex approximation of Oscar matrix
- compute_parity: extract sign ε in S² = ε·C (charge conjugation sign)
- numeric_order: order of matrix T (= conductor of its eigenvalues)
"""

using Oscar
using LinearAlgebra: diag

# ============================================================
#  Atomic irrep data structure
# ============================================================

"""
    AtomicIrrep

A single irreducible representation of SL(2, ℤ/N), realized as matrices
over the cyclotomic field Q(ζ_N).

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
#  Helper functions (reimplementation of v5's mtc_types.jl)
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
    matrix_to_complex(M_oscar, ζ::ComplexF64, deg::Int) -> Matrix{ComplexF64}

Convert an Oscar matrix over Q(ζ_N) to a numerical complex matrix by
substituting ζ_N = e^(2πi/N).

`deg` is φ(N), the dimension of Q(ζ_N) over Q.
"""
function matrix_to_complex(M_oscar, ζ::ComplexF64, deg::Int)
    r = nrows(M_oscar)
    r_c = ncols(M_oscar)
    result = zeros(ComplexF64, r, r_c)
    for i in 1:r, j in 1:r_c
        coeffs = Oscar.coordinates(M_oscar[i, j])
        # coordinates returns coefficients in integral basis {1, ζ, ζ², ..., ζ^(deg-1)}
        val = zero(ComplexF64)
        for k in 1:length(coeffs)
            val += Float64(Rational{BigInt}(coeffs[k])) * ζ^(k - 1)
        end
        result[i, j] = val
    end
    return result
end

"""
    compute_parity(S_num::Matrix{ComplexF64}) -> Int

Compute ε ∈ {+1, -1} such that S² = ε · C for some permutation matrix C.
The parity distinguishes two classes of SL(2, ℤ/N) irreducibles.

Returns +1 if S² is close to a permutation matrix with +1 entries,
-1 if entries are -1. Throws if neither.
"""
function compute_parity(S_num::Matrix{ComplexF64}; tol::Float64 = 1e-8)
    r = size(S_num, 1)
    S2 = S_num * S_num
    # Each row should have exactly one non-zero entry, value ±1
    for i in 1:r
        non_zero_idx = 0
        non_zero_val = 0.0 + 0.0im
        for j in 1:r
            if abs(S2[i, j]) > tol
                if non_zero_idx != 0
                    error("S² row $i has multiple non-zero entries; not a (signed) permutation")
                end
                non_zero_idx = j
                non_zero_val = S2[i, j]
            end
        end
        non_zero_idx == 0 && error("S² row $i is all zero")
        # Check value is ±1
        if abs(non_zero_val - 1) < tol
            # +1 entry
        elseif abs(non_zero_val + 1) < tol
            # -1 entry
        else
            error("S² entry ($i, $non_zero_idx) = $non_zero_val; not ±1")
        end
    end
    # All entries should have same sign
    first_entry = S2[1, findfirst(j -> abs(S2[1, j]) > tol, 1:r)]
    return real(first_entry) > 0 ? +1 : -1
end

"""
    numeric_order(T_num::Matrix{ComplexF64}, r::Int, N::Int) -> Int

Compute the order of the diagonal matrix T (= lcm of orders of its diagonal
entries as N-th roots of unity).
"""
function numeric_order(T_num::Matrix{ComplexF64}, r::Int, N::Int; tol::Float64 = 1e-8)
    ord = 1
    for i in 1:r
        θ = T_num[i, i]
        # θ is an N-th root of unity: θ = exp(2πi k / N) for some k
        angle_val = angle(θ)
        k_float = angle_val * N / (2π)
        k = round(Int, k_float)
        k = mod(k, N)
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
    ok == Oscar.GAP.Globals.true || error(
        "Failed to load GAP SL2Reps package. Install once via:\n" *
        "    julia> using Oscar\n" *
        "    julia> Oscar.GAP.Packages.install(\"SL2Reps\")"
    )
    return nothing
end

"""
    build_atomic_catalog(N::Int; max_rank::Int = 20, verbose::Bool = true) -> Vector{AtomicIrrep}

Build the catalog of all SL(2, ℤ/N) irreducible representations of dimension
≤ max_rank, realized as Oscar matrices over Q(ζ_N).

Enumerates irreps at all levels (divisors of N), requiring that both
cond(S) | N and cond(T) | N so the irrep genuinely lives in Q(ζ_N).
Deduplicates by T-spectrum + dimension.

Port of v5's build_atomic_catalog (mtc_pipeline.jl:192-232).
"""
function build_atomic_catalog(N::Int; max_rank::Int = 20, verbose::Bool = true)
    load_sl2reps_package()

    K, z = cyclotomic_field(N)
    ζ_num = exp(2π * im / N)
    deg = degree(K)

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
            S_num = matrix_to_complex(S, ζ_num, deg)
            T_num = matrix_to_complex(T, ζ_num, deg)

            # Dedup signature: (dim, sorted T-phases)
            t_phases = sort([round(angle(T_num[k, k]); digits = 8) for k in 1:r])
            sig = hash((r, t_phases))
            sig in seen && continue
            push!(seen, sig)

            par = compute_parity(S_num)
            level = numeric_order(T_num, r, N)
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
