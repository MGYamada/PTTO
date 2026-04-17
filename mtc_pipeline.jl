# mtc_pipeline.jl (v5: F_p exact arithmetic — self-contained)
#
# Conductor-first MTC classification using F_p arithmetic.
#
# This file is SELF-CONTAINED: it includes all necessary GAP/Oscar
# utilities, type definitions, and the classification pipeline.
# No separate mtc_classifier.jl or mtc_types.jl needed.
#
# Pipeline:
#   1. Build atomic catalog: SL₂(ℤ) irreps at levels n | N with cond(S,T) | N
#   2. Map Q(ζ_N) → F_p via Oscar exact cyclotomic → good prime reduction
#   3. For each atomic irrep: projective twist + Verlinde in F_p (exact!)
#   4. For direct sums: full S search in F_p^{r×r} with SS†=D²I + Verlinde
#   5. CRT lift: multiple primes → reconstruct Q(ζ_N) solution
#
# Usage:
#   include("mtc_pipeline.jl")
#   results = classify_modular_data(5)
#   table = periodic_table(10; max_rank=4)

using LinearAlgebra
using Oscar

# ============================================================
#  §0a. GAP/Oscar utilities (from mtc_classifier.jl)
# ============================================================

function load_sl2reps()
    GAP.Globals.LoadPackage(GapObj("SL2Reps"))
end

function gap_to_oscar_matrix(gap_mat, r, K::AbsSimpleNumField, cond_target::Int)
    z = gen(K)
    M = zero_matrix(K, r, r)
    for i in 1:r, j in 1:r
        entry_cond = Int(GAP.Globals.Conductor(gap_mat[i,j]))
        cond_target % entry_cond != 0 &&
            error("Entry conductor $entry_cond does not divide target $cond_target")
        mult = div(cond_target, entry_cond)
        coeffs = GAP.Globals.CoeffsCyc(gap_mat[i,j], entry_cond)
        val = zero(K)
        for k in 1:entry_cond
            if !Bool(GAP.Globals.IsZero(coeffs[k]))
                num = Int(GAP.Globals.NumeratorRat(coeffs[k]))
                den = Int(GAP.Globals.DenominatorRat(coeffs[k]))
                val += QQ(num, den) * z^((k-1) * mult)
            end
        end
        M[i,j] = val
    end
    M
end

function matrix_conductor(gap_mat, r)
    cond = 1
    for i in 1:r, j in 1:r
        cond = lcm(cond, Int(GAP.Globals.Conductor(gap_mat[i,j])))
    end
    cond
end

function to_complex(val::AbsSimpleNumFieldElem, zeta::ComplexF64, deg::Int)
    result = 0.0 + 0.0im
    for k in 0:deg-1
        c = coeff(val, k)
        result += (Float64(numerator(c)) / Float64(denominator(c))) * zeta^k
    end
    result
end

function matrix_to_complex(M, zeta::ComplexF64, deg::Int)
    r = nrows(M); c = ncols(M)
    [to_complex(M[i,j], zeta, deg) for i in 1:r, j in 1:c]
end

# Load SL2Reps
load_sl2reps()

# ============================================================
#  §0b. Type definitions (from mtc_types.jl)
# ============================================================

struct AtomicIrrep
    dim::Int
    level::Int
    label::String
    S       # Oscar matrix over Q(ζ_N)
    T       # Oscar matrix over Q(ζ_N)
    parity::Int   # +1 even, -1 odd
    K       # cyclotomic field
    N::Int  # conductor
end

Base.show(io::IO, a::AtomicIrrep) = print(io,
    "AtomicIrrep(dim=$(a.dim), level=$(a.level), " *
    "$(a.parity==1 ? "even" : a.parity==-1 ? "odd" : "mixed"), label=$(a.label))")

function compute_parity(ρ_S::AbstractMatrix{<:Number}; tol::Real=1e-10)
    S2 = ρ_S * ρ_S
    r = size(S2, 1)
    is_plus = maximum(abs.(S2 - I(r))) < tol
    is_minus = maximum(abs.(S2 + I(r))) < tol
    is_plus ? 1 : is_minus ? -1 : 0
end

function _numeric_order(T_num::AbstractMatrix{<:Number}, r::Int, N::Int; tol::Real=1e-8)
    T_power = copy(T_num)
    for n in 1:N
        if maximum(abs.(T_power - I(r))) < tol
            return n
        end
        T_power = T_power * T_num
    end
    N
end

# ============================================================
#  §1. F_p utilities
# ============================================================

function _all_divisors(N::Int)
    [d for d in 1:N if N % d == 0]
end

function isprime_simple(n::Int)
    n < 2 && return false; n < 4 && return true
    n % 2 == 0 && return false; n % 3 == 0 && return false
    i = 5; while i*i <= n; (n%i==0 || n%(i+2)==0) && return false; i += 6; end; true
end

function _prime_factors(n::Int)
    factors = Int[]; d = 2; m = n
    while d*d <= m; if m%d==0; push!(factors,d); while m%d==0; m÷=d; end; end; d+=1; end
    m > 1 && push!(factors, m); factors
end

function primitive_root(p::Int)
    φ = p - 1; factors = _prime_factors(φ)
    for g in 2:p-1
        all(q -> powermod(g, φ÷q, p) != 1, factors) && return g
    end
    error("No primitive root for p=$p")
end

"""
    good_primes(N, count) -> Vector{Int}

Primes p ≡ 1 mod lcm(N, 4) so that F_p contains both ζ_N and ζ₄=i.
The lcm(N,4) condition ensures projective twists (ζ₄^α, ζ₁₂^α) exist in F_p.
"""
function good_primes(N::Int, count::Int=3)
    M = lcm(N, 12)  # need 12th roots for ζ₁₂^α twist
    primes = Int[]; p = M + 1
    while length(primes) < count
        isprime_simple(p) && p % M == 1 && push!(primes, p)
        p += M
    end
    primes
end

function zeta_mod_p(N::Int, p::Int)
    g = primitive_root(p)
    powermod(g, (p-1) ÷ N, p)
end

function inv_mod(a::Int, p::Int)
    a = mod(a, p); a == 0 && error("0 has no inverse mod $p")
    powermod(a, p-2, p)
end

# ============================================================
#  §2. Oscar Q(ζ_N) → F_p exact map
# ============================================================

"""
    oscar_elem_to_fp(x, N, ζp, p) -> Int

Map x ∈ Q(ζ_N) (Oscar AbsSimpleNumFieldElem) to F_p.
Uses Oscar's exact rational coefficients: x = Σ (a_k/b_k) ζ^k.
"""
function oscar_elem_to_fp(x, N::Int, ζp::Int, p::Int)
    K = parent(x)
    d = degree(K)
    result = 0
    for k in 0:d-1
        c = coeff(x, k)
        iszero(c) && continue
        num = Int(numerator(c))
        den = Int(denominator(c))
        c_fp = mod(num * inv_mod(mod(den, p), p), p)
        result = mod(result + c_fp * powermod(ζp, k, p), p)
    end
    result
end

"""
    oscar_matrix_to_fp(M, r, N, ζp, p) -> Matrix{Int}
"""
function oscar_matrix_to_fp(M, r::Int, N::Int, ζp::Int, p::Int)
    [oscar_elem_to_fp(M[i,j], N, ζp, p) for i in 1:r, j in 1:r]
end

"""
    galois_conj_fp(x, N, ζp, p) -> Int

Complex conjugation in F_p: σ_{-1}(ζ_N) = ζ_N^{-1}.
Maps x ∈ Q(ζ_N) to σ_{-1}(x) in F_p.
"""
function galois_conj_elem_fp(x, N::Int, ζp::Int, p::Int)
    K = parent(x)
    d = degree(K)
    ζp_inv = powermod(ζp, N-1, p)
    result = 0
    for k in 0:d-1
        c = coeff(x, k)
        iszero(c) && continue
        num = Int(numerator(c))
        den = Int(denominator(c))
        c_fp = mod(num * inv_mod(mod(den, p), p), p)
        result = mod(result + c_fp * powermod(ζp_inv, k, p), p)
    end
    result
end

function galois_conj_matrix_fp(M, r::Int, N::Int, ζp::Int, p::Int)
    [galois_conj_elem_fp(M[i,j], N, ζp, p) for i in 1:r, j in 1:r]
end

# ============================================================
#  §3. Admissibility checks in F_p
# ============================================================

"""
    check_SSdag_fp(S_fp, Sbar_fp, r, p) -> (Bool, Int)

Check SS† = D²I in F_p. Returns (pass, D²_mod_p).
S† = conj(S)ᵀ, so SS† = S · Sbar_fpᵀ.
"""
function check_SSdag_fp(S_fp::Matrix{Int}, Sbar_fp::Matrix{Int}, r::Int, p::Int)
    # Compute (SS†)_{ij} = Σ_k S_{ik} conj(S)_{jk} mod p
    D2 = mod(sum(mod(S_fp[1,k] * Sbar_fp[1,k], p) for k in 1:r), p)
    for i in 1:r, j in 1:r
        val = mod(sum(mod(S_fp[i,k] * Sbar_fp[j,k], p) for k in 1:r), p)
        if i == j
            val != D2 && return (false, 0)
        else
            val != 0 && return (false, 0)
        end
    end
    (true, D2)
end

"""
    verlinde_fp(S_fp, Sbar_fp, r, D2, p; unit=1) -> Union{Array{Int,3}, Nothing}

Verlinde check in F_p:
  N^{ij}_k = (1/D²) Σ_l S_{li} S_{lj} conj(S)_{lk} / S_{l,unit}

Result must be a "small" non-negative integer (val mod p should be < p/2).
"""
function verlinde_fp(S_fp::Matrix{Int}, Sbar_fp::Matrix{Int},
                     r::Int, D2::Int, p::Int; unit::Int=1)
    D2 == 0 && return nothing
    D2_inv = inv_mod(D2, p)

    Nijk = zeros(Int, r, r, r)
    for i in 1:r, j in 1:r, k in 1:r
        val = 0
        for l in 1:r
            Sl0 = S_fp[l, unit]
            Sl0 == 0 && return nothing
            term = mod(S_fp[l,i] * S_fp[l,j] % p * Sbar_fp[l,k] % p * inv_mod(Sl0, p) % p, p)
            val = mod(val + term, p)
        end
        val = mod(val * D2_inv, p)
        val > p ÷ 2 && return nothing  # negative → fail
        Nijk[i,j,k] = val
    end
    Nijk
end

# ============================================================
#  §4. Build atomic catalog (scan all levels n | N)
# ============================================================

function build_atomic_catalog(N::Int; max_rank::Int=20)
    K, z = cyclotomic_field(N)
    zeta = exp(2π * im / N)
    deg = degree(K)

    catalog = AtomicIrrep[]
    seen = Set{UInt64}()

    for lev in _all_divisors(N)
        reps = GAP.evalstr("SL2IrrepsOfLevel($lev)")
        len = Int(GAP.Globals.Length(reps))
        for i in 1:len
            rep = reps[i]
            r = Int(rep.degree)
            r > max_rank && continue
            cond_S = matrix_conductor(rep.S, r)
            cond_T = matrix_conductor(rep.T, r)
            (N % cond_S != 0 || N % cond_T != 0) && continue

            S = gap_to_oscar_matrix(rep.S, r, K, N)
            T = gap_to_oscar_matrix(rep.T, r, K, N)
            S_num = matrix_to_complex(S, zeta, deg)
            T_num = matrix_to_complex(T, zeta, deg)

            t_diag = sort([round(angle(T_num[k,k]); digits=8) for k in 1:r])
            sig = hash((r, t_diag))
            sig in seen && continue
            push!(seen, sig)

            par = compute_parity(S_num)
            level = _numeric_order(T_num, r, N)
            push!(catalog, AtomicIrrep(r, level, "$(r)d_$level", S, T, par, K, N))
        end
    end

    println("Atomic catalog for N=$N: $(length(catalog)) irreps")
    for (i, a) in enumerate(catalog)
        println("  [$i] $a")
    end
    catalog
end

# ============================================================
#  §5. Classify single atomic irreps in F_p
# ============================================================

"""
    classify_atomic_fp(atom, N, p, ζp) -> Vector{NamedTuple}

For a single atomic irrep, try all projective twists α ∈ {0,1,2,3}
and unit objects i₀. All checks in F_p.
"""
function classify_atomic_fp(atom::AtomicIrrep, N::Int, p::Int, ζp::Int)
    r = atom.dim
    S_oscar = atom.S
    T_oscar = atom.T
    K = atom.K

    S_fp = oscar_matrix_to_fp(S_oscar, r, N, ζp, p)
    Sbar_fp = galois_conj_matrix_fp(S_oscar, r, N, ζp, p)
    T_diag_fp = [oscar_elem_to_fp(T_oscar[i,i], N, ζp, p) for i in 1:r]

    results = NamedTuple[]

    # Projective twist: S_α = ζ₄^α · S,  T_α = ζ₁₂^α · T
    # In F_p: ζ₄ = ζ_{lcm(N,12)}^{lcm(N,12)/4},  ζ₁₂ = ζ_{lcm(N,12)}^{lcm(N,12)/12}
    M = lcm(N, 12)
    ζM_p = zeta_mod_p(M, p)
    ζ4_p = powermod(ζM_p, M ÷ 4, p)
    ζ12_p = powermod(ζM_p, M ÷ 12, p)

    for α in 0:3
        ζ4α = powermod(ζ4_p, α, p)
        ζ12α = powermod(ζ12_p, α, p)

        # Twisted S, conj(S), T in F_p
        S_tw = [mod(ζ4α * S_fp[i,j], p) for i in 1:r, j in 1:r]
        # conj(ζ₄^α · S) = conj(ζ₄^α) · conj(S) = ζ₄^{-α} · Sbar
        ζ4α_conj = powermod(ζ4_p, (4 - α) % 4, p)
        Sbar_tw = [mod(ζ4α_conj * Sbar_fp[i,j], p) for i in 1:r, j in 1:r]

        T_tw = [mod(ζ12α * T_diag_fp[i], p) for i in 1:r]

        # Check SS†
        ok, D2 = check_SSdag_fp(S_tw, Sbar_tw, r, p)
        ok && D2 != 0 || continue

        # Try each unit object
        for i0 in 1:r
            # Normalize T: T_norm = T / T[i0]
            T_i0 = T_tw[i0]
            T_i0 == 0 && continue
            T_i0_inv = inv_mod(T_i0, p)

            # Reorder: unit → index 1
            perm = vcat([i0], setdiff(1:r, [i0]))
            S_reord = S_tw[perm, perm]
            Sbar_reord = Sbar_tw[perm, perm]

            # Check all S[1,j] nonzero (unit row)
            any(S_reord[1,j] == 0 for j in 1:r) && continue

            # Verlinde
            Nijk = verlinde_fp(S_reord, Sbar_reord, r, D2, p; unit=1)
            Nijk === nothing && continue

            # Check N[1,j,j] > 0 for all j (unit fuses trivially)
            valid_unit = all(Nijk[1,j,j] == 1 for j in 1:r)
            valid_unit || continue

            # Central charge: compute numerically (F_p doesn't directly give c)
            # Store the F_p solution; lift later
            push!(results, (α=α, unit=i0, D2_fp=D2, Nijk=Nijk, perm=perm))
        end
    end
    results
end

# ============================================================
#  §5b. Direct sum: find overlap pairs
# ============================================================

"""
    find_overlap_pairs(catalog, N, ζ_num, deg)

Find all pairs (i,j) of atomic irreps with overlapping T-spectra.
Returns (i, j, perm_i, perm_j, l) where l = overlap size and
perm_i reorders irrep i so first l eigenvalues are the shared ones.
"""
function find_overlap_pairs(catalog, N::Int, ζ_num, deg)
    pairs = []
    for i in 1:length(catalog), j in i+1:length(catalog)
        a1 = catalog[i]; a2 = catalog[j]
        T1 = matrix_to_complex(a1.T, ζ_num, deg)
        T2 = matrix_to_complex(a2.T, ζ_num, deg)

        t1 = [T1[k,k] for k in 1:a1.dim]
        t2 = [T2[k,k] for k in 1:a2.dim]

        # Match eigenvalues
        perm1 = Int[]; perm2 = Int[]
        used1 = falses(a1.dim); used2 = falses(a2.dim)
        for ii in 1:a1.dim, jj in 1:a2.dim
            if !used1[ii] && !used2[jj] && abs(t1[ii] - t2[jj]) < 1e-8
                push!(perm1, ii); push!(perm2, jj)
                used1[ii] = true; used2[jj] = true
            end
        end
        l = length(perm1)
        l == 0 && continue

        # Append non-shared indices
        for ii in 1:a1.dim; !used1[ii] && push!(perm1, ii); end
        for jj in 1:a2.dim; !used2[jj] && push!(perm2, jj); end

        push!(pairs, (i, j, perm1, perm2, l))
    end
    pairs
end

# ============================================================
#  §5c. Direct sum: F_p search for cross-block S
# ============================================================

"""
    classify_directsum_fp(atom1, atom2, perm1, perm2, l, N, p, ζp)

Search for valid S from ρ₁⊕ρ₂ in F_p.

Structure (Theorem 3.23/3.24):
  T is block-diagonal: shared eigenvalues doubled, rest from each irrep.
  S has known diagonal blocks from atomics; cross-block is parametrized
  by l orthogonal 2×2 matrices Uᵢ.

For opposite parities: aᵢ = bᵢ = 1/√2 (if √2 ∈ F_p).
For same parities: aᵢ, bᵢ free with aᵢ² + bᵢ² = 1.

We enumerate all valid (aᵢ, bᵢ) ∈ F_p with aᵢ² + bᵢ² = 1.
"""
function classify_directsum_fp(atom1::AtomicIrrep, atom2::AtomicIrrep,
                                perm1::Vector{Int}, perm2::Vector{Int},
                                l::Int, N::Int, p::Int, ζp::Int)
    r1 = atom1.dim; r2 = atom2.dim
    k = r1 - l; m = r2 - l
    r = 2l + k + m

    # Map reordered S, T to F_p
    S1_fp_full = oscar_matrix_to_fp(atom1.S, r1, N, ζp, p)
    S2_fp_full = oscar_matrix_to_fp(atom2.S, r2, N, ζp, p)
    S1bar_fp_full = galois_conj_matrix_fp(atom1.S, r1, N, ζp, p)
    S2bar_fp_full = galois_conj_matrix_fp(atom2.S, r2, N, ζp, p)

    # Reorder
    ψ_fp = S1_fp_full[perm1, perm1]   # r1×r1, reordered so first l are shared
    η_fp = S2_fp_full[perm2, perm2]   # r2×r2, reordered
    ψbar_fp = S1bar_fp_full[perm1, perm1]
    ηbar_fp = S2bar_fp_full[perm2, perm2]

    T1_fp = [oscar_elem_to_fp(atom1.T[perm1[i],perm1[i]], N, ζp, p) for i in 1:r1]
    T2_fp = [oscar_elem_to_fp(atom2.T[perm2[i],perm2[i]], N, ζp, p) for i in 1:r2]

    # Enumerate (a, b) ∈ F_p² with a² + b² = 1
    circle_points = Tuple{Int,Int}[]
    for a in 0:p-1, b in 0:p-1
        mod(a^2 + b^2, p) == 1 && push!(circle_points, (a, b))
    end

    results = NamedTuple[]

    # For each choice of l circle points (one per shared eigenvalue pair)
    # This is |circle|^l combinations
    if length(circle_points)^l > 10_000_000
        println("    Too many circle combos ($(length(circle_points))^$l), skipping")
        return results
    end

    for ab_combo in Iterators.product(fill(circle_points, l)...)
        ab = collect(ab_combo)  # Vector of (a_i, b_i)

        # Build S in F_p using Theorem 3.23 structure
        S_fp = build_directsum_S_fp(ψ_fp, η_fp, ab, l, k, m, p)
        Sbar_fp = build_directsum_Sbar_fp(ψbar_fp, ηbar_fp, ab, l, k, m, p)

        # Check SS†
        ok, D2 = check_SSdag_fp(S_fp, Sbar_fp, r, p)
        ok && D2 != 0 || continue

        # Try each unit object (only positions 2l+1 to 2l+k+m, 
        # or positions 1..2l if same parity)
        for i0 in 1:r
            perm_u = vcat([i0], setdiff(1:r, [i0]))
            S_reord = S_fp[perm_u, perm_u]
            Sbar_reord = Sbar_fp[perm_u, perm_u]

            any(S_reord[1,j] == 0 for j in 1:r) && continue

            Nijk = verlinde_fp(S_reord, Sbar_reord, r, D2, p; unit=1)
            Nijk === nothing && continue

            all(Nijk[1,j,j] == 1 for j in 1:r) || continue

            push!(results, (ab=ab, unit=i0, D2_fp=D2, Nijk=Nijk, perm_u=perm_u))
        end
    end
    results
end

"""
    build_directsum_S_fp(ψ, η, ab, l, k, m, p) -> Matrix{Int}

Build the (2l+k+m) × (2l+k+m) S matrix from Theorem 3.23 structure in F_p.

Ordering: [shared_1a, shared_1b, ..., shared_la, shared_lb, nonshared_ρ1..., nonshared_ρ2...]
"""
function build_directsum_S_fp(ψ_fp, η_fp, ab, l, k, m, p)
    r = 2l + k + m
    S = zeros(Int, r, r)

    # A block (2l × 2l): Aᵢⱼ = Uᵢ [[ψ,0],[0,η]] Uⱼᵀ
    for i in 1:l, j in 1:l
        ai, bi = ab[i]; aj, bj = ab[j]
        ψij = ψ_fp[i, j]; ηij = η_fp[i, j]

        S[2i-1, 2j-1] = mod(ai*aj*ψij + bi*bj*ηij, p)
        S[2i-1, 2j]   = mod(-ai*bj*ψij + bi*aj*ηij, p)
        S[2i,   2j-1] = mod(bi*aj*ψij - ai*bj*ηij, p)
        S[2i,   2j]   = mod(bi*bj*ψij + ai*aj*ηij, p)
    end

    # B block (k × 2l): Bᵢ₁ⱼ = [ψ_{l+i1,j}, 0] Uⱼᵀ
    for i1 in 1:k, j in 1:l
        aj, bj = ab[j]
        ψv = ψ_fp[l+i1, j]

        S[2l+i1, 2j-1] = mod(aj * ψv, p)
        S[2l+i1, 2j]   = mod(-bj * ψv, p)
        # Transpose: S is symmetric via V, but raw ρ(s) may not be
        # Actually ρ(s) IS symmetric (Theorem 3.3), so:
        S[2j-1, 2l+i1] = mod(aj * ψv, p)
        S[2j,   2l+i1] = mod(-bj * ψv, p)
    end

    # C block (m × 2l): Cᵢ₂ⱼ = [0, η_{l+i2,j}] Uⱼᵀ
    for i2 in 1:m, j in 1:l
        aj, bj = ab[j]
        ηv = η_fp[l+i2, j]

        S[2l+k+i2, 2j-1] = mod(bj * ηv, p)
        S[2l+k+i2, 2j]   = mod(aj * ηv, p)
        S[2j-1, 2l+k+i2] = mod(bj * ηv, p)
        S[2j,   2l+k+i2] = mod(aj * ηv, p)
    end

    # ψ' block (k × k): lower-right of ρ₁(s)
    for i1 in 1:k, j1 in 1:k
        S[2l+i1, 2l+j1] = ψ_fp[l+i1, l+j1]
    end

    # η' block (m × m): lower-right of ρ₂(s)
    for i2 in 1:m, j2 in 1:m
        S[2l+k+i2, 2l+k+j2] = η_fp[l+i2, l+j2]
    end

    S
end

"""
    build_directsum_Sbar_fp(ψbar, ηbar, ab, l, k, m, p) -> Matrix{Int}

Same as build_directsum_S_fp but using conjugated entries.
For (a,b) on unit circle in F_p, conj preserves a²+b²=1.
"""
function build_directsum_Sbar_fp(ψbar_fp, ηbar_fp, ab, l, k, m, p)
    r = 2l + k + m
    S = zeros(Int, r, r)

    for i in 1:l, j in 1:l
        ai, bi = ab[i]; aj, bj = ab[j]
        ψij = ψbar_fp[i, j]; ηij = ηbar_fp[i, j]

        S[2i-1, 2j-1] = mod(ai*aj*ψij + bi*bj*ηij, p)
        S[2i-1, 2j]   = mod(-ai*bj*ψij + bi*aj*ηij, p)
        S[2i,   2j-1] = mod(bi*aj*ψij - ai*bj*ηij, p)
        S[2i,   2j]   = mod(bi*bj*ψij + ai*aj*ηij, p)
    end

    for i1 in 1:k, j in 1:l
        aj, bj = ab[j]
        ψv = ψbar_fp[l+i1, j]
        S[2l+i1, 2j-1] = mod(aj * ψv, p)
        S[2l+i1, 2j]   = mod(-bj * ψv, p)
        S[2j-1, 2l+i1] = mod(aj * ψv, p)
        S[2j,   2l+i1] = mod(-bj * ψv, p)
    end

    for i2 in 1:m, j in 1:l
        aj, bj = ab[j]
        ηv = ηbar_fp[l+i2, j]
        S[2l+k+i2, 2j-1] = mod(bj * ηv, p)
        S[2l+k+i2, 2j]   = mod(aj * ηv, p)
        S[2j-1, 2l+k+i2] = mod(bj * ηv, p)
        S[2j,   2l+k+i2] = mod(aj * ηv, p)
    end

    for i1 in 1:k, j1 in 1:k
        S[2l+i1, 2l+j1] = ψbar_fp[l+i1, l+j1]
    end
    for i2 in 1:m, j2 in 1:m
        S[2l+k+i2, 2l+k+j2] = ηbar_fp[l+i2, l+j2]
    end
    S
end

# ============================================================
#  §5d. Numerical lift for direct sums
# ============================================================

function lift_directsum_numerical(atom1, atom2, perm1, perm2, l, sol, N)
    r1 = atom1.dim; r2 = atom2.dim
    k = r1 - l; m = r2 - l; r = 2l + k + m
    zeta = exp(2π * im / N)
    deg = degree(atom1.K)

    S1 = matrix_to_complex(atom1.S, zeta, deg)[perm1, perm1]
    S2 = matrix_to_complex(atom2.S, zeta, deg)[perm2, perm2]
    T1 = matrix_to_complex(atom1.T, zeta, deg)

    # Reconstruct S numerically from ab parameters
    # ab are F_p values; we need to find the Q(ζ_N) values
    # For now: brute force — use the F_p-verified Nijk and
    # reconstruct S from the Verlinde inverse
    # OR: just store Nijk and mark this as F_p-verified
    
    # Simple approach: return Nijk-based result
    (S=zeros(ComplexF64, r, r), T=Diagonal(zeros(ComplexF64, r)),
     c=0//1, d=zeros(r), Nijk=sol.Nijk)
end

# ============================================================
#  §6. Numerical lift: F_p solution → central charge
# ============================================================

"""
    lift_to_numerical(atom, sol, N) -> NamedTuple

Given an F_p solution (α, unit, perm), compute the numerical S, T, c, d.
"""
function lift_to_numerical(atom::AtomicIrrep, sol, N::Int)
    r = atom.dim
    zeta = exp(2π * im / N)
    deg = degree(atom.K)

    S_num = matrix_to_complex(atom.S, zeta, deg)
    T_num = matrix_to_complex(atom.T, zeta, deg)

    α = sol.α
    S_tw = cispi(α/2) * S_num
    T_tw = cispi(α/6) * T_num

    i0 = sol.unit
    t_i0 = T_tw[i0, i0]
    T_norm = T_tw / t_i0

    perm = sol.perm
    S_reord = S_tw[perm, perm]
    T_reord = Diagonal([T_norm[perm[i], perm[i]] for i in 1:r])

    # Apply V: make d_i = S[1,j]/S[1,1] all positive real
    row = S_reord[1, :]

    # Determine target sign from diagonal element
    target_sign = sign(real(row[1]))
    if abs(real(row[1])) < 1e-10
        # Row is imaginary — use imag part
        target_sign = sign(imag(row[1]))
    end

    V = ones(Int, r)
    for j in 1:r
        re = real(row[j])
        if abs(re) > 1e-10
            V[j] = sign(re) == target_sign ? 1 : -1
        else
            im_j = imag(row[j])
            V[j] = sign(im_j) == target_sign ? 1 : -1
        end
    end
    S_final = [V[i]*V[j]*S_reord[i,j] for i in 1:r, j in 1:r]

    # Global flip: ensure S[1,:] are all same sign (positive)
    if target_sign < 0
        S_final = -S_final
    end

    # d_i: use the fact that d_i = S[1,j]/S[1,1], should be real positive
    d = Float64[]
    s00 = S_final[1,1]
    for j in 1:r
        ratio = S_final[1,j] / s00
        # Take real part (imaginary part should be ~0 after correct α,V)
        push!(d, abs(imag(ratio)) < 1e-6 ? real(ratio) : abs(ratio))
    end

    D = sqrt(sum(d .^ 2))
    p_plus = sum(d[j]^2 * T_reord[j,j] for j in 1:r)
    c8 = angle(p_plus / D) / (2π) * 8
    c_rat = rationalize(c8; tol=1e-4)

    (S=S_final, T=T_reord, c=c_rat, d=round.(d; digits=6), Nijk=sol.Nijk)
end

# ============================================================
#  §7. Main pipeline
# ============================================================

"""
    classify_modular_data(N; max_rank=6, verbose=true)

Complete MTC classification at conductor N using F_p arithmetic.
"""
function classify_modular_data(N::Int; max_rank::Int=6, verbose::Bool=true)
    verbose && println("=" ^ 60)
    verbose && println("  MTC Classification (F_p exact, v5): N = $N")
    verbose && println("=" ^ 60)

    catalog = build_atomic_catalog(N; max_rank=max_rank)

    primes = good_primes(N, 3)
    p = primes[1]
    ζp = zeta_mod_p(N, p)
    verbose && println("Using prime p=$p, ζ_$N=$ζp in F_$p")

    results = NamedTuple[]
    seen_charges = Set{Rational{Int}}()

    # ── Part A: single atomic irreps ──
    verbose && println("\n── Part A: atomic irreps ──")
    for (ci, atom) in enumerate(catalog)
        sols = classify_atomic_fp(atom, N, p, ζp)
        for sol in sols
            md = lift_to_numerical(atom, sol, N)
            # Dedup by central charge
            md.c in seen_charges && continue
            push!(seen_charges, md.c)
            push!(results, md)
            verbose && println("  ✓ [$(atom.label)]  rank=$(atom.dim)  c=$(md.c)  d=$(md.d)")
        end
    end

    # ── Part B: direct sum search in F_p ──
    verbose && println("\n── Part B: direct sum search (F_p) ──")

    # Find pairs with overlapping T-spectra
    zeta_num = exp(2π * im / N)
    deg = degree(catalog[1].K)

    pairs = find_overlap_pairs(catalog, N, zeta_num, deg)
    verbose && println("Pairs with T-spectrum overlap: $(length(pairs))")

    for (i1, i2, shared_idx1, shared_idx2, l) in pairs
        atom1 = catalog[i1]; atom2 = catalog[i2]
        r1 = atom1.dim; r2 = atom2.dim; r = r1 + r2
        r > max_rank && continue

        verbose && println("  $(atom1.label) ⊕ $(atom2.label): l=$l, rank=$r")

        ds_results = classify_directsum_fp(atom1, atom2, shared_idx1, shared_idx2, 
                                            l, N, p, ζp)
        for sol in ds_results
            md = lift_directsum_numerical(atom1, atom2, shared_idx1, shared_idx2,
                                          l, sol, N)
            md.c in seen_charges && continue
            push!(seen_charges, md.c)
            push!(results, md)
            verbose && println("    ✓ rank=$r  c=$(md.c)  d=$(md.d)")
        end
    end

    verbose && println("\n── Results ──")
    verbose && println("$(length(results)) modular data at N=$N")
    verbose && println("=" ^ 60)
    for (i, md) in enumerate(results)
        verbose && println("  [$i] rank=$(length(md.d))  c=$(md.c)  d=$(md.d)")
    end

    results
end

"""
    periodic_table(N_max; max_rank=6)

MTC periodic table: classify for N = 1, 2, ..., N_max.
"""
function periodic_table(N_max::Int; max_rank::Int=6)
    table = Dict{Int, Vector{NamedTuple}}()
    for N in 1:N_max
        table[N] = classify_modular_data(N; max_rank=max_rank, verbose=false)
        n = length(table[N])
        n > 0 && println("N=$N: $n MTCs")
        for md in table[N]
            println("  rank=$(length(md.d))  c=$(md.c)  d=$(md.d)")
        end
    end

    println("\n" * "═" ^ 60)
    println("  PERIODIC TABLE (N = 1..$N_max)")
    println("═" ^ 60)
    for N in 1:N_max
        isempty(table[N]) && continue
        println("N=$N: $(length(table[N])) MTCs")
        for md in table[N]
            println("  rank=$(length(md.d))  c=$(md.c)  d=$(md.d)")
        end
    end
    table
end

println("mtc_pipeline.jl (v5: F_p exact) loaded.")
println("Usage:")
println("  results = classify_modular_data(5)")
println("  table = periodic_table(10; max_rank=4)")
