# mtc_classifier.jl
# MTC Classification via Conductor Decomposition
# April 14-15, 2026 — Arithmetic Condensed Matter Geometry
#
# Level I (exact over Q(ζ_N)):
#   N = ∏ pᵢ^aᵢ → SL₂(Z/pᵢ^aᵢ Z) irreps → external tensor products
#                → DIRECT SUM enumeration with multiplicity (⊕_α m_α ρ_α)
#   Admissibility: Verlinde, (ST)³ = ζ_24^c · S² with c ∈ 0..23, Cauchy, FS
# Level II (numerical, ComplexF64):
#   Pentagon — two solver backends:
#     :homotopy (default) — HomotopyContinuation.jl, polyhedral start system
#                           gives (near-)exhaustive coverage; optional certify
#     :newton             — damped Newton with KrylovKit (legacy, kept for comparison)
# Level III (numerical, ComplexF64):
#   Hexagon — HomotopyContinuation.jl, R-symbols only (F fixed from Level II)
#
# Usage:
#   include("mtc_classifier.jl")
#   classify_mtc(5)                           # homotopy by default
#   classify_mtc(9; solver=:newton)           # legacy damped Newton
#   classify_mtc(5; certify_solutions=true)   # rigorous certification

using LinearAlgebra
using SparseArrays
using Oscar
using TensorCategories
using KrylovKit
# Use `import` (not `using`) so that names like `degree`, `nvariables`, etc.
# from HomotopyContinuation/DynamicPolynomials do NOT clash with Oscar's.
import HomotopyContinuation
const HC = HomotopyContinuation

# ============================================================
#  Level I: Exact modular data over Q(ζ_N)
# ============================================================

function load_sl2reps()
    GAP.Globals.LoadPackage(GapObj("SL2Reps"))
end

# --- GAP → Oscar ---

function gap_to_oscar_matrix(gap_mat, r, K::AbsSimpleNumField, cond_target::Int)
    """Convert GAP matrix entries into given cyclotomic field K = Q(ζ_cond_target)."""
    z = gen(K)
    M = zero_matrix(K, r, r)
    
    for i in 1:r, j in 1:r
        entry_cond = Int(GAP.Globals.Conductor(gap_mat[i,j]))
        # entry_cond must divide cond_target
        if cond_target % entry_cond != 0
            error("Entry conductor $entry_cond does not divide target $cond_target")
        end
        mult = div(cond_target, entry_cond)
        coeffs = GAP.Globals.CoeffsCyc(gap_mat[i,j], entry_cond)
        val = zero(K)
        for k in 1:entry_cond
            if !Bool(GAP.Globals.IsZero(coeffs[k]))
                num = Int(GAP.Globals.NumeratorRat(coeffs[k]))
                den = Int(GAP.Globals.DenominatorRat(coeffs[k]))
                # ζ_{entry_cond}^(k-1) = ζ_{cond_target}^((k-1)*mult)
                val += QQ(num, den) * z^((k-1) * mult)
            end
        end
        M[i,j] = val
    end
    return M
end

function matrix_conductor(gap_mat, r)
    cond = 1
    for i in 1:r, j in 1:r
        cond = lcm(cond, Int(GAP.Globals.Conductor(gap_mat[i,j])))
    end
    cond
end

# --- Embed K = Q(ζ_N) into ComplexF64 ---

function to_complex(val::AbsSimpleNumFieldElem, zeta::ComplexF64, deg::Int)
    result = 0.0 + 0.0im
    for k in 0:deg-1
        c = coeff(val, k)
        result += (Float64(numerator(c)) / Float64(denominator(c))) * zeta^k
    end
    result
end

function matrix_to_complex(M, zeta::ComplexF64, deg::Int)
    r = nrows(M)
    c = ncols(M)
    result = zeros(ComplexF64, r, c)
    for i in 1:r, j in 1:c
        result[i,j] = to_complex(M[i,j], zeta, deg)
    end
    result
end

# --- CRT decomposition ---

function prime_power_factors(N::Int)
    result = Tuple{Int,Int}[]
    for (p, a) in factor(N)
        push!(result, (Int(p), Int(a)))
    end
    result
end

function irreps_for_prime_power(p::Int, a::Int)
    """Get all irreps of SL₂(Z/p^a Z) (levels dividing p^a)."""
    irreps = []
    for lam in 1:a
        lev = p^lam
        reps = GAP.evalstr("SL2IrrepsOfLevel($lev)")
        len = Int(GAP.Globals.Length(reps))
        for i in 1:len
            push!(irreps, reps[i])
        end
    end
    reps_1 = GAP.evalstr("SL2IrrepsOfLevel(1)")
    len_1 = Int(GAP.Globals.Length(reps_1))
    for i in 1:len_1
        push!(irreps, reps_1[i])
    end
    irreps
end

# --- Kronecker product in K ---

function kron_oscar(A, B, K)
    ra = nrows(A); ca = ncols(A)
    rb = nrows(B); cb = ncols(B)
    C = zero_matrix(K, ra*rb, ca*cb)
    for i in 1:ra, j in 1:ca, k in 1:rb, l in 1:cb
        C[(i-1)*rb + k, (j-1)*cb + l] = A[i,j] * B[k,l]
    end
    C
end

# --- Verlinde over K ---

function verlinde_exact(S, r::Int, K)
    """Compute fusion coefficients N_{ij}^k = Σ_l S[i,l] S[j,l] conj(S[k,l]) / S[1,l].
    Returns Nijk as Array{Int,3} if all entries are non-negative integers, else nothing.
    conj in Q(ζ_N) = complex conjugation = ζ_N → ζ_N^{-1}."""
    
    # Check S[1,l] nonzero for all l
    for l in 1:r
        iszero(S[1,l]) && return nothing
    end
    
    # Complex conjugation on K = Q(ζ_N): σ(ζ) = ζ^{-1} = ζ^{N-1}
    N = Int(get_attribute(K, :cyclo, -1))
    if N < 0
        # Fallback: infer N from degree via Euler's totient
        # but cyclotomic_field stores the order
        N = degree(K)  # Not quite right; this is φ(N). Need to set attribute.
    end
    
    function conj_K(x)
        # If x = Σ c_k ζ^k, then σ(x) = Σ c_k ζ^{-k}
        deg = degree(K)
        z = gen(K)
        # ζ^{-1} = ζ^{N-1} where N is the conductor
        # We'll compute it directly from the minimal polynomial
        # Easier: use the Hecke function if available
        return Hecke.conjugate(x)
    end
    
    # Try to use complex_conjugation from Hecke
    Nijk = zeros(Int, r, r, r)
    for i in 1:r, j in 1:r, k in 1:r
        val = zero(K)
        for l in 1:r
            val += S[i,l] * S[j,l] * conj_K(S[k,l]) // S[1,l]
        end
        # val must be an integer
        if !is_integer(val)
            return nothing
        end
        n = ZZ(val)
        if n < 0 || n > 1000
            return nothing
        end
        Nijk[i,j,k] = Int(n)
    end
    Nijk
end

function is_integer(x::AbsSimpleNumFieldElem)
    deg = degree(parent(x))
    # x = c_0 + c_1 ζ + ... must have c_0 ∈ Q and c_1 = c_2 = ... = 0
    for k in 1:deg-1
        !iszero(coeff(x, k)) && return false
    end
    # c_0 must be integer
    c0 = coeff(x, 0)
    return isone(denominator(c0))
end

function to_ZZ(x::AbsSimpleNumFieldElem)
    c0 = coeff(x, 0)
    ZZ(numerator(c0))
end

# Simplified Verlinde: compute complex, check closeness to integer
function verlinde_float_from_exact(S_exact, K, N::Int, r::Int; tol=1e-6)
    zeta = exp(2π * im / N)
    deg = degree(K)
    S_num = matrix_to_complex(S_exact, zeta, deg)
    
    for l in 1:r
        abs(S_num[1,l]) < 1e-12 && return nothing
    end
    
    Nijk = zeros(Int, r, r, r)
    for i in 1:r, j in 1:r, k in 1:r
        val = sum(S_num[i,l] * S_num[j,l] * conj(S_num[k,l]) / S_num[1,l] for l in 1:r)
        n = round(Int, real(val))
        if abs(real(val) - n) > tol || abs(imag(val)) > tol
            return nothing
        end
        if n < 0
            return nothing
        end
        Nijk[i,j,k] = n
    end
    Nijk
end

# --- Sign fix ---

function fix_signs_exact(S_exact, K, N::Int, r::Int)
    """Determine sign of each quantum dimension numerically, then fix S exactly."""
    zeta = exp(2π * im / N)
    deg = degree(K)
    
    signs = ones(Int, r)
    S11 = to_complex(S_exact[1,1], zeta, deg)
    for i in 1:r
        S1i = to_complex(S_exact[1,i], zeta, deg)
        d = S1i / S11
        if real(d) < -1e-10
            signs[i] = -1
        end
    end
    
    S_fixed = zero_matrix(K, r, r)
    for i in 1:r, j in 1:r
        S_fixed[i,j] = signs[i] * signs[j] * S_exact[i,j]
    end
    return S_fixed, signs
end

# --- Admissibility: (ST)³ = λ·S² with λ an 8-th root of unity ---

"""
    check_ST_exact(S, T, K, N, r) -> (ok::Bool, c::Union{Int,Nothing})

Check whether there exists a central charge c ∈ {0, 1, ..., 7} such that
`(ST)³ = ζ_8^c · S²` over K = Q(ζ_N). Returns `(true, c)` on the first
match, `(false, nothing)` otherwise.

Physically, c mod 8 is the central charge of the MTC. The old version only
tested c ∈ {0, 4} (i.e. λ = ±1) and thereby rejected all MTCs with c ≢ 0
mod 4 (e.g. Ising has c = 1/2, so 16 · c = 8, which corresponds to
λ = ζ_16 — actually the 16-th root of unity case; see note below).

Note: some conventions use 24-th roots (for SL₂(Z) projective class in
central extensions), others use 8-th roots (for the MTC multiplicative
central charge e^{2πic/8}). We scan both 8 and 24 to be safe. If N is
divisible by 8 (resp 24), ζ_8 (resp ζ_24) lives in K and can be constructed;
if not, we still try to solve for λ as a root of `λ · S² - (ST)³ = 0`
entry-wise, which only works when (ST)³ and S² are proportional.
"""
function check_ST_exact(S, T, K, N::Int, r::Int)
    ST = S * T
    ST3 = ST * ST * ST
    S2 = S * S

    # (ST)³ = λ · S² must hold entry-wise with one common scalar λ.
    # Find the first nonzero entry of S² and compute λ candidate from it.
    λ_candidate = nothing
    for i in 1:r, j in 1:r
        if !iszero(S2[i, j])
            # λ = ST3[i,j] / S2[i,j]
            λ_candidate = ST3[i, j] // S2[i, j]
            # Verify λ is consistent across ALL entries
            consistent = true
            for a in 1:r, b in 1:r
                if !iszero(S2[a, b])
                    if !iszero(ST3[a, b] - λ_candidate * S2[a, b])
                        consistent = false
                        break
                    end
                else
                    if !iszero(ST3[a, b])
                        consistent = false
                        break
                    end
                end
            end
            if !consistent
                return (false, nothing)
            end
            break
        end
    end

    if λ_candidate === nothing
        # S² is zero → pathological; reject
        return (false, nothing)
    end

    # Check if λ_candidate is a root of unity of order dividing 24.
    # We try ζ_8^c for c = 0..7 first (most common), then ζ_24^c for full scan.
    # Since ζ_n ∈ Q(ζ_N) iff n | N (when we embed canonically), we can only
    # directly test roots of unity whose order divides N.
    # Pragmatic approach: check λ^8 = 1, then return c such that λ = ζ_8^c.

    # Compute λ^8 - 1 and test exactness.
    λ8 = λ_candidate^8
    if !isone(λ8)
        # Try λ^24 = 1 as a fallback
        λ24 = λ_candidate^24
        if !isone(λ24)
            return (false, nothing)
        end
    end

    # λ is a root of unity; find c ∈ 0..23 such that λ = ζ_24^c in K (if 24 | N),
    # else fall back to matching c ∈ 0..7 via numeric embedding.
    # Numeric match (always works):
    zeta_N_num = exp(2π * im / N)
    deg = degree(K)
    λ_num = to_complex(λ_candidate, zeta_N_num, deg)
    for c in 0:23
        if abs(λ_num - exp(2π * im * c / 24)) < 1e-8
            return (true, c)
        end
    end
    return (false, nothing)
end

# --- Complex conjugation on Q(ζ_N) ---

function conj_cyclotomic(x, K, N::Int)
    """Complex conjugation: ζ_N → ζ_N^{-1}."""
    deg = degree(K)
    z = gen(K)
    # ζ_N^{-1} = ζ_N^{N-1}
    zinv = z^(N-1)
    result = zero(K)
    for k in 0:deg-1
        c = coeff(x, k)
        result += c * zinv^k
    end
    result
end

function conj_matrix(M, K, N::Int)
    r = nrows(M); c = ncols(M)
    result = zero_matrix(K, r, c)
    for i in 1:r, j in 1:c
        result[i,j] = conj_cyclotomic(M[i,j], K, N)
    end
    result
end

# --- Cauchy theorem: Galois norm check ---

function check_cauchy_exact(S, K, N::Int, r::Int)
    """Cauchy theorem via Galois norm of D²."""
    
    S11 = S[1,1]
    S11_conj = conj_cyclotomic(S11, K, N)
    # |S[1,1]|² ∈ Q(ζ_N) ∩ R = Q(ζ_N + ζ_N^{-1})、有理数とは限らない
    D_sq_inv = S11 * S11_conj
    
    # D² = 1 / D_sq_inv ∈ K
    D_sq = inv(D_sq_inv)
    
    # Galois norm over Q: N(D²) = ∏_{σ ∈ Gal(K/Q)} σ(D²)
    # Gal(Q(ζ_N)/Q) = (Z/NZ)× で、σ_a: ζ_N → ζ_N^a
    
    norm_val = one(K)
    for a in 1:N-1
        gcd(a, N) != 1 && continue
        norm_val *= galois_action(D_sq, K, N, a)
    end
    
    # norm_val should be rational
    if !is_rational(norm_val, K)
        return false  # shouldn't happen
    end
    norm_Q = extract_rational(norm_val)
    
    # Prime factors of numerator and denominator ⊆ primes of N
    N_primes = Set(Int(p) for (p, _) in factor(ZZ(N)))
    
    num = numerator(norm_Q)
    den = denominator(norm_Q)
    
    if !iszero(num) && abs(num) != 1
        for (p, _) in factor(abs(num))
            !(Int(p) in N_primes) && return false
        end
    end
    if den != 1
        for (p, _) in factor(abs(den))
            !(Int(p) in N_primes) && return false
        end
    end
    
    return true
end

function galois_action(x, K, N::Int, a::Int)
    """σ_a: ζ_N → ζ_N^a applied to x ∈ Q(ζ_N)."""
    deg = degree(K)
    z = gen(K)
    za = z^a  # ζ^a
    
    # x = Σ c_k ζ^k → Σ c_k (ζ^a)^k = Σ c_k ζ^{ak}
    result = zero(K)
    for k in 0:deg-1
        c = coeff(x, k)
        result += c * za^k
    end
    result
end

function is_rational(x, K)
    deg = degree(K)
    for k in 1:deg-1
        !iszero(coeff(x, k)) && return false
    end
    return true
end

function extract_rational(x)
    coeff(x, 0)
end

# --- Frobenius-Schur indicator (exact) ---

function check_frobenius_schur_exact(S, T, Nijk, K, N::Int, r::Int; n=2)
    """ν_n(V_k) = (1/D²) Σ_{i,j} N_{ij}^k d_i d_j (θ_i/θ_j)^n must be integer in Z[ζ_N]."""
    S11 = S[1,1]
    d = [S[1,Y] // S11 for Y in 1:r]
    D_sq = sum(d[Y] * d[Y] for Y in 1:r)
    
    for k_idx in 1:r
        nu_n = zero(K)
        for i in 1:r, j in 1:r
            theta_ratio_n = T[i,i]^n // T[j,j]^n
            nu_n += Nijk[i,j,k_idx] * d[i] * d[j] * theta_ratio_n
        end
        nu_n = nu_n // D_sq
        
        # ν_n(k) must be in Z[ζ_N], i.e., an algebraic integer of K
        # For ν_2: ν_2(k) = 0 if k ≠ k*, ν_2(k) = ±1 if k = k*
        
        # Check algebraic integrality: easier is to check rationality for n=2 self-dual case
        # Full check: ν_n should be an algebraic integer in Z[ζ_N]
        # For simplicity, check that it's integer for self-dual labels
        
        if !is_algebraic_integer_in_Q_zeta_N(nu_n, K)
            return false
        end
    end
    return true
end

function is_algebraic_integer_in_Q_zeta_N(x, K)
    deg = degree(K)
    for k in 0:deg-1
        c = coeff(x, k)
        !isone(denominator(c)) && return false
    end
    return true
end

# --- Level I main ---

"""
    _block_diag_matrices(M_list, K) -> big_matrix

Given a list of square matrices over K, produce the block-diagonal matrix.
Used to assemble ⊕_α m_α · (S_α) into a single (S) matrix.
"""
function _block_diag_matrices(M_list, K)
    total_size = sum(nrows(M) for M in M_list)
    big = zero_matrix(K, total_size, total_size)
    offset = 0
    for M in M_list
        sz = nrows(M)
        for i in 1:sz, j in 1:sz
            big[offset + i, offset + j] = M[i, j]
        end
        offset += sz
    end
    return big
end

"""
    _enumerate_multiplicity_vectors(rank_list, max_total_rank) -> iterator

Enumerate non-negative integer vectors (m_1, ..., m_K) such that
∑ m_k · rank_list[k] ≤ max_total_rank, excluding the zero vector.
"""
function _enumerate_multiplicity_vectors(rank_list::Vector{Int}, max_total_rank::Int)
    result = Vector{Vector{Int}}()
    K = length(rank_list)
    function recurse(idx::Int, remaining::Int, current::Vector{Int})
        if idx > K
            if sum(current) > 0  # at least one nonzero multiplicity
                push!(result, copy(current))
            end
            return
        end
        max_m = div(remaining, rank_list[idx])
        for m in 0:max_m
            current[idx] = m
            recurse(idx + 1, remaining - m * rank_list[idx], current)
        end
        current[idx] = 0
    end
    recurse(1, max_total_rank, zeros(Int, K))
    return result
end

function enumerate_modular_data(N::Int; max_rank::Int=20,
                                 verbose::Bool=false,
                                 enumerate_direct_sums::Bool=true,
                                 max_multiplicity::Int=max_rank)
    pf = prime_power_factors(N)
    println("N = $N = ", join(["$(p)^$(a)" for (p,a) in pf], " × "))
    
    # The target cyclotomic field is Q(ζ_N)
    K, z = cyclotomic_field(N)
    
    # Step 1: Get (S, T) for each prime power factor as matrices over K
    factor_irreps = []
    for (p, a) in pf
        gap_irreps = irreps_for_prime_power(p, a)
        st_pairs = []
        for rep in gap_irreps
            r = Int(rep.degree)
            # Determine common conductor for this rep (S and T entries)
            cond_S = matrix_conductor(rep.S, r)
            cond_T = matrix_conductor(rep.T, r)
            rep_cond = lcm(cond_S, cond_T)
            # rep_cond must divide N
            if N % rep_cond != 0
                continue
            end
            S = gap_to_oscar_matrix(rep.S, r, K, N)
            T = gap_to_oscar_matrix(rep.T, r, K, N)
            push!(st_pairs, (S, T, r))
        end
        push!(factor_irreps, st_pairs)
        println("  SL₂(Z/$(p^a)Z): $(length(st_pairs)) irreps")
    end
    
    # Step 2: Build catalog of "atomic" irreps of SL₂(Z/NZ) via external tensor product.
    # Each entry is (S_α, T_α, r_α, name_α).
    atomic_catalog = []
    if length(factor_irreps) == 1
        for (S, T, r) in factor_irreps[1]
            push!(atomic_catalog, (S=S, T=T, rank=r, name="$(r)d"))
        end
    else
        indices = [1:length(fi) for fi in factor_irreps]
        for idx_tuple in Iterators.product(indices...)
            total_rank = prod(factor_irreps[k][idx_tuple[k]][3] for k in 1:length(pf))
            total_rank > max_rank && continue
            S = factor_irreps[1][idx_tuple[1]][1]
            T = factor_irreps[1][idx_tuple[1]][2]
            for k in 2:length(factor_irreps)
                S = kron_oscar(S, factor_irreps[k][idx_tuple[k]][1], K)
                T = kron_oscar(T, factor_irreps[k][idx_tuple[k]][2], K)
            end
            name = join(["$(factor_irreps[k][idx_tuple[k]][3])d" for k in 1:length(pf)], "⊗")
            push!(atomic_catalog, (S=S, T=T, rank=total_rank, name=name))
        end
    end
    println("  Atomic irreps of SL₂(Z/$(N)Z) (≤ rank $max_rank): $(length(atomic_catalog))")

    # Step 3: Enumerate direct sums with multiplicity.
    valid_candidates = []
    reject_count = Dict("Verlinde"=>0, "(ST)³"=>0, "Cauchy"=>0, "FS"=>0)

    function process_candidate(S, T, r, name="")
        r > max_rank && return
        
        # Fix signs using float heuristic
        S_fixed, signs = fix_signs_exact(S, K, N, r)
        
        # Verlinde (use float for speed, then verify integers)
        Nijk = verlinde_float_from_exact(S_fixed, K, N, r)
        if Nijk === nothing
            reject_count["Verlinde"] += 1
            return
        end
        if !all(Nijk .>= 0)
            reject_count["Verlinde"] += 1
            return
        end
        
        # (ST)³ = ζ_24^c · S² for some c ∈ 0..23
        ok_st, central_charge = check_ST_exact(S_fixed, T, K, N, r)
        if !ok_st
            reject_count["(ST)³"] += 1
            verbose && println("  rank=$r rejected: (ST)³")
            return
        end
        
        # Cauchy
        if !check_cauchy_exact(S_fixed, K, N, r)
            reject_count["Cauchy"] += 1
            verbose && println("  rank=$r rejected: Cauchy")
            return
        end
        
        # Frobenius-Schur
        if !check_frobenius_schur_exact(S_fixed, T, Nijk, K, N, r)
            reject_count["FS"] += 1
            verbose && println("  rank=$r rejected: FS")
            return
        end
        
        # Float versions for Level II
        zeta = exp(2π * im / N)
        deg = degree(K)
        S_float = matrix_to_complex(S_fixed, zeta, deg)
        T_float = matrix_to_complex(T, zeta, deg)
        
        push!(valid_candidates, (
            S=S_float, T=T_float, S_exact=S_fixed, T_exact=T,
            Nijk=Nijk, rank=r, name=name,
            central_charge=central_charge,
        ))
    end

    rank_list = [cand.rank for cand in atomic_catalog]

    if !enumerate_direct_sums
        # Legacy behavior: try each atomic irrep alone.
        for (α, cand) in enumerate(atomic_catalog)
            process_candidate(cand.S, cand.T, cand.rank, cand.name)
        end
    else
        mult_vectors = _enumerate_multiplicity_vectors(rank_list, max_rank)
        # Filter out too-large per-irrep multiplicities
        mult_vectors = filter(v -> all(m -> m <= max_multiplicity, v), mult_vectors)
        println("  Multiplicity vectors to try: $(length(mult_vectors))")
        for mv in mult_vectors
            # Build block diagonal (S, T) from the multiplicities.
            blocks_S = AbstractAlgebra.Generic.MatSpaceElem{AbsSimpleNumFieldElem}[]
            blocks_T = AbstractAlgebra.Generic.MatSpaceElem{AbsSimpleNumFieldElem}[]
            name_parts = String[]
            total_rank = 0
            for (α, m) in enumerate(mv)
                m == 0 && continue
                for _ in 1:m
                    push!(blocks_S, atomic_catalog[α].S)
                    push!(blocks_T, atomic_catalog[α].T)
                    total_rank += atomic_catalog[α].rank
                end
                push!(name_parts, m == 1 ? atomic_catalog[α].name : "$(m)·$(atomic_catalog[α].name)")
            end
            total_rank > max_rank && continue
            S_big = _block_diag_matrices(blocks_S, K)
            T_big = _block_diag_matrices(blocks_T, K)
            name = "⊕ " * join(name_parts, " ⊕ ")
            process_candidate(S_big, T_big, total_rank, name)
        end
    end
    
    println("\nRejections: Verlinde=$(reject_count["Verlinde"]), " *
            "(ST)³=$(reject_count["(ST)³"]), " *
            "Cauchy=$(reject_count["Cauchy"]), " *
            "FS=$(reject_count["FS"])")
    println("$(length(valid_candidates)) valid modular data candidates.")
    for (i, cand) in enumerate(valid_candidates)
        name = hasproperty(cand, :name) ? cand.name : ""
        cc = hasproperty(cand, :central_charge) ? " (c=$(cand.central_charge)/8)" : ""
        println("  [$i] rank=$(cand.rank) $name$cc")
    end
    
    return valid_candidates
end

# ============================================================
#  Level II: Pentagon solver (damped Newton over C, sparse Krylov)
# ============================================================

function get_pentagon_system(Nijk::Array{Int,3}, r::Int)
    one_vec = zeros(Int, r)
    one_vec[1] = 1
    C, eqs_raw = pentagon_equations(Nijk, one_vec)
    eqs = filter(eq -> !(eq isa Integer) && !iszero(eq), eqs_raw)
    R = parent(eqs[1])
    n = nvars(R)
    return R, eqs, n
end

function eval_poly_complex(f, vals::Vector{ComplexF64})
    isa(f, Integer) && return ComplexF64(f)
    iszero(f) && return ComplexF64(0.0)
    result = 0.0 + 0.0im
    for (c, m) in zip(coefficients(f), monomials(f))
        degs = degrees(m)
        term = ComplexF64(Float64(numerator(c)) / Float64(denominator(c)))
        for i in 1:length(degs)
            degs[i] > 0 && (term *= vals[i]^degs[i])
        end
        result += term
    end
    result
end

function sparse_jacobian(eqs, derivs, x, n)
    m = length(eqs)
    I = Int[]
    J_idx = Int[]
    V = ComplexF64[]
    for i in 1:m, j in 1:n
        v = eval_poly_complex(derivs[i][j], x)
        if abs(v) > 1e-15
            push!(I, i)
            push!(J_idx, j)
            push!(V, v)
        end
    end
    sparse(I, J_idx, V, m, n)
end

function solve_pentagon_newton(eqs, n, zeta; max_trials=20, max_iter=200)
    derivs = [[derivative(eq, j) for j in 1:n] for eq in eqs]
    
    solutions = []
    for trial in 1:max_trials
        x = ComplexF64[rand(-3:3) + rand(-3:3)*zeta for _ in 1:n]
        for i in 1:n
            abs(x[i]) < 0.1 && (x[i] = 1.0 + zeta)
        end
        
        for iter in 1:max_iter
            F_val = ComplexF64[eval_poly_complex(eq, x) for eq in eqs]
            res = maximum(abs.(F_val))
            if res < 1e-12
                println("  Trial $trial: CONVERGED at iter $iter")
                push!(solutions, copy(x))
                break
            end
            
            J = sparse_jacobian(eqs, derivs, x, n)
            delta, info = linsolve(v -> J' * (J * v), J' * F_val;
                                   ishermitian=true, isposdef=true, verbosity=0)
            
            alpha = 1.0
            for _ in 1:20
                x_new = x - alpha * delta
                F_new = ComplexF64[eval_poly_complex(eq, x_new) for eq in eqs]
                if maximum(abs.(F_new)) < res
                    x = x_new
                    break
                end
                alpha *= 0.5
            end
        end
    end
    solutions
end

# ============================================================
#  Level II (alt): Pentagon solver via HomotopyContinuation.jl
# ============================================================

"""
    oscar_poly_to_hc(f, hc_vars) -> HC.Expression

Convert an Oscar multivariate polynomial `f` (with QQ coefficients) into a
HomotopyContinuation.jl Expression in the given `hc_vars` (Vector{HC.Variable}).
Integer constants are accepted as well (returned as constant expressions).
"""
function oscar_poly_to_hc(f, hc_vars::Vector{HC.Variable})
    isa(f, Integer) && return HC.Expression(f)
    iszero(f) && return HC.Expression(0)
    expr = HC.Expression(0)
    for (c, m) in zip(coefficients(f), monomials(f))
        # coefficient: QQFieldElem -> Rational{BigInt} -> Float via ComplexF64
        num = BigInt(numerator(c))
        den = BigInt(denominator(c))
        coef_expr = HC.Expression(num // den)
        degs = degrees(m)
        term = coef_expr
        for i in 1:length(degs)
            if degs[i] > 0
                term *= hc_vars[i]^degs[i]
            end
        end
        expr += term
    end
    expr
end

"""
    build_hc_system(eqs, n) -> HC.System

Build a HomotopyContinuation.jl `System` from Oscar polynomial equations `eqs`
in `n` variables. Filters out trivial (zero/integer) equations.
"""
function build_hc_system(eqs, n::Int)
    # Build variables directly without the @var macro (safer under `import`).
    hc_vars = [HC.Variable(Symbol("x", i)) for i in 1:n]
    hc_exprs = HC.Expression[]
    for eq in eqs
        (isa(eq, Integer) || iszero(eq)) && continue
        push!(hc_exprs, oscar_poly_to_hc(eq, hc_vars))
    end
    return HC.System(hc_exprs; variables=hc_vars), hc_vars
end

"""
    solve_pentagon_homotopy(eqs, n; kwargs...) -> Vector{Vector{ComplexF64}}

Solve the pentagon equations via polynomial homotopy continuation.
Returns finite solutions as Vector{ComplexF64}.

Pentagon equations as returned by TensorCategories.jl typically carry a
**gauge symmetry** (diagonal F-basis change), so the raw variety is
positive-dimensional and HC reports all endpoints as "singular". To reduce
to a zero-dimensional system we support two strategies:

- `slice = 0` (default): no gauge fixing. Useful if you know the system is
  already gauge-fixed, or you want the raw behavior to inspect.
- `slice = k > 0`: append `k` random complex-linear equations `aᵀx = b` to
  cut k gauge dimensions. A standard algebraic-geometry trick for slicing
  positive-dimensional components. Each run uses a fresh random slice so
  re-running probes different representatives.
- `include_singular = true`: also return solutions flagged as singular
  (last-resort when gauge slicing is not enough).

Keyword arguments:
- `slice = 0`: number of random linear slices to append
- `include_singular = false`: if true, return singular endpoints too
- `certify_solutions = false`: run `HC.certify` for rigorous validation (only meaningful for nonsingular)
- `threading = true`: parallel path tracking
- `start_system = :polyhedral`: :polyhedral (sparse) or :total_degree
- `show_progress = false`
"""
function solve_pentagon_homotopy(eqs, n::Int;
                                 slice::Int=0,
                                 include_singular::Bool=false,
                                 certify_solutions::Bool=false,
                                 threading::Bool=true,
                                 start_system::Symbol=:polyhedral,
                                 show_progress::Bool=false)
    sys, hc_vars = build_hc_system(eqs, n)

    # Optionally append random linear slices to cut gauge directions.
    if slice > 0
        extra = HC.Expression[]
        for _ in 1:slice
            # Random complex-Gaussian linear form aᵀx - b
            a = randn(ComplexF64, n)
            b = randn(ComplexF64)
            lf = sum(HC.Expression(a[i]) * hc_vars[i] for i in 1:n) - HC.Expression(b)
            push!(extra, lf)
        end
        all_eqs = vcat(HC.expressions(sys), extra)
        sys = HC.System(all_eqs; variables=hc_vars)
    end

    result = HC.solve(sys;
                      start_system=start_system,
                      threading=threading,
                      show_progress=show_progress)

    # Collect the requested slice of solutions
    sols_raw = if include_singular
        HC.solutions(result)
    else
        HC.solutions(result; only_nonsingular=true)
    end

    # Optional certification (nonsingular only)
    if certify_solutions && !isempty(sols_raw)
        cert = HC.certify(sys, sols_raw)
        certs = HC.certificates(cert)
        sols_raw = [HC.solution_approximation(c) for c in certs if HC.is_certified(c)]
        println("  Certified: $(length(sols_raw)) / $(length(certs)) solutions")
    end

    # Brief summary
    nsols_total = HC.nsolutions(result)
    nsing = HC.nsingular(result)
    natinf = HC.nat_infinity(result)
    nfail = HC.nfailed(result)
    slicestr = slice > 0 ? " [+$slice slice]" : ""
    println("  HC$slicestr: $(length(sols_raw)) returned " *
            "(total=$nsols_total, singular=$nsing, at_infinity=$natinf, failed=$nfail)")

    return [ComplexF64.(s) for s in sols_raw]
end

"""
    refine_solution_newton(eqs, x0; tol=1e-14, max_iter=50) -> Vector{ComplexF64}

Polish a HC solution with a few damped-Newton steps at ComplexF64 precision.
Useful when HC's double-precision output needs tightening before PSLQ /
algebraic recognition downstream.
"""
function refine_solution_newton(eqs, x0::Vector{ComplexF64};
                                tol::Float64=1e-14, max_iter::Int=50)
    n = length(x0)
    derivs = [[derivative(eq, j) for j in 1:n] for eq in eqs]
    x = copy(x0)
    for _ in 1:max_iter
        F_val = ComplexF64[eval_poly_complex(eq, x) for eq in eqs]
        res = maximum(abs.(F_val))
        res < tol && return x
        J = sparse_jacobian(eqs, derivs, x, n)
        delta, _ = linsolve(v -> J' * (J * v), J' * F_val;
                            ishermitian=true, isposdef=true, verbosity=0)
        alpha = 1.0
        for _ in 1:20
            x_new = x - alpha * delta
            F_new = ComplexF64[eval_poly_complex(eq, x_new) for eq in eqs]
            if maximum(abs.(F_new)) < res
                x = x_new
                break
            end
            alpha *= 0.5
        end
    end
    return x
end

# ============================================================
#  Level III: Hexagon solver (F fixed from Level II, solve for R)
# ============================================================

"""
    _coerce_complex(c, K) -> K element

Convert a Julia `Number` (typically ComplexF64) into an element of the ring K.
If K is an AcbField, goes through real/imag parts. If K is a polynomial ring
over AcbField, first goes to the scalar AcbField then lifts to K.
"""
function _coerce_complex(c::Number, K)
    scalar = K isa AcbField ? K : base_ring(K)
    # Build scalar(real) + i * scalar(imag) in the AcbField
    re = scalar(Float64(real(c)))
    im_part = scalar(Float64(imag(c)))
    acb_val = re + onei(scalar) * im_part
    # Lift into K (identity if K is already the scalar field)
    return K(acb_val)
end

"""
    invert_associator_numeric(F_values, mult, one) -> Vector{ComplexF64}

Compute the inverse F-symbol values. Given `F_values` that would produce the
associator via `assign_F_to_associator!`, return a vector `F_inv_values` such
that `assign_F_to_associator!(_, F_inv_values)` produces the inverse
associator (block-wise inverse over ComplexF64).

Implementation: simulate `assign_F_to_associator!` exactly to extract each
block in-order, invert, and rebuild a fresh vector that would regenerate
those inverted blocks under the same `assign` convention.
"""
function invert_associator_numeric(F_values::Vector{ComplexF64},
                                   mult::Array{Int,3}, one_vec::Vector{Int})
    dummy = TensorCategories.six_j_category(QQ, mult)
    dummy.one = one_vec
    m = dummy.simples

    # First pass: extract blocks in the same order `assign_F_to_associator!` uses.
    # Each block is built from `entries = [pop!(y) for _ in 1:(r*t)]` and then
    # `matrix(K, r, t, entries)` fills row-major.
    y_stack = copy(F_values)
    extracted_blocks = Vector{Tuple{Int,Int,Matrix{ComplexF64}}}()  # (r, t, M)
    for i in 1:m, j in 1:m, k in 1:m, o in 1:m
        sum(one_vec[[i, j, k]]) > 0 && continue
        (r, t) = size(dummy.ass[i, j, k, o])
        @assert r == t "F-block at ($i,$j,$k,$o) is $r×$t — not square"
        # Mirror assign exactly:
        entries = ComplexF64[pop!(y_stack) for _ in 1:(r * t)]
        # `matrix(K, r, t, entries)` in Oscar fills row-major, so:
        # M[a, b] = entries[(a-1)*t + b]
        M = Matrix{ComplexF64}(undef, r, t)
        for a in 1:r, b in 1:t
            M[a, b] = entries[(a - 1) * t + b]
        end
        push!(extracted_blocks, (r, t, M))
    end
    @assert isempty(y_stack) "F_values had leftover entries"

    # Invert each block
    inverted_blocks = [(r, t, inv(M)) for (r, t, M) in extracted_blocks]

    # Rebuild a vector that, when passed to assign_F_to_associator!, reproduces
    # the inverted blocks. We need to "un-pop" the pop! sequence.
    # For each block (r, t, Minv):
    #   entries_new = [Minv[a,b] in row-major order] = flat_rowmajor(Minv)
    #   These entries must appear in F_inv_values such that `pop!` consumes them
    #   in REVERSE, i.e., F_inv_values[end] -> entries_new[1], etc.
    # So for this block alone, we append `reverse(entries_new)` to the OUTPUT
    # END... but blocks are processed in-order and each block is placed at
    # increasingly deeper positions. Concretely:
    # At start, y_stack = F_values. First block pops r1*t1 from the END.
    # → Those entries are F_values[end-r1*t1+1 : end] (in original order).
    # In pop order: F_values[end], F_values[end-1], ..., F_values[end-r1*t1+1]
    # become entries_new[1], entries_new[2], ..., entries_new[r1*t1].
    # So F_values[end - (p-1)] = entries_new[p]  for p = 1..r1*t1
    # Therefore to ENCODE entries_new into a slice of F_inv_values at the end:
    # F_inv_values[end - (p-1)] = inv_entries[p]  (same relation for the inverted block)

    F_inv = Vector{ComplexF64}(undef, length(F_values))
    cursor_end = length(F_values)
    for (r, t, Minv) in inverted_blocks
        inv_entries = Vector{ComplexF64}(undef, r * t)
        for a in 1:r, b in 1:t
            inv_entries[(a - 1) * t + b] = Minv[a, b]
        end
        # Place into F_inv[cursor_end - r*t + 1 : cursor_end] such that
        # F_inv[cursor_end - (p-1)] = inv_entries[p]
        for p in 1:(r * t)
            F_inv[cursor_end - (p - 1)] = inv_entries[p]
        end
        cursor_end -= r * t
    end
    @assert cursor_end == 0
    return F_inv
end

"""
    assign_F_to_associator!(poly_C, F_values)

Write a ComplexF64 F-symbol solution into `poly_C.ass` following the **same
traversal order** used by `TensorCategories.pentagon_equations`, so that the
i-th component of `F_values` lands in the slot that was called `x_i` during
Level II. The `pop!` convention matches `pentagon_equations` (pops from the
end of a copy of the variable vector).
"""
function assign_F_to_associator!(poly_C, F_values::Vector{<:Number})
    K = base_ring(poly_C)
    m = poly_C.simples
    one_vec = poly_C.one
    y = copy(F_values)     # treat as a stack; pop! takes from the end
    for i in 1:m, j in 1:m, k in 1:m, o in 1:m
        sum(one_vec[[i, j, k]]) > 0 && continue
        (r, t) = size(poly_C.ass[i, j, k, o])
        entries = [_coerce_complex(pop!(y), K) for _ in 1:(r * t)]
        poly_C.ass[i, j, k, o] = matrix(K, r, t, entries)
    end
    @assert isempty(y) "F_values length ($(length(F_values))) does not match the expected number of F-symbol slots"
    return poly_C
end

"""
    _number_of_variables_in_hexagon_equations(poly_C) -> Int

Count independent R-symbol entries. For each (i, j), the braiding lives in
`C.braiding[i, j, k]` as a matrix whose size is determined by multiplicity.
We count the total number of matrix entries across all (i, j, k) with
N_{ij}^k > 0.

Multiplicity-free case: each nonzero fusion coefficient contributes 1.
"""
function _number_of_variables_in_hexagon_equations(poly_C)
    m = poly_C.simples
    mult = poly_C.tensor_product  # N_{ij}^k  (TensorCategories stores this as tensor_product)
    n_r = 0
    for i in 1:m, j in 1:m, k in 1:m
        N_ijk = mult[i, j, k]
        N_ijk == 0 && continue
        # Multiplicity-free: N_ijk is 0 or 1, contributes 1
        # General case: would contribute N_ijk^2 (full N×N block)
        n_r += N_ijk * N_ijk
    end
    return n_r
end

"""
    hexagon_equations(mult, one, F_values) -> (poly_C_fwd, poly_C_rev, eqs)

Build the hexagon equations with F-symbols fixed to numerical values. To
avoid symbolic matrix inversion (which fails on polynomial coefficients),
we treat the forward braiding `c` and the reverse braiding `c'` as
INDEPENDENT variable families, with three kinds of equations:

  (A) Left hexagon, using forward braiding c and forward associator α:
      (c_{X,Y} ⊗ id_Z) ∘ α_{Y,X,Z} ∘ (id_Y ⊗ c_{X,Z})
          = α_{Y,Z,X} ∘ c_{X, Y⊗Z} ∘ α_{X,Y,Z}

  (B) Right hexagon, rewritten with reverse braiding c' and inverse
      associator α⁻¹ (both as separate data):
      (id_Z ⊗ c'_{X,Y}) ∘ α'_{Z,X,Y} ∘ (c'_{X,Z} ⊗ id_Y)
          = α'_{X,Y,Z} ∘ c'_{X⊗Y, Z} ∘ α'_{Y,Z,X}
      where α' = α⁻¹ is precomputed NUMERICALLY (no symbolic inv).

  (C) Inverse consistency: c_{X,Y} ∘ c'_{X,Y} = id  (equivalently c ⋅ c' = 1
      for each multiplicity-free R-symbol component).

Returns the two categories (forward and reverse) and the combined equation
list. The polynomial ring has 2 * r_var_count variables (first half = r,
second half = s = reverse braiding).
"""
function hexagon_equations(mult::Array{Int,3}, one_vec::Vector{Int},
                           F_values::Vector{ComplexF64})
    m = size(mult, 1)

    # Sanity check F_values length
    _dummy = TensorCategories.six_j_category(QQ, mult)
    _dummy.one = one_vec
    expected_n_F = TensorCategories._number_of_variables_in_pentagon_equations(_dummy)
    @assert length(F_values) == expected_n_F "F_values length $(length(F_values)) != expected $expected_n_F"

    # Precompute F⁻¹ values numerically (block-wise)
    Finv_values = invert_associator_numeric(F_values, mult, one_vec)

    # Count R variables (using the dummy category for consistent structure)
    _dummy_cc = TensorCategories.six_j_category(AcbField(), mult)
    _dummy_cc.one = one_vec
    r_var_count = _number_of_variables_in_hexagon_equations(_dummy_cc)

    # Build a single polynomial ring with 2 * r_var_count variables:
    # first half -> forward braiding (r), second half -> reverse braiding (s)
    R_ring, xs = polynomial_ring(AcbField(), 2 * r_var_count)
    r_vars = xs[1:r_var_count]
    s_vars = xs[r_var_count + 1 : 2 * r_var_count]

    # Two categories sharing the same R_ring:
    # poly_C_fwd: associator = F, braiding = r
    # poly_C_rev: associator = F⁻¹, braiding = s
    poly_C_fwd = TensorCategories.six_j_category(R_ring, mult)
    poly_C_fwd.one = one_vec
    assign_F_to_associator!(poly_C_fwd, F_values)

    poly_C_rev = TensorCategories.six_j_category(R_ring, mult)
    poly_C_rev.one = one_vec
    assign_F_to_associator!(poly_C_rev, Finv_values)

    # Fill braidings
    function _fill_braiding!(poly_C, vars)
        m_ = poly_C.simples
        braid_arr = Array{MatElem, 3}(undef, m_, m_, m_)
        y = copy(vars)  # stack
        for i in 1:m_, j in 1:m_, k in 1:m_
            N_ijk = mult[i, j, k]
            if N_ijk == 0
                braid_arr[i, j, k] = zero_matrix(R_ring, 0, 0)
                continue
            end
            entries = [pop!(y) for _ in 1:(N_ijk * N_ijk)]
            braid_arr[i, j, k] = matrix(R_ring, N_ijk, N_ijk, entries)
        end
        @assert isempty(y)
        TensorCategories.set_braiding!(poly_C, braid_arr)
    end
    _fill_braiding!(poly_C_fwd, r_vars)
    _fill_braiding!(poly_C_rev, s_vars)

    # Write down the equations categorically — NO symbolic inv is used.
    eqs = elem_type(R_ring)[]
    Ss_fwd = simples(poly_C_fwd)
    Ss_rev = simples(poly_C_rev)

    # (A) Left hexagon in poly_C_fwd
    for X in Ss_fwd, Y in Ss_fwd, Z in Ss_fwd
        lhs = (braiding(X, Y) ⊗ id(Z)) ∘ associator(Y, X, Z) ∘ (id(Y) ⊗ braiding(X, Z))
        rhs = associator(Y, Z, X) ∘ braiding(X, Y ⊗ Z) ∘ associator(X, Y, Z)
        append!(eqs, collect(matrix(lhs - rhs))[:])
    end

    # (B) Left hexagon SHAPE in poly_C_rev (which has α' = α⁻¹ and c' = reverse braiding).
    # The RIGHT hexagon of a braided category, when you treat c' and α' as new forward
    # data, becomes a left hexagon in the "reverse" category. Concretely:
    #   (c'_{X,Y} ⊗ id_Z) ∘ α'_{Y,X,Z} ∘ (id_Y ⊗ c'_{X,Z})
    #       = α'_{Y,Z,X} ∘ c'_{X, Y⊗Z} ∘ α'_{X,Y,Z}
    # This is equivalent to the standard right hexagon with original α and c⁻¹.
    for X in Ss_rev, Y in Ss_rev, Z in Ss_rev
        lhs = (braiding(X, Y) ⊗ id(Z)) ∘ associator(Y, X, Z) ∘ (id(Y) ⊗ braiding(X, Z))
        rhs = associator(Y, Z, X) ∘ braiding(X, Y ⊗ Z) ∘ associator(X, Y, Z)
        append!(eqs, collect(matrix(lhs - rhs))[:])
    end

    # (C) Inverse consistency: c_{X,Y} ∘ c'_{X,Y} = id_{X⊗Y}
    # Using the forward-category objects, braiding(X,Y) goes X⊗Y -> Y⊗X.
    # For c ∘ c' to make sense as X⊗Y -> X⊗Y, we compose:
    #   (c'_{Y,X} from poly_C_rev) ∘ (c_{X,Y} from poly_C_fwd) : X⊗Y -> Y⊗X -> X⊗Y
    # but that requires moving between the two categories, which share R_ring
    # but are different SixJCategory instances. We instead impose the equation
    # DIRECTLY at the level of R-symbols and S-symbols:
    # For each (i, j, k) with N_{ij}^k > 0:
    #   R^{ij}_k · S^{ji}_k = 1   (multiplicity-free case; treating each as a scalar)
    # Or in general as matrix equation: R^{ij}_k · S^{ji}_k = I
    # Rationale: c_{X_i, X_j}: X_i ⊗ X_j -> X_j ⊗ X_i has component X_k equal to
    # R^{ij}_k; the reverse braiding c'_{X_j, X_i}: X_j ⊗ X_i -> X_i ⊗ X_j has
    # component X_k equal to S^{ji}_k. Their composite on the X_k component
    # must be identity.
    #
    # We need to index r_vars and s_vars consistently. Recall _fill_braiding!
    # assigns variables by popping; we reproduce that indexing here.
    function _var_index(i::Int, j::Int, k::Int, a::Int, b::Int,
                        var_pool_size::Int, m_::Int)
        # Compute the position of entry (a, b) of braid[i,j,k] inside the vars
        # used by _fill_braiding!. Since pop! consumes from the end, the LAST
        # assigned block is (m_, m_, m_) and within it entries[1] = pop = vars[end].
        # General formula: count how many (i', j', k', a', b') come AFTER (i,j,k,a,b)
        # in the assignment order (= pop order).
        # Total positions BEFORE = count of earlier (i',j',k') blocks plus earlier
        # (a',b') within current block.
        # Easier: enumerate and memoize.
        error("use the precomputed lookup `_rs_index` instead")
    end

    # Build lookup: (i, j, k) -> (start_pos, size) for the R variable indices.
    rs_block_info = Dict{Tuple{Int,Int,Int}, Tuple{Int,Int}}()
    let pos = r_var_count
        for i in 1:m, j in 1:m, k in 1:m
            N_ijk = mult[i, j, k]
            N_ijk == 0 && continue
            sq = N_ijk * N_ijk
            # In pop order, the LAST block iterated is (m, m, m). It receives the
            # LAST sq variables (positions pos-sq+1 .. pos).
            # We walk in normal order and need the positions that pop! would assign.
            # Since pop! consumes from the end, the FIRST block in our loop
            # receives the LAST sq variables of the pool. But our loop order is
            # i in 1:m, j in 1:m, k in 1:m (nested), which means the LAST block
            # visited is (m, m, m), which pops FIRST from the remaining stack,
            # i.e., it gets the highest position.
            # Wait — `pop!(y)` where y = copy(vars) with vars=xs[1:r_var_count].
            # y starts as [x1, x2, ..., x_{r_var_count}]. pop! returns x_{r_var_count} first.
            # Loop order visits (1,1,1), (1,1,2), ..., (m,m,m). So (1,1,1) block
            # gets x_{r_var_count}, x_{r_var_count - 1}, ... (first sq entries in pop order).
            # So block (1,1,1) occupies positions r_var_count - sq + 1 .. r_var_count (but in reverse).
            # Hmm — this is getting complex. Let me just do it carefully via simulation.
        end
    end

    # Simulate _fill_braiding! indexing exactly to build lookup tables.
    # For forward (r_vars) positions within r_vars (1..r_var_count):
    r_block_positions = Dict{Tuple{Int,Int,Int}, Vector{Int}}()
    let y_positions = collect(1:r_var_count)  # positions in r_vars
        for i in 1:m, j in 1:m, k in 1:m
            N_ijk = mult[i, j, k]
            N_ijk == 0 && continue
            # Each block pops N_ijk*N_ijk positions from the end of y_positions.
            # entries[p] = pop! = y_positions[end - (p-1)] (original), then list shrinks.
            popped = Int[]
            for _ in 1:(N_ijk * N_ijk)
                push!(popped, pop!(y_positions))
            end
            # popped[p] is the r_vars-index that becomes entries[p], which fills
            # matrix[a,b] with (a-1)*N_ijk + b = p, i.e., row-major.
            # We store the list so we can map (a,b) -> global index later.
            r_block_positions[(i, j, k)] = popped
        end
        @assert isempty(y_positions)
    end

    # Same structure for s_vars (indices within s_vars are 1..r_var_count;
    # global index in xs is r_var_count + local_index).
    s_block_positions = Dict{Tuple{Int,Int,Int}, Vector{Int}}()
    let y_positions = collect(1:r_var_count)
        for i in 1:m, j in 1:m, k in 1:m
            N_ijk = mult[i, j, k]
            N_ijk == 0 && continue
            popped = Int[]
            for _ in 1:(N_ijk * N_ijk)
                push!(popped, pop!(y_positions))
            end
            s_block_positions[(i, j, k)] = popped
        end
        @assert isempty(y_positions)
    end

    # Now add the inverse-consistency equations R^{ij}_k · S^{ji}_k = I.
    # For multiplicity-free (N_ijk ∈ {0, 1}) this is simply r * s - 1 = 0.
    # In general: matrix R^{ij}_k (N_ijk × N_ijk) times S^{ji}_k (= N_jik × N_jik
    # = N_ijk × N_ijk since N_ijk = N_jik for symmetric tensor product, but we
    # should be careful).
    for i in 1:m, j in 1:m, k in 1:m
        N_ijk = mult[i, j, k]
        N_ijk == 0 && continue

        # Forward block r-matrix R^{ij}_k: lives in r_vars at positions r_block_positions[(i,j,k)]
        r_pos = r_block_positions[(i, j, k)]
        R_block = [r_vars[r_pos[(a - 1) * N_ijk + b]] for a in 1:N_ijk, b in 1:N_ijk]

        # Reverse block s-matrix S^{ji}_k: lives in s_vars at positions s_block_positions[(j,i,k)]
        # (note the swap of i and j)
        if !haskey(s_block_positions, (j, i, k))
            continue  # fusion rule is not symmetric for this (i,j,k)? skip defensively
        end
        s_pos = s_block_positions[(j, i, k)]
        N_jik = mult[j, i, k]
        @assert N_jik == N_ijk "Fusion multiplicity not symmetric: N^$k_{$i,$j}=$N_ijk vs N^$k_{$j,$i}=$N_jik"
        S_block = [s_vars[s_pos[(a - 1) * N_jik + b]] for a in 1:N_jik, b in 1:N_jik]

        # Product R_block * S_block = I
        for a in 1:N_ijk, c in 1:N_ijk
            expr = sum(R_block[a, b] * S_block[b, c] for b in 1:N_ijk)
            target = (a == c) ? one(R_ring) : zero(R_ring)
            push!(eqs, expr - target)
        end
    end

    return poly_C_fwd, poly_C_rev, filter(e -> !iszero(e), unique(eqs))
end

"""
    get_hexagon_system(Nijk, r, F_values) -> (R_ring, eqs, n_vars)

Wrapper analogous to `get_pentagon_system`. Returns polynomial ring, equation
list, and number of variables (= 2 * r_var_count, since both forward and
reverse braidings are variables).
"""
function get_hexagon_system(Nijk::Array{Int,3}, r::Int, F_values::Vector{<:Number})
    one_vec = zeros(Int, r)
    one_vec[1] = 1
    # F_values may arrive as Vector{Complex{Float64}}, make sure it's ComplexF64
    F_values_cf64 = ComplexF64.(F_values)
    poly_C_fwd, poly_C_rev, eqs = hexagon_equations(Nijk, one_vec, F_values_cf64)
    eqs_filt = filter(eq -> !(eq isa Integer) && !iszero(eq), eqs)
    @assert !isempty(eqs_filt) "Hexagon produced no equations — system may be trivially satisfied"
    R_ring = parent(eqs_filt[1])
    n_vars = nvars(R_ring)
    return R_ring, eqs_filt, n_vars
end

"""
    solve_hexagon_homotopy(eqs, n; kwargs...) -> Vector{Vector{ComplexF64}}

Solve hexagon equations via HomotopyContinuation.jl. Hexagon has NO gauge
freedom (gauge was already exhausted by Pentagon), so the system should be
zero-dimensional with all solutions nonsingular.

Keyword arguments are the same as `solve_pentagon_homotopy`, minus `slice`
and `include_singular` (which are pentagon-specific concerns).
"""
function solve_hexagon_homotopy(eqs, n::Int;
                                certify_solutions::Bool=false,
                                threading::Bool=true,
                                start_system::Symbol=:polyhedral,
                                show_progress::Bool=false)
    # The hexagon polynomials have AcbField coefficients. For HC we need
    # ComplexF64 coefficients, so conversion happens in oscar_poly_to_hc_complex.
    sys, hc_vars = build_hc_system_complex(eqs, n)

    result = HC.solve(sys;
                      start_system=start_system,
                      threading=threading,
                      show_progress=show_progress)
    sols_raw = HC.solutions(result; only_nonsingular=true)

    if certify_solutions && !isempty(sols_raw)
        cert = HC.certify(sys, sols_raw)
        certs = HC.certificates(cert)
        sols_raw = [HC.solution_approximation(c) for c in certs if HC.is_certified(c)]
        println("  Certified: $(length(sols_raw)) / $(length(certs)) solutions")
    end

    nsols_total = HC.nsolutions(result)
    nsing = HC.nsingular(result)
    natinf = HC.nat_infinity(result)
    nfail = HC.nfailed(result)
    println("  HC hexagon: $(length(sols_raw)) nonsingular " *
            "(total=$nsols_total, singular=$nsing, at_infinity=$natinf, failed=$nfail)")
    return [ComplexF64.(s) for s in sols_raw]
end

"""
    oscar_poly_to_hc_complex(f, hc_vars) -> HC.Expression

Convert an Oscar polynomial with **AcbField (complex)** coefficients into a
HomotopyContinuation Expression. Coefficients are converted via ComplexF64.
"""
function oscar_poly_to_hc_complex(f, hc_vars::Vector{HC.Variable})
    iszero(f) && return HC.Expression(0)
    expr = HC.Expression(0)
    for (c, m) in zip(coefficients(f), monomials(f))
        # c is an AcbFieldElem; convert to ComplexF64 via real/imag midpoints
        cre = Float64(real(c))
        cim = Float64(imag(c))
        coef = ComplexF64(cre, cim)
        term = HC.Expression(coef)
        degs = degrees(m)
        for i in 1:length(degs)
            if degs[i] > 0
                term *= hc_vars[i]^degs[i]
            end
        end
        expr += term
    end
    expr
end

function build_hc_system_complex(eqs, n::Int)
    hc_vars = [HC.Variable(Symbol("r", i)) for i in 1:n]
    hc_exprs = HC.Expression[]
    for eq in eqs
        iszero(eq) && continue
        push!(hc_exprs, oscar_poly_to_hc_complex(eq, hc_vars))
    end
    return HC.System(hc_exprs; variables=hc_vars), hc_vars
end

# ============================================================
#  Main
# ============================================================

function classify_mtc(N::Int; max_rank::Int=20, max_trials::Int=20,
                      solver::Symbol=:homotopy,
                      certify_solutions::Bool=false,
                      refine::Bool=true,
                      slice::Int=0,
                      auto_slice::Bool=true,
                      include_singular::Bool=false,
                      solve_hexagon::Bool=true,
                      max_hexagon_per_candidate::Int=3)
    println("=" ^ 60)
    println("MTC Classification: N = $N  (pentagon solver = :$solver)")
    println("=" ^ 60)
    
    if N == 1
        println("N=1: trivial MTC (rank 1, F = 1)")
        println("=" ^ 60)
        return
    end
    
    candidates = enumerate_modular_data(N; max_rank=max_rank)
    
    for (idx, cand) in enumerate(candidates)
        println("\n--- Candidate $idx: rank=$(cand.rank) ---")
        
        if cand.rank == 1
            println("  Pointed (rank 1). F = trivial.")
            continue
        end
        
        R, eqs, n = get_pentagon_system(cand.Nijk, cand.rank)
        println("  Pentagon: $n variables, $(length(eqs)) equations")
        
        if n == 0
            println("  No free F-symbols. Trivial solution.")
            continue
        end
        
        zeta = exp(2π * im / N)
        solutions = if solver === :newton
            solve_pentagon_newton(eqs, n, zeta; max_trials=max_trials)
        elseif solver === :homotopy
            # First attempt with the user-specified slice count.
            sols = solve_pentagon_homotopy(eqs, n;
                                           slice=slice,
                                           include_singular=include_singular,
                                           certify_solutions=certify_solutions)
            # Auto slicing: if no nonsingular solutions came back and
            # include_singular is false, try progressively more slices.
            if auto_slice && isempty(sols) && !include_singular && slice == 0
                for s in 1:min(n-1, 5)
                    println("  [auto] retrying with slice=$s")
                    sols = solve_pentagon_homotopy(eqs, n;
                                                   slice=s,
                                                   certify_solutions=certify_solutions)
                    !isempty(sols) && break
                end
            end
            if refine && !isempty(sols)
                [refine_solution_newton(eqs, s) for s in sols]
            else
                sols
            end
        else
            error("Unknown solver :$solver (use :newton or :homotopy)")
        end
        println("  $(length(solutions)) pentagon solutions found")

        # Level III: Hexagon for each pentagon solution (up to a cap)
        if solve_hexagon && !isempty(solutions)
            println("\n  === Level III (Hexagon) ===")
            n_try = min(length(solutions), max_hexagon_per_candidate)
            for (ki, F_sol) in enumerate(solutions[1:n_try])
                println("\n  Pentagon solution $ki / $(length(solutions)):")
                try
                    R_ring, hex_eqs, n_r = get_hexagon_system(cand.Nijk, cand.rank, F_sol)
                    println("    Hexagon: $n_r R-variables, $(length(hex_eqs)) equations")
                    if n_r == 0
                        println("    No R-symbols free. Trivial braiding.")
                        continue
                    end
                    R_sols = solve_hexagon_homotopy(hex_eqs, n_r;
                                                   certify_solutions=certify_solutions)
                    println("    → $(length(R_sols)) braiding(s) found")
                catch e
                    println("    Hexagon failed: ", e)
                    # Do not abort — continue with other pentagon solutions
                end
            end
        end
    end
    
    println("\n" * "=" ^ 60)
    println("N=$N classification complete.")
    println("=" ^ 60)
end

load_sl2reps()
println("MTC Classifier (Level I exact, Level II+III: HomotopyContinuation) loaded.")
println("Usage: classify_mtc(N; solver=:homotopy|:newton, solve_hexagon=true)")
