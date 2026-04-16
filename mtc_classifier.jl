# mtc_classifier.jl
# MTC Classification via Conductor Decomposition
# April 14-15, 2026 — Arithmetic Condensed Matter Geometry
#
# Level I (exact over Q(ζ_N)):
#   N = ∏ pᵢ^aᵢ → SL₂(Z/pᵢ^aᵢ Z) irreps → tensor products → modular data
#   Admissibility checks: Verlinde, (ST)³, Cauchy (via norm), Frobenius-Schur
# Level II (numerical, ComplexF64):
#   Pentagon — damped Newton with KrylovKit sparse solver
#
# Usage:
#   include("mtc_classifier.jl")
#   classify_mtc(5)
#   classify_mtc(9)

using LinearAlgebra
using SparseArrays
using Oscar
using TensorCategories
using KrylovKit

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

# --- Admissibility: (ST)³ = ±S² in Q(ζ_N) ---

function check_ST_exact(S, T, K, r::Int)
    """Check (ST)³ = ±S² over Q(ζ_N). Projective SL₂(Z) representation."""
    ST = S * T
    ST3 = ST * ST * ST
    S2 = S * S
    
    # Check ST3 == S2 or ST3 == -S2
    diff_plus = ST3 - S2
    diff_minus = ST3 + S2
    
    iszero(diff_plus) && return true
    iszero(diff_minus) && return true
    
    # Also try scalar multiples ST3 = λ S2 for λ a root of unity
    # For simplicity, only check ±1
    return false
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

function enumerate_modular_data(N::Int; max_rank::Int=20, verbose::Bool=false)
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
    
    # Step 2: Tensor products
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
        
        # (ST)³ = ±S²
        if !check_ST_exact(S_fixed, T, K, r)
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
            Nijk=Nijk, rank=r, name=name
        ))
    end
    
    if length(factor_irreps) == 1
        for (S, T, r) in factor_irreps[1]
            process_candidate(S, T, r)
        end
    else
        indices = [1:length(fi) for fi in factor_irreps]
        for idx_tuple in Iterators.product(indices...)
            total_rank = prod(factor_irreps[k][idx_tuple[k]][3] for k in 1:length(pf))
            total_rank > max_rank && continue
            
            # Tensor product
            S = factor_irreps[1][idx_tuple[1]][1]
            T = factor_irreps[1][idx_tuple[1]][2]
            for k in 2:length(factor_irreps)
                S = kron_oscar(S, factor_irreps[k][idx_tuple[k]][1], K)
                T = kron_oscar(T, factor_irreps[k][idx_tuple[k]][2], K)
            end
            
            name = join(["$(factor_irreps[k][idx_tuple[k]][3])d" for k in 1:length(pf)], "⊗")
            process_candidate(S, T, total_rank, name)
        end
    end
    
    println("\nRejections: Verlinde=$(reject_count["Verlinde"]), " *
            "(ST)³=$(reject_count["(ST)³"]), " *
            "Cauchy=$(reject_count["Cauchy"]), " *
            "FS=$(reject_count["FS"])")
    println("$(length(valid_candidates)) valid modular data candidates.")
    for (i, cand) in enumerate(valid_candidates)
        name = hasproperty(cand, :name) ? cand.name : ""
        println("  [$i] rank=$(cand.rank) $name")
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
#  Main
# ============================================================

function classify_mtc(N::Int; max_rank::Int=20, max_trials::Int=20)
    println("=" ^ 60)
    println("MTC Classification: N = $N")
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
        solutions = solve_pentagon_newton(eqs, n, zeta; max_trials=max_trials)
        println("  $(length(solutions)) solutions found")
    end
    
    println("\n" * "=" ^ 60)
    println("N=$N classification complete.")
    println("=" ^ 60)
end

load_sl2reps()
println("MTC Classifier (exact Level I) loaded.")
println("Usage: classify_mtc(N)")
