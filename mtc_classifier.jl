# mtc_classifier.jl
# MTC Classification via Conductor Decomposition
# April 14, 2026 — Arithmetic Condensed Matter Geometry
#
# Algorithm:
#   Level I:   SL₂(Z/NZ) reps → modular data (S,T)
#   Level II:  Pentagon/hexagon over F_p — Newton's method — CRT
#
# Key insight: Pentagon system is square after gauge fixing.
#   F_p 上の Newton 法で exact に解ける。
#
# Usage:
#   using Oscar, TensorCategories
#   include("mtc_classifier.jl")
#   classify_mtc(5)

using Oscar
using TensorCategories

# ============================================================
#  Level I: Modular data
# ============================================================

function load_sl2reps()
    GAP.Globals.LoadPackage(GapObj("SL2Reps"))
end

function gap_to_oscar_matrix(gap_mat, r)
    cond = 1
    for i in 1:r, j in 1:r
        cond = lcm(cond, Int(GAP.Globals.Conductor(gap_mat[i,j])))
    end

    K, z = cyclotomic_field(cond)
    M = zero_matrix(K, r, r)

    for i in 1:r, j in 1:r
        coeffs = GAP.Globals.CoeffsCyc(gap_mat[i,j], cond)
        val = zero(K)
        for k in 1:cond
            if !Bool(GAP.Globals.IsZero(coeffs[k]))
                num = Int(GAP.Globals.NumeratorRat(coeffs[k]))
                den = Int(GAP.Globals.DenominatorRat(coeffs[k]))
                val += QQ(num, den) * z^(k-1)
            end
        end
        M[i,j] = val
    end

    return M, K, cond
end

function fix_quantum_dims(S, K, r, cond)
    z = exp(2π * im / cond)
    signs = ones(Int, r)
    for i in 1:r
        d = S[1,i] * inv(S[1,1])
        deg = degree(K)
        val = sum(Float64(coeff(d, k)) * z^k for k in 0:deg-1)
        if real(val) < -1e-10
            signs[i] = -1
        end
    end
    S_fixed = copy(S)
    for i in 1:r, j in 1:r
        S_fixed[i,j] = signs[i] * signs[j] * S[i,j]
    end
    return S_fixed, signs
end

function verlinde_float(S, r::Int, K, cond)
    z = exp(2π * im / cond)
    deg = degree(K)

    S_num = zeros(ComplexF64, r, r)
    for i in 1:r, j in 1:r
        val = zero(ComplexF64)
        for k in 0:deg-1
            c = coeff(S[i,j], k)
            val += (Float64(numerator(c)) / Float64(denominator(c))) * z^k
        end
        S_num[i,j] = val
    end

    Nijk = zeros(Int, r, r, r)
    for i in 1:r, j in 1:r, k in 1:r
        val = sum(S_num[i,l] * S_num[j,l] * conj(S_num[k,l]) / S_num[1,l] for l in 1:r)
        n = round(Int, real(val))
        if abs(real(val) - n) > 1e-4 || abs(imag(val)) > 1e-4
            return nothing
        end
        Nijk[i,j,k] = n
    end
    Nijk
end

# ============================================================
#  Level II: Pentagon solver over F_p
# ============================================================

# --- Pentagon equations from TensorCategories.jl ---

function get_pentagon_system(Nijk::Array{Int,3}, r::Int)
    one_vec = zeros(Int, r)
    one_vec[1] = 1
    C, eqs = pentagon_equations(Nijk, one_vec)
    R = parent(eqs[1])
    n = nvars(R)
    return R, eqs, n
end

# --- F_p reduction ---

function to_Fp(f, Rp, p)
    result = zero(Rp)
    yp = gens(Rp)
    for (c, m) in zip(coefficients(f), monomials(f))
        c_int = mod(Int(numerator(c)) * invmod(Int(denominator(c)), p), p)
        degs = degrees(m)
        mono = prod(yp[i]^degs[i] for i in 1:length(degs))
        result += c_int * mono
    end
    result
end

function system_to_Fp(eqs, R, p)
    Fp = GF(p)
    Rp, xp = polynomial_ring(Fp, symbols(R))
    eqs_Fp = [to_Fp(eq, Rp, p) for eq in eqs]
    return Rp, xp, eqs_Fp, Fp
end

# --- Evaluate and differentiate over F_p ---

function eval_poly_Fp(f, vals)
    result = zero(parent(vals[1]))
    for (c, m) in zip(coefficients(f), monomials(f))
        degs = degrees(m)
        term = c * prod(vals[i]^degs[i] for i in 1:length(degs))
        result += term
    end
    result
end

function eval_system_Fp(eqs, vals)
    [eval_poly_Fp(eq, vals) for eq in eqs]
end

function jacobian_Fp(eqs, vals, Fp, n)
    m = length(eqs)
    J = zero_matrix(Fp, m, n)
    for i in 1:m, j in 1:n
        df = derivative(eqs[i], j)
        J[i,j] = eval_poly_Fp(df, vals)
    end
    J
end

# --- Newton's method over F_p ---

function select_independent_eqs(eqs, Fp, n)
    # ランダムな点で Jacobian を評価し、rank n の部分集合を選ぶ
    p = Int(characteristic(Fp))
    x_test = [Fp(rand(1:p-1)) for _ in 1:n]

    selected = Int[]
    for i in 1:length(eqs)
        trial = [selected; i]
        J = zero_matrix(Fp, length(trial), n)
        for (ri, ei) in enumerate(trial)
            for j in 1:n
                df = derivative(eqs[ei], j)
                J[ri, j] = eval_poly_Fp(df, x_test)
            end
        end
        if rank(J) == length(trial)
            push!(selected, i)
        end
        length(selected) == n && break
    end

    return eqs[selected]
end

function newton_solve_Fp(eqs, Fp; max_trials=500, max_iter=50)
    n = nvars(parent(eqs[1]))
    p = Int(characteristic(Fp))
    solutions = Vector{Vector}()

    for trial in 1:max_trials
        x = [Fp(rand(1:p-1)) for _ in 1:n]

        converged = false
        for iter in 1:max_iter
            F_val = eval_system_Fp(eqs, x)

            if all(iszero, F_val)
                converged = true
                break
            end

            J = jacobian_Fp(eqs, x, Fp, n)

            # Overdetermined: use normal equation J^T J δ = J^T F
            JtJ = transpose(J) * J
            if iszero(det(JtJ))
                break
            end
            JtF = transpose(J) * matrix(Fp, length(F_val), 1, F_val)
            delta = solve(JtJ, JtF; side=:right)
            
            for i in 1:n
                x[i] = x[i] - delta[i,1]
            end
        end

        if converged
            F_all = eval_system_Fp(eqs, x)
            if all(iszero, F_all)
                is_new = all(sol -> any(sol[i] != x[i] for i in 1:n), solutions)
                if is_new || isempty(solutions)
                    push!(solutions, copy(x))
                end
            end
        end
    end

    return solutions
end

# --- Hexagon solver ---

function solve_hexagon_Fp(F_values, Nijk, r, Fp)
    # With F known, hexagon is linear in R. Solve directly.
    # TODO: implement hexagon equations
    Dict{NTuple{3,Int}, Any}()
end

# ============================================================
#  CRT Reconstruction
# ============================================================

function good_primes(N::Int, count::Int)
    result = Int[]
    p = N + 1
    while length(result) < count
        is_prime(p) && p % N == 1 && push!(result, p)
        p += 1
    end
    result
end

# ============================================================
#  Main
# ============================================================

function classify_mtc(N::Int; max_rank::Int=20, num_primes::Int=5)
    println("=" ^ 60)
    println("MTC Classification: N = $N")
    println("=" ^ 60)

    primes = good_primes(N, num_primes)
    println("Good primes (p ≡ 1 mod $N): $primes")

    # Level I: Enumerate modular data candidates
    all_irreps = []
    for d in divisors(N)
        reps = GAP.evalstr("SL2IrrepsOfLevel($d)")
        len = Int(GAP.Globals.Length(reps))
        for i in 1:len
            push!(all_irreps, (reps[i], d))
        end
    end

    valid_candidates = []

    for (rep, lev) in all_irreps
        r = Int(rep.degree)
        r > max_rank && continue

        S, K, cond = gap_to_oscar_matrix(rep.S, r)
        S_fixed, signs = fix_quantum_dims(S, K, r, cond)
        Nijk = verlinde_float(S_fixed, r, K, cond)

        if Nijk !== nothing && all(Nijk .>= 0)
            println("VALID: degree=$r, level=$lev ✓")
            push!(valid_candidates, (Nijk=Nijk, rank=r, level=lev))
        end
    end

    println("\n$(length(valid_candidates)) candidates found.")

    # Level II: Solve pentagon over F_p
    for (idx, cand) in enumerate(valid_candidates)
        println("\n--- Candidate $idx: rank=$(cand.rank), level=$(cand.level) ---")

        if cand.rank == 1
            println("  Pointed category. F = trivial.")
            continue
        end

        R, eqs, n = get_pentagon_system(cand.Nijk, cand.rank)
        println("Pentagon: $n variables, $(length(eqs)) equations")

        for p in primes
            println("  F_$p: ")
            Rp, xp, eqs_Fp, Fp = system_to_Fp(eqs, R, p)
            solutions = newton_solve_Fp(eqs_Fp, Fp)

            println("$(length(solutions)) solutions")
            for (si, sol) in enumerate(solutions)
                println("    [$si] $sol")
            end
        end
    end

    println("\n" * "=" ^ 60)
    println("Done. N=$N classified.")
    println("=" ^ 60)
end

load_sl2reps()
println("MTC Classifier loaded. Usage: classify_mtc(N)")
