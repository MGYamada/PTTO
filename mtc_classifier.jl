# mtc_classifier.jl
# MTC Classification via Conductor Decomposition
# April 14, 2026 — Arithmetic Condensed Matter Geometry
#
# Algorithm:
#   Level I:   SL₂(Z/NZ) reps → modular data (S,T)
#   Level II:  Pentagon/hexagon — Newton's method
#
# Usage:
#   include("mtc_classifier.jl")
#   classify_mtc(5)

using LinearAlgebra
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
#  Level II: Pentagon solver
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

function newton_damped(eqs, n, zeta; max_trials=20, max_iter=200)
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
                println("Trial $trial: CONVERGED at iter $iter")
                push!(solutions, copy(x))
                break
            end
            
            J = jacobian_complex(eqs, x, n)
            delta = pinv(J) * F_val
            
            # Line search: find best step size
            alpha = 1.0
            for _ in 1:20
                x_new = x - alpha * delta
                F_new = ComplexF64[eval_poly_complex(eq, x_new) for eq in eqs]
                res_new = maximum(abs.(F_new))
                if res_new < res
                    x = x_new
                    break
                end
                alpha *= 0.5
            end
        end
    end
    solutions
end

function eval_poly_complex(f, vals::Vector{ComplexF64})
    # f が定数 0 の場合
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

function jacobian_complex(eqs, vals::Vector{ComplexF64}, n)
    m = length(eqs)
    J = zeros(ComplexF64, m, n)
    for i in 1:m
        isa(eqs[i], Integer) && continue
        iszero(eqs[i]) && continue
        for j in 1:n
            df = derivative(eqs[i], j)
            J[i,j] = eval_poly_complex(df, vals)
        end
    end
    J
end

# --- Hexagon solver ---

function solve_hexagon_Fp(F_values, Nijk, r, Fp)
    # With F known, hexagon is linear in R. Solve directly.
    # TODO: implement hexagon equations
    Dict{NTuple{3,Int}, Any}()
end

# ============================================================
#  Main
# ============================================================

function classify_mtc(N::Int; max_rank::Int=20)
    println("=" ^ 60)
    println("MTC Classification: N = $N")
    println("=" ^ 60)

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

    # Level II: Solve pentagon
    for (idx, cand) in enumerate(valid_candidates)
        println("\n--- Candidate $idx: rank=$(cand.rank), level=$(cand.level) ---")

        if cand.rank == 1
            println("  Pointed category. F = trivial.")
            continue
        end

        R, eqs, n = get_pentagon_system(cand.Nijk, cand.rank)
        println("Pentagon: $n variables, $(length(eqs)) equations")

        solutions = newton_damped(eqs, n, exp(2π * im / N))
        println("$(length(solutions)) solutions")
        for (si, sol) in enumerate(solutions)
            println("    [$si] $sol")
        end
    end

    println("\n" * "=" ^ 60)
    println("Done. N=$N classified.")
    println("=" ^ 60)
end

load_sl2reps()
println("MTC Classifier loaded. Usage: classify_mtc(N)")
