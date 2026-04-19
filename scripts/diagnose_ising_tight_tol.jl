"""
scripts/diagnose_ising_tight_tol.jl

Decisive diagnostic: redo Newton+slice with much tighter tolerance and
more iterations. If the 4 perturbation solutions are genuine pentagon
solutions, they should stay; if they were numerical artifacts, they
should collapse onto the base F (= solution 1) or disappear.
"""

using LinearAlgebra
using Random
using ACMG
using ACMG.Phase4
using Oscar

const KC  = Phase4.KitaevComplex
const SPS = Phase4.SlicedPentagonSolver
const PS  = Phase4.PentagonSolver

function ising_Nijk()
    r = 3
    N = zeros(Int, r, r, r)
    for j in 1:r
        N[1, j, j] = 1
        N[j, 1, j] = 1
    end
    N[2, 2, 1] = 1
    N[2, 3, 3] = 1
    N[3, 2, 3] = 1
    N[3, 3, 1] = 1
    N[3, 3, 2] = 1
    return N
end

function ising_F_func()
    F_sss_s = (1.0 / sqrt(2.0)) * [1.0  1.0;
                                    1.0 -1.0]
    fr = FusionRule(ising_Nijk())
    function F(a::Int, b::Int, c::Int, d::Int, e::Int, f::Int)
        (fr.N[a, b, e] == 1 && fr.N[e, c, d] == 1 &&
         fr.N[b, c, f] == 1 && fr.N[a, f, d] == 1) || return 0.0
        1 in (a, b, c, d) && return 1.0
        nsig = count(==(3), (a, b, c, d))
        nsig == 0 && return 1.0
        if nsig == 2
            (a, b, c, d) == (2, 3, 2, 3) && return -1.0
            (a, b, c, d) == (3, 2, 3, 1) && return  1.0
            (a, b, c, d) == (3, 2, 3, 2) && return -1.0
            return 1.0
        end
        if nsig == 4
            d != 3 && return 0.0
            return F_sss_s[e, f]
        end
        return 1.0
    end
    return F
end

function main()
    Nijk = ising_Nijk()
    F_fn = ising_F_func()
    r = 3

    println("="^68)
    println("Test 1: tight tolerance, many iterations")
    println("="^68)
    Random.seed!(20260418)
    sols_tight = SPS.solve_pentagon_newton_with_slice(Nijk, r, F_fn;
        max_trials    = 5,
        max_iter      = 2000,
        perturb_scale = 0.05,
        tol           = 1e-14,
        verbose       = true)
    println("sols (tol 1e-14): $(length(sols_tight))")

    println("\n" * "="^68)
    println("Test 2: loose tolerance (test's original)")
    println("="^68)
    Random.seed!(20260418)
    sols_loose = SPS.solve_pentagon_newton_with_slice(Nijk, r, F_fn;
        max_trials    = 5,
        max_iter      = 200,
        perturb_scale = 0.05,
        tol           = 1e-10,
        verbose       = true)
    println("sols (tol 1e-10): $(length(sols_loose))")

    println("\n" * "="^68)
    println("Test 3: zero perturbation (same initial point)")
    println("="^68)
    Random.seed!(20260418)
    sols_zero = SPS.solve_pentagon_newton_with_slice(Nijk, r, F_fn;
        max_trials    = 5,
        max_iter      = 2000,
        perturb_scale = 0.0,
        tol           = 1e-14,
        verbose       = true)
    println("sols (perturb=0): $(length(sols_zero))  (should be 1)")

    println("\n" * "="^68)
    println("Test 4: moderate perturbation with tight tol")
    println("="^68)
    Random.seed!(20260418)
    sols_mod = SPS.solve_pentagon_newton_with_slice(Nijk, r, F_fn;
        max_trials    = 10,
        max_iter      = 2000,
        perturb_scale = 0.3,
        tol           = 1e-14,
        verbose       = true)
    println("sols (perturb=0.3, tol 1e-14): $(length(sols_mod))")

    # Diagnostics on sols_mod
    R, eqs, n = Phase4.get_pentagon_system(Nijk, r)
    one_vec = [1, 0, 0]
    fkey_map = SPS.build_fkey_to_xvar_map(Nijk, r, one_vec)
    F_base = zeros(ComplexF64, n)
    for (key, pidx) in fkey_map
        F_base[pidx] = ComplexF64(F_fn(key...))
    end

    println("\nPentagon residuals at each sols_mod solution (at tol 1e-14):")
    for (i, s) in enumerate(sols_mod)
        F_pent = ComplexF64[PS.eval_poly_complex(eq, s) for eq in eqs]
        pent_res = maximum(abs.(F_pent))
        dist = norm(s - F_base)
        println("  sol $i: pent_res=$pent_res,  ‖s − F_base‖=$dist")
    end
end

main()
