"""
scripts/diagnose_ising_newton_solutions.jl

Diagnose what the 5 solutions from `solve_pentagon_newton_with_slice`
on Ising actually are.

For each solution s in sols:
  1. Pentagon residual at s (sanity: should be < tol)
  2. Slice residual at s (sanity: should be ~ 0)
  3. Key F-entries:
        F^{σσσ}_{σ; 0,0}, F^{σσσ}_{σ; 0,1}, F^{σσσ}_{σ; 1,0}, F^{σσσ}_{σ; 1,1}
        F^{ψσψ}_σ        = F^{(2,3,2,3; 2,2)} — the 1-element block
        F^{σψσ}_ψ        = F^{(3,2,3,2; 2,2)}
        F^{σψσ}_1        = F^{(3,2,3,1; 2,2)}
  4. Pairwise differences ‖s_i − s_j‖

Compare against the base Ising F to see which solutions are close to it
(gauge-equivalent) vs far (potentially different F-class).

Usage:
    julia --project=. scripts/diagnose_ising_newton_solutions.jl
"""

using LinearAlgebra
using Random
using ACMG
using ACMG.Phase4
using Oscar

const KC  = Phase4.KitaevComplex
const SPS = Phase4.SlicedPentagonSolver
const PS  = Phase4.PentagonSolver

# ---------- Ising setup ----------

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
    Nijk  = ising_Nijk()
    F_fn  = ising_F_func()
    r     = 3

    # Re-run the Newton+slice with the same settings as the test,
    # but with verbose=true and collect all solutions.
    Random.seed!(20260418)

    sols = SPS.solve_pentagon_newton_with_slice(Nijk, r, F_fn;
                                                 max_trials   = 5,
                                                 max_iter     = 200,
                                                 perturb_scale = 0.05,
                                                 tol          = 1e-10,
                                                 verbose      = true)

    println("\nSolutions found: $(length(sols))")

    # Pentagon eqs + fkey map for residual evaluation
    R, eqs, n = Phase4.get_pentagon_system(Nijk, r)
    one_vec = [1, 0, 0]
    fkey_map = SPS.build_fkey_to_xvar_map(Nijk, r, one_vec)
    fr = FusionRule(Nijk)
    fcs = KC.build_F_coord_space(fr, F_fn)
    ga  = KC.analyze_gauge(fcs)
    A_slice, b_slice = SPS.build_slice_linear_system(ga, fkey_map, n)

    # Base F as ordered vector
    F_base = zeros(ComplexF64, n)
    for (key, pidx) in fkey_map
        F_base[pidx] = ComplexF64(F_fn(key...))
    end

    # Keys of interest
    key_sss = [(3, 3, 3, 3, e, f) for e in 1:2, f in 1:2]   # F^{σσσ}_σ 2×2 block
    key_2323 = (2, 3, 2, 3, 2, 2)  # F^{ψσψ}_σ
    key_3232 = (3, 2, 3, 2, 2, 2)  # F^{σψσ}_ψ
    key_3231 = (3, 2, 3, 1, 2, 2)  # F^{σψσ}_1

    println("\n" * "="^68)
    println("Base Ising F-symbol (reference)")
    println("="^68)
    println("  F^{σσσ}_{σ;e,f} block:")
    for e in 1:2, f in 1:2
        k = (3, 3, 3, 3, e, f)
        if haskey(fkey_map, k)
            p = fkey_map[k]
            println("    (e=$e,f=$f)  F = $(F_base[p])   x_$p")
        end
    end
    for (label, k) in (("F^{ψσψ}_σ", key_2323),
                       ("F^{σψσ}_ψ", key_3232),
                       ("F^{σψσ}_1", key_3231))
        if haskey(fkey_map, k)
            p = fkey_map[k]
            println("  $label = $(F_base[p])   x_$p")
        end
    end

    println("\n" * "="^68)
    println("Solution diagnostics")
    println("="^68)

    for (i, s) in enumerate(sols)
        println("\n--- Solution $i ---")

        # Pentagon residual
        F_pent = ComplexF64[PS.eval_poly_complex(eq, s) for eq in eqs]
        pent_res = maximum(abs.(F_pent))
        println("  pentagon residual max : $pent_res")

        # Slice residual
        slice_res = maximum(abs.(A_slice * s - b_slice))
        println("  slice residual max    : $slice_res")

        # Distance from base F
        println("  ‖s − F_base‖          : $(norm(s - F_base))")

        # σσσσ block
        println("  F^{σσσ}_{σ;e,f} block:")
        for e in 1:2, f in 1:2
            k = (3, 3, 3, 3, e, f)
            if haskey(fkey_map, k)
                p = fkey_map[k]
                print("    (e=$e,f=$f) = $(round(s[p], digits=5))")
                print("   (base was $(round(F_base[p], digits=5)))")
                println()
            end
        end

        # Key 2-σ entries
        for (label, k) in (("F^{ψσψ}_σ", key_2323),
                           ("F^{σψσ}_ψ", key_3232),
                           ("F^{σψσ}_1", key_3231))
            if haskey(fkey_map, k)
                p = fkey_map[k]
                print("  $label = $(round(s[p], digits=5))")
                println("   (base $(round(F_base[p], digits=5)))")
            end
        end
    end

    # Pairwise distances
    println("\n" * "="^68)
    println("Pairwise ‖s_i − s_j‖ (should exceed ~0.05 if genuinely distinct)")
    println("="^68)
    for i in 1:length(sols), j in (i+1):length(sols)
        d = norm(sols[i] - sols[j])
        println("  s_$i ↔ s_$j : $d")
    end
end

main()
