"""
scripts/inspect_ising_distant_solutions.jl

Inspect the 3 genuine pentagon solutions found by Newton+slice with
perturb_scale=0.3 on Ising:
  sol 1 = base F (reference)
  sol 2 = dist 0.61 from base
  sol 3 = dist 0.41 from base

For each solution, report:
  1. All 14 F-values side-by-side with base
  2. F^{σσσ}_σ as a 2×2 matrix with its determinant (gauge invariant)
  3. Characteristic polynomial of F^{σσσ}_σ block (gauge invariant)
  4. Which sign of Hadamard it corresponds to

Ising gauge-invariant quantities that should distinguish F-classes:
  - det(F^{σσσ}_σ) = ±(sign depending on the Hadamard orientation)
  - tr(F^{σσσ}_σ) = (gauge invariant under unit-preserving vertex gauge)
  - The signs {F^{ψσψ}_σ, F^{σψσ}_ψ} as a set of ±1
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

    # Reproduce Test 4 of the previous script
    Random.seed!(20260418)
    sols = SPS.solve_pentagon_newton_with_slice(Nijk, r, F_fn;
        max_trials    = 10,
        max_iter      = 2000,
        perturb_scale = 0.3,
        tol           = 1e-14,
        verbose       = false)
    println("Found $(length(sols)) solutions.")

    # Setup for displaying
    one_vec = [1, 0, 0]
    fkey_map = SPS.build_fkey_to_xvar_map(Nijk, r, one_vec)
    n = length(fkey_map)
    F_base = zeros(ComplexF64, n)
    for (key, pidx) in fkey_map
        F_base[pidx] = ComplexF64(F_fn(key...))
    end

    # Build an inverse map: pentagon index -> fkey
    pidx_to_key = Dict(v => k for (k, v) in fkey_map)

    println("\n" * "="^78)
    println("All F-values for every solution (side-by-side)")
    println("="^78)

    # Column headers
    println(rpad("idx", 4), rpad("F-key (a,b,c,d;e,f)", 25),
            rpad("base", 20),
            rpad("sol 2", 20),
            rpad("sol 3", 20))

    # Sort by pentagon index
    for p in 1:n
        key = pidx_to_key[p]
        key_str = "($(key[1]),$(key[2]),$(key[3]),$(key[4]);$(key[5]),$(key[6]))"
        print(rpad(p, 4), rpad(key_str, 25))
        print(rpad(round(F_base[p], digits=5), 20))
        if length(sols) ≥ 2
            print(rpad(round(sols[2][p], digits=5), 20))
        end
        if length(sols) ≥ 3
            print(rpad(round(sols[3][p], digits=5), 20))
        end
        println()
    end

    # Extract F^{σσσ}_σ blocks and gauge invariants
    println("\n" * "="^78)
    println("F^{σσσ}_σ block (2×2 matrix): gauge invariants")
    println("="^78)

    function extract_sss_block(F_vec)
        M = zeros(ComplexF64, 2, 2)
        for e in 1:2, f in 1:2
            k = (3, 3, 3, 3, e, f)
            if haskey(fkey_map, k)
                M[e, f] = F_vec[fkey_map[k]]
            end
        end
        return M
    end

    for (i, s) in enumerate(sols)
        M = extract_sss_block(s)
        println("\n--- Solution $i ---")
        println("  matrix:")
        for e in 1:2
            print("    ")
            for f in 1:2
                print(rpad(round(M[e, f], digits=5), 22))
            end
            println()
        end
        println("  trace (gauge invariant): $(round(tr(M), digits=6))")
        println("  det   (gauge invariant): $(round(det(M), digits=6))")
        # eigenvalues
        evs = eigvals(M)
        println("  eigenvalues: $(round.(evs, digits=5))")
        # M^2 (Ising should satisfy M^2 = I for Hadamard)
        M2 = M * M
        err_id = opnorm(M2 - I)
        println("  ‖M² - I‖: $err_id    (Ising standard: M² = I)")
    end

    println("\n" * "="^78)
    println("2-σ entries (Ising has these as ±1)")
    println("="^78)
    for (label, k) in (("F^{ψσψ}_σ         ", (2, 3, 2, 3, 2, 2)),
                       ("F^{σψσ}_ψ         ", (3, 2, 3, 2, 2, 2)),
                       ("F^{σψσ}_1         ", (3, 2, 3, 1, 2, 2)),
                       ("F^{ψψσ}_σ;1,2    ", (2, 2, 3, 3, 1, 2)),
                       ("F^{σσψ}_ψ;2,1    ", (3, 3, 2, 2, 2, 1)))
        if haskey(fkey_map, k)
            p = fkey_map[k]
            print(rpad(label, 20))
            for s in sols
                print(rpad(round(s[p], digits=5), 22))
            end
            println()
        end
    end

    println("\n" * "="^78)
    println("Gauge-equivalence test: try to find vertex gauge that maps sol → base")
    println("="^78)
    println("For Ising with unit-preserving vertex gauge, pentagon solutions in the")
    println("same gauge class must share:")
    println("  • det(F^{σσσ}_σ)")
    println("  • tr(F^{σσσ}_σ) modulo reordering e/f swaps")
    println("  • signs of the 1×1 blocks {F^{ψσψ}_σ, F^{σψσ}_ψ, F^{σψσ}_1}")
end

main()
