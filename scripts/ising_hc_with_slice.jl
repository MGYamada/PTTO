"""
scripts/ising_hc_with_slice.jl

Benchmark the Ising pentagon HC solve with and without the Kitaev/Hodge
gauge slice.

Reports:
- Mixed volume of the augmented polynomial system (cheap upper bound)
- Time to solve via HC polyhedral homotopy (may be long WITHOUT slice)
- Number of non-singular solutions
- Residuals at the returned solutions

Usage:
    julia --project=. scripts/ising_hc_with_slice.jl [--with-unsliced]

By default only the sliced version is run (it's the interesting case).
Pass `--with-unsliced` to also attempt the unsliced HC (may take hours
or fail; the purpose is to see what mixed volume HC reports for baseline).
"""

using LinearAlgebra
using ACMG
using ACMG.Phase4
import HomotopyContinuation as HC

const KC  = Phase4.KitaevComplex
const SPS = Phase4.SlicedPentagonSolver

# ---------------- Ising setup (same as tests) ----------------

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

# ---------------- Benchmark helpers ----------------

function report_system(label, sys, hc_vars)
    println("="^68)
    println(label)
    println("="^68)
    n_vars = length(hc_vars)
    n_eqs = length(HC.expressions(sys))
    println("  variables           : $n_vars")
    println("  equations           : $n_eqs")
    # Polyhedral (mixed volume) bound — only meaningful for square systems.
    # For over-determined systems we still attempt the path count estimate.
    if n_eqs == n_vars
        try
            mv = HC.mixed_volume(sys)
            println("  mixed volume        : $mv")
        catch e
            println("  mixed volume        : (not computed; $e)")
        end
    else
        println("  mixed volume        : (system is $n_eqs×$n_vars, not square)")
        # For over-determined systems, we can still pass to HC.solve, which
        # will internally slice the system with random linear equations to
        # make it square, or report nonsquare behaviour.
    end
end

function bench_solve(label, sys; show_progress::Bool = true,
                                 start_system::Symbol = :polyhedral)
    println("-"^68)
    println("Solving: $label")
    println("-"^68)
    t0 = time()
    local result
    try
        result = HC.solve(sys;
                          start_system = start_system,
                          show_progress = show_progress,
                          threading = true)
    catch e
        elapsed = time() - t0
        println("  SOLVE FAILED after $(round(elapsed, digits=2))s: $e")
        return nothing
    end
    elapsed = time() - t0
    println("  time                : $(round(elapsed, digits=2))s")
    println("  result              : $result")
    sols = HC.solutions(result)
    println("  solutions found     : $(length(sols))")
    real_sols = HC.real_solutions(result)
    println("  real solutions      : $(length(real_sols))")
    return result
end

# ---------------- Main ----------------

function main()
    with_unsliced = "--with-unsliced" in ARGS

    Nijk = ising_Nijk()
    F_fn = ising_F_func()

    println("\n" * "="^68)
    println("Ising pentagon HC benchmark")
    println("="^68)

    # ----- Sliced system -----
    res = SPS.get_sliced_pentagon_system_hc(Nijk, 3, F_fn)
    report_system("Ising pentagon WITH slice", res.sys, res.hc_vars)
    sliced_result = bench_solve("sliced Ising pentagon", res.sys)

    if sliced_result !== nothing
        # Evaluate pentagon residual at each returned solution (only the first
        # 14 equations = pentagon, last 1 = slice).
        sols = HC.solutions(sliced_result)
        n_eqs_total = length(HC.expressions(res.sys))
        n_slice = res.n_slice
        n_pent = n_eqs_total - n_slice
        println("\n  Residuals at returned solutions (pentagon only):")
        for (k, s) in enumerate(sols[1:min(end, 10)])
            all_res = res.sys(s)
            pent_res = all_res[1:n_pent]
            println("    sol $k: pentagon max = $(maximum(abs.(pent_res)))")
        end
    end

    # ----- Unsliced system (optional, may take hours) -----
    if with_unsliced
        println("\n" * "="^68)
        println("BASELINE: Ising pentagon WITHOUT slice")
        println("="^68)
        R, eqs, n = Phase4.get_pentagon_system(Nijk, 3)
        # Build HC system from Oscar pentagon eqs only
        hc_vars = [HC.Variable(Symbol("x", i)) for i in 1:n]
        hc_exprs = HC.Expression[]
        for eq in eqs
            (isa(eq, Integer) || iszero(eq)) && continue
            push!(hc_exprs, Phase4.PentagonSolver.oscar_poly_to_hc(eq, hc_vars))
        end
        sys_unsliced = HC.System(hc_exprs; variables = hc_vars)
        report_system("Ising pentagon WITHOUT slice", sys_unsliced, hc_vars)
        bench_solve("unsliced Ising pentagon", sys_unsliced)
    else
        println("\n(Run with --with-unsliced to also attempt the baseline.)")
    end
end

main()
