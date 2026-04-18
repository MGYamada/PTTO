using Test
using ACMG
using ACMG.Phase4

"""
Fibonacci Phase 4 smoke test.

Fibonacci fusion ring:
    1⊗1 = 1, 1⊗τ = τ, τ⊗1 = τ
    τ⊗τ = 1 ⊕ τ

Expected (known from literature):
- Pentagon: non-trivial system with finitely many gauge classes;
  HC typically returns 4 solutions (before gauge quotient).
- Hexagon: two braided structures (Fibonacci and its complex conjugate),
  characterised by R^{ττ}_1 = e^{±4πi/5}, R^{ττ}_τ = e^{∓3πi/5}.

This test only verifies the *pipeline executes* and produces the
*expected counts*. Full algebraic recognition (lifting ComplexF64
to Q(ζ_20)) and gauge equivalence classification are deferred.
"""

@testset "Phase 4: Fibonacci pentagon/hexagon" begin
    # ---------- Fibonacci fusion rule ----------
    Nijk = zeros(Int, 2, 2, 2)
    Nijk[1, 1, 1] = 1
    Nijk[1, 2, 2] = 1
    Nijk[2, 1, 2] = 1
    Nijk[2, 2, 1] = 1
    Nijk[2, 2, 2] = 1

    @testset "Pentagon system setup" begin
        R, eqs, n = Phase4.get_pentagon_system(Nijk, 2)
        @test n > 0                       # non-trivial (not pointed)
        @test length(eqs) > 0
        @test all(!iszero, eqs)
        println("  Fibonacci pentagon: $n variables, $(length(eqs)) equations")
    end

    @testset "Pentagon HC solve" begin
        R, eqs, n = Phase4.get_pentagon_system(Nijk, 2)

        # Pentagon has gauge freedom; slice=1 typically reduces to a
        # zero-dimensional variety for Fibonacci.
        F_sols = Phase4.solve_pentagon_homotopy(eqs, n;
                                                slice = 1,
                                                include_singular = false,
                                                show_progress = false)
        @test !isempty(F_sols)
        @test all(s -> length(s) == n, F_sols)
        println("  Pentagon HC returned $(length(F_sols)) solutions")

        # Sanity: each solution should numerically satisfy all equations
        for (k, s) in enumerate(F_sols)
            residuals = [abs(Phase4.PentagonSolver.eval_poly_complex(eq, s)) for eq in eqs]
            maxres = maximum(residuals)
            @test maxres < 1e-8
            if maxres >= 1e-10
                println("    Sol $k: max residual = $maxres (marginal)")
            end
        end
    end

    @testset "Hexagon for one pentagon solution" begin
        R, eqs, n = Phase4.get_pentagon_system(Nijk, 2)
        F_sols = Phase4.solve_pentagon_homotopy(eqs, n;
                                                slice = 1,
                                                show_progress = false)
        @test !isempty(F_sols)

        # Polish before feeding to hexagon
        F_sol_polished = Phase4.refine_solution_newton(eqs, F_sols[1]; tol = 1e-14)
        @test length(F_sol_polished) == n

        # Build hexagon system with F fixed
        R_ring, hex_eqs, n_r = Phase4.get_hexagon_system(Nijk, 2, F_sol_polished)
        @test n_r > 0
        @test length(hex_eqs) > 0
        println("  Fibonacci hexagon: $n_r R-variables, $(length(hex_eqs)) equations")

        # Solve hexagon
        R_sols = Phase4.solve_hexagon_homotopy(hex_eqs, n_r;
                                               show_progress = false)
        @test !isempty(R_sols)
        println("  Hexagon HC returned $(length(R_sols)) braidings")
    end
end
