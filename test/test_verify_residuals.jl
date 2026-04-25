using Test
using ACMG

"""
Phase 4 Verify tests.

Important finding: Pentagon HC with slice=1 returns 4 F-solutions that
are NOT all gauge-equivalent. Empirically:
- 2 of them (F[1] ≈ -1/φ, F[4] ≈ +1/φ) yield hexagons that produce
  Fibonacci MTC and its complex conjugate.
- 2 of them (F[1] ≈ +φ, F[4] ≈ -φ) yield hexagons for a DIFFERENT MTC
  with the same fusion ring (likely the Yang-Lee-type variant or a
  different F gauge-class giving inequivalent braided categories).

Strategy: iterate over every (F, R) pair and find the one matching
T_expected via ribbon. This reflects how Phase 4 will be used in
practice: the modular data (S, T) from Phase 3 pins down a specific
MTC, and Phase 4's job is to find a (F, R) realising it.

Fibonacci data:
- rank 2, N=5, τ ⊗ τ = 1 ⊕ τ
- twists: θ_1 = 1, θ_τ = exp(4πi/5)
"""

"""
    find_matching_FR(Nijk, T_expected; atol=1e-8)
        -> (F, R, meta) or nothing

Iterate over all pentagon × hexagon combinations and return the first
(F, R) pair whose ribbon residual against `T_expected` is below `atol`.
`meta` is a NamedTuple with diagnostic info.
"""
function find_matching_FR(Nijk, T_expected; atol = 1e-8)
    r = size(Nijk, 1)
    _R, eqs, n = get_pentagon_system(Nijk, r)
    F_sols = solve_pentagon_homotopy(eqs, n;
                                            slice = 1, show_progress = false)

    best = (ribbon_max = Inf, F = nothing, R = nothing,
            f_idx = 0, r_idx = 0, n_pentagon = length(F_sols),
            n_tried = 0, n_matches = 0)

    for (fi, F_raw) in enumerate(F_sols)
        F = refine_solution_newton(eqs, F_raw; tol = 1e-14)
        local R_sols
        try
            _, hex_eqs, n_r = get_hexagon_system(Nijk, r, F)
            R_sols = solve_hexagon_homotopy(hex_eqs, n_r;
                                                   show_progress = false)
        catch
            continue
        end

        for (ri, R) in enumerate(R_sols)
            best = (; best..., n_tried = best.n_tried + 1)
            rib_max = maximum(ribbon_residuals(R, T_expected, Nijk))
            if rib_max < atol
                best = (; best..., n_matches = best.n_matches + 1)
                if rib_max < best.ribbon_max
                    best = (; best..., ribbon_max = rib_max,
                            F = F, R = R, f_idx = fi, r_idx = ri)
                end
            end
        end
    end
    return best
end

@testset "Phase 4 Verify" begin
    # ----- Fibonacci setup -----
    Nijk = zeros(Int, 2, 2, 2)
    Nijk[1, 1, 1] = 1
    Nijk[1, 2, 2] = 1
    Nijk[2, 1, 2] = 1
    Nijk[2, 2, 1] = 1
    Nijk[2, 2, 2] = 1

    T_expected = ComplexF64[1.0, exp(4π * im / 5)]

    best = find_matching_FR(Nijk, T_expected)

    println("  Pentagon solutions: $(best.n_pentagon)")
    println("  (F, R) pairs tried: $(best.n_tried)")
    println("  Fibonacci matches:  $(best.n_matches)")
    if best.F !== nothing
        println("  Best match at F[$(best.f_idx)] × R[$(best.r_idx)], " *
                "ribbon_max = $(best.ribbon_max)")
    end

    @testset "A Fibonacci realisation exists" begin
        # Whatever random slice HC used, SOME (F, R) pair must realise
        # Fibonacci. Historically pentagon returns 4 F-solutions and at
        # least 2 of them give Fibonacci hexagons (another 2 give a
        # different MTC), so we expect >= 2 matches, but we only *require*
        # >= 1 to keep the test robust.
        @test best.n_matches >= 1
        @test best.F !== nothing
        @test best.R !== nothing
    end

    # From here on use the matching pair
    F = best.F
    R = best.R

    @testset "Pentagon residuals" begin
        pent_res = pentagon_residuals(F, Nijk)
        max_pent = maximum(pent_res)
        println("  Pentagon max residual: $max_pent")
        @test max_pent < 1e-8
    end

    @testset "Hexagon residuals" begin
        hex_res = hexagon_residuals(F, R, Nijk)
        max_hex = maximum(hex_res)
        println("  Hexagon max residual: $max_hex")
        @test max_hex < 1e-8
    end

    @testset "R-block extraction" begin
        R11_1 = extract_R_block(R, Nijk, 1, 1, 1)
        @test size(R11_1) == (1, 1)

        R22_1 = extract_R_block(R, Nijk, 2, 2, 1)
        @test size(R22_1) == (1, 1)
        R22_2 = extract_R_block(R, Nijk, 2, 2, 2)
        @test size(R22_2) == (1, 1)

        # No such block: Nijk[1,1,2] = 0
        empty_block = extract_R_block(R, Nijk, 1, 1, 2)
        @test size(empty_block) == (0, 0)
    end

    @testset "verify_mtc roundtrip (all residuals near zero)" begin
        report = verify_mtc(F, R, Nijk; T = T_expected)
        @test report.rank == 2
        @test report.pentagon_max < 1e-8
        @test report.hexagon_max < 1e-8
        @test maximum(ribbon_residuals(R, T_expected, Nijk)) < 1e-8
        println("  $report")
    end

    @testset "compute_FR_from_ST solves from fusion data only" begin
        fr = compute_FR_from_ST(Nijk;
                                return_all = true,
                                pentagon_slice = 1,
                                show_progress = false,
                                verbose = false)
        @test fr.F !== nothing
        @test fr.R !== nothing
        @test fr.report !== nothing
        @test length(fr.candidates) >= 1
    end
end
