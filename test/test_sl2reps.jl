using Test
using ACMG

"""
Tests for SL2Reps module.

NOTE: These tests require Oscar.jl with the GAP SL2Reps package installed.
Install the SL2Reps GAP package once via:

    julia> using Oscar
    julia> Oscar.GAP.Packages.install("SL2Reps")

Tests will skip gracefully if Oscar is not available.

Run with:  julia --project=. test/test_sl2reps.jl
"""

@testset "SL2Reps helpers (no dependencies)" begin
    @test all_divisors(1) == [1]
    @test all_divisors(6) == [1, 2, 3, 6]
    @test all_divisors(12) == [1, 2, 3, 4, 6, 12]
    @test all_divisors(24) == [1, 2, 3, 4, 6, 8, 12, 24]
end

# The following tests require Oscar.jl (which bundles GAP + SL2Reps package).
# Wrapped in try/catch so the test suite passes in environments without Oscar.

@testset "Atomic catalog construction (requires Oscar + SL2Reps)" begin
    oscar_available = try
        @eval using Oscar
        true
    catch e
        @info "Oscar.jl not available, skipping catalog tests" exception=e
        false
    end

    if oscar_available
        # Small conductor test: N = 5 (Fibonacci / Yang-Lee related)
        @testset "N=5 catalog" begin
            catalog = build_atomic_catalog(5; max_rank = 10, verbose = false)
            @test length(catalog) >= 1
            # SL(2, Z/5) has well-known irreps; at least the 2d one should appear
            dims = [a.dim for a in catalog]
            @test 2 in dims  # Fibonacci irrep is 2-dimensional
        end

        # N = 8 (semion / Ising related)
        @testset "N=8 catalog" begin
            catalog = build_atomic_catalog(8; max_rank = 10, verbose = false)
            @test length(catalog) >= 1
            dims = [a.dim for a in catalog]
            @test 2 in dims  # semion
        end

        # N = 16 (Ising 3d irrep)
        @testset "N=16 catalog" begin
            catalog = build_atomic_catalog(16; max_rank = 10, verbose = false)
            dims = [a.dim for a in catalog]
            @test 3 in dims  # Ising 3-dimensional irrep
        end
    end
end
