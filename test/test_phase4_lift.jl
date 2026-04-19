using Test
using Oscar
using ACMG
using ACMG.Phase4

"""
Phase 4 ModularDataLift tests.

Round-trip tests: take known ℂ-valued (S, T), reduce to F_p / ℤ[√d],
lift back, and compare.

Two focused cases:
1. Fibonacci (N=5, simple T lift only — S doesn't live in ℤ[√d])
2. SU(2)_4 (N=24, both T and S lifts)
"""

@testset "Phase 4 ModularDataLift" begin

    @testset "DiscreteLogTable basic" begin
        # F_41, N=5: ζ_5 = 10 (from Fibonacci test)
        tbl = Phase4.DiscreteLogTable(5, 41, 10)
        @test tbl.lookup[1] == 0                      # ζ^0 = 1
        @test tbl.lookup[10] == 1                     # ζ^1 = 10
        @test tbl.lookup[mod(10^2, 41)] == 2          # ζ^2 = 18
        @test tbl.lookup[mod(10^3, 41)] == 3          # ζ^3
        @test tbl.lookup[mod(10^4, 41)] == 4          # ζ^4
        @test length(tbl.lookup) == 5                 # exactly 5 powers

        # A value NOT in ⟨ζ_5⟩ should fail on discrete_log
        # (pick something explicitly outside: 2 is not a 5th power in F_41
        #  since F_41* has order 40 and 5 | 40, so only 8 elements are 5th
        #  powers. Let's check: 2^8 mod 41 = 256 mod 41 = 10, 10 ≠ 1, so
        #  2 is not a 5th power residue.)
        non_5th_power = 2
        @test !haskey(tbl.lookup, non_5th_power)
    end

    @testset "T lift: Fibonacci" begin
        # From test_fibonacci.jl: N=5, p=41, ζ_5=10, T = [1, 18] (18 = ζ_5²)
        N, p = 5, 41
        zeta_Fp = ACMG.find_zeta_in_Fp(N, p)
        @test zeta_Fp == 10

        T_Fp = [1, 18]
        T_complex = Phase4.lift_T_Fp_to_complex(T_Fp, N, p, zeta_Fp)
        @test length(T_complex) == 2
        @test isapprox(T_complex[1], ComplexF64(1.0); atol = 1e-12)
        # Expected: ζ_5^2 = exp(4πi/5)
        expected = exp(2π * im * 2 / 5)
        @test isapprox(T_complex[2], expected; atol = 1e-12)
    end

    @testset "S lift: ℤ[√3] entries" begin
        # Toy example: suppose reconstruct_S_matrix returned 2√3 · S with
        # one entry equal to (a=3, b=2), representing 3 + 2√3.
        # Then S_entry = (3 + 2√3) / (2·√3) = 3/(2√3) + 1 = √3/2 + 1
        d = 3
        recon = [(3, 2);;]  # 1×1 matrix via the recent comprehension syntax
        S = Phase4.lift_S_sqrtd_to_complex(recon, d; scale = 2)
        @test size(S) == (1, 1)
        expected = (3 + 2 * sqrt(3.0)) / (2 * sqrt(3.0))
        @test isapprox(real(S[1, 1]), expected; atol = 1e-12)
        @test isapprox(imag(S[1, 1]), 0.0; atol = 1e-12)

        # 2×2: identity-like
        recon2 = [(0, 1) (2, 0); (2, 0) (0, 1)]
        S2 = Phase4.lift_S_sqrtd_to_complex(recon2, d; scale = 2)
        # (0 + 1·√3) / (2√3) = 1/2
        @test isapprox(real(S2[1, 1]), 0.5; atol = 1e-12)
        # (2 + 0·√3) / (2√3) = 1/√3
        @test isapprox(real(S2[1, 2]), 1.0 / sqrt(3.0); atol = 1e-12)
    end

    @testset "T lift roundtrip via v0.2 cyclotomic_to_Fp" begin
        # Create a cyclotomic element, reduce to F_p, lift back, compare.
        N, p = 24, 73
        zeta_Fp = ACMG.find_zeta_in_Fp(N, p)

        K, zeta = cyclotomic_field(N)

        # Test values: ζ_24^k for k = 0, 1, 3, 7, 11
        ks = [0, 1, 3, 7, 11]
        T_Q = [zeta^k for k in ks]
        T_Fp = [ACMG.cyclotomic_to_Fp(t, zeta_Fp, p) for t in T_Q]
        T_recovered = Phase4.lift_T_Fp_to_complex(T_Fp, N, p, zeta_Fp)

        for (i, k) in enumerate(ks)
            expected = exp(2π * im * k / N)
            @test isapprox(T_recovered[i], expected; atol = 1e-10)
        end
    end

    @testset "Non-root-of-unity T_Fp value errors" begin
        # If T_Fp contains a value NOT in ⟨ζ_N⟩, we should error loudly.
        N, p = 5, 41
        zeta_Fp = 10
        # 2 is not in ⟨ζ_5⟩ = {1, 10, 18, 16, 37}
        @test_throws Exception Phase4.lift_T_Fp_to_complex([2], N, p, zeta_Fp)
    end
end
