using Test
using Oscar
using ACMG

@testset "Exact ModularDataLift" begin
    @testset "DiscreteLogTable basic" begin
        tbl = DiscreteLogTable(5, 41, 10)
        @test tbl.lookup[1] == 0
        @test tbl.lookup[10] == 1
        @test tbl.lookup[mod(10^2, 41)] == 2
        @test tbl.lookup[mod(10^3, 41)] == 3
        @test tbl.lookup[mod(10^4, 41)] == 4
        @test length(tbl.lookup) == 5
        @test !haskey(tbl.lookup, 2)
    end

    @testset "T lift lands in Q(zeta_N)" begin
        N, p = 20, 41
        zeta_Fp = ACMG.find_zeta_in_Fp(N, p)
        K, zeta = cyclotomic_field(N)
        ks = [0, 1, 7, 11]
        T_Q = [zeta^k for k in ks]
        T_Fp = [ACMG.cyclotomic_to_Fp(t, zeta_Fp, p) for t in T_Q]

        T_recovered = lift_T_Fp_to_cyclotomic(T_Fp, N, p, zeta_Fp)

        @test all(t -> parent(t) === K, T_recovered)
        @test T_recovered == T_Q
    end

    @testset "S lift lands in Q(zeta_N)" begin
        N, d = 20, 5
        K, zeta = cyclotomic_field(N)
        recon = [(0, 1) (1, 0); (1, 0) (0, -1)]

        S = lift_S_sqrtd_to_cyclotomic(recon, d, N; scale = 2)
        sqrt5 = one(K) + 2 * (zeta^(N ÷ 5) + zeta^(4N ÷ 5))

        @test base_ring(S) === K
        @test S[1, 1] == K(1) // K(2)
        @test S[1, 2] == inv(K(2) * sqrt5)
        @test S[2, 2] == -K(1) // K(2)
    end

    @testset "Non-root-of-unity T_Fp value errors" begin
        N, p = 5, 41
        @test_throws Exception lift_T_Fp_to_cyclotomic([2], N, p, 10)
    end
end
