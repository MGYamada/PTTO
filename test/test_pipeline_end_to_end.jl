using Test
using ACMG

@testset "Pipeline end-to-end" begin
    @testset "classify_mtcs_at_conductor rank-1 smoke (skip_FR)" begin
        classified = ACMG.classify_mtcs_at_conductor(1;
                                                     max_rank = 1,
                                                     primes = [73, 97],
                                                     skip_FR = true,
                                                     verbose = false)
        @test classified isa Vector{ACMG.ClassifiedMTC}
    end

    @testset "N=20 exact layer is single-pass and roundtrips" begin
        classified = ACMG.classify_mtcs_at_conductor(20;
                                                     max_rank = 2,
                                                     primes = [41, 61],
                                                     skip_FR = false,
                                                     verbose = false)
        @test length(classified) == 6
        @test all(m -> m.verify_report !== nothing && m.verify_report.ok, classified)
        @test all(m -> iszero(m.verify_report.S_max), classified)
        @test all(m -> iszero(m.verify_report.T_max), classified)
        @test !occursin("pent=? hex=?", sprint(show, classified))
    end
end
