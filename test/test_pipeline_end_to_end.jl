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

    @testset "classify_mtcs_at_conductor rank-1 auto-prime exact path" begin
        classified = ACMG.classify_mtcs_at_conductor(1;
                                                     max_rank = 1,
                                                     min_primes = 2,
                                                     verbose = false)
        @test length(classified) == 1
        @test classified[1].verify_fresh
        @test classified[1].verify_report !== nothing
        @test classified[1].verify_report.ok
    end

    @testset "N=7 fixed atomic strata lift exactly" begin
        classified = ACMG.classify_mtcs_at_conductor(7;
                                                     max_rank = 3,
                                                     primes = [29, 43],
                                                     skip_FR = true,
                                                     verbose = false)
        @test length(classified) == 3
        @test any(m -> m.rank == 3, classified)
        @test all(m -> m.verify_fresh, classified)
    end

    @testset "N=8 rank-3 groups without F/R do not abort pipeline" begin
        classified = ACMG.classify_mtcs_at_conductor(8;
                                                     max_rank = 3,
                                                     primes = [41, 73],
                                                     verbose = false)
        @test length(classified) == 6
        @test count(m -> m.rank == 3, classified) == 2
        @test all(m -> m.verify_report === nothing, filter(m -> m.rank == 3, classified))
        @test all(m -> m.verify_report !== nothing && m.verify_report.ok,
                  filter(m -> m.rank <= 2, classified))
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
