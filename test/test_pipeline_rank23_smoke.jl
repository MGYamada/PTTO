using Test
using Oscar
using ACMG

@testset "Rank 2-3 pipeline smoke tests" begin
    cases = [
        (name = "semion", N = 8, primes = [17, 41], max_rank = 2, expected_rank = 2,
         skip_FR = false),
        (name = "Fibonacci", N = 20, primes = [41, 61], max_rank = 2, expected_rank = 2,
         skip_FR = false),
        (name = "Ising", N = 16, primes = [17, 97], max_rank = 3, expected_rank = 3,
         skip_FR = false),
    ]

    for case in cases
        @testset "$(case.name): N=$(case.N)" begin
            classified = ACMG.classify_mtcs_at_conductor(case.N;
                                                         max_rank = case.max_rank,
                                                         primes = case.primes,
                                                         skip_FR = case.skip_FR,
                                                         verbose = false)

            @test classified isa Vector{ACMG.ClassifiedMTC}
            @test !isempty(classified)
            @test any(m -> m.rank == case.expected_rank, classified)
            @test all(m -> m.N == case.N, classified)
            @test all(m -> parent(m.T_cyclotomic[1]) == cyclotomic_field(case.N)[1],
                      classified)
            if case.name == "Ising"
                rank3 = filter(m -> m.rank == 3, classified)
                @test !isempty(rank3)
                @test all(m -> m.verify_report !== nothing && m.verify_report.ok, rank3)
                @test all(m -> iszero(m.verify_report.T_error), rank3)
                @test all(m -> m.T_cyclotomic[3]^8 == -one(parent(m.T_cyclotomic[3])), rank3)
            end
        end
    end

    @testset "toric gauge fixing can be disabled" begin
        classified = ACMG.classify_mtcs_at_conductor(8;
                                                     max_rank = 2,
                                                     primes = [17, 41],
                                                     skip_FR = false,
                                                     toric_gauge_fixing = false,
                                                     verbose = false)
        @test !isempty(classified)
        @test any(m -> m.rank == 2 && m.verify_report !== nothing && m.verify_report.ok,
                  classified)
    end
end
