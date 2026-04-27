using Test
using Oscar
using ACMG

@testset "Rank 2-3 pipeline smoke tests" begin
    cases = [
        (name = "semion", N = 8, primes = [17, 41], max_rank = 2, expected_rank = 2),
        (name = "Fibonacci", N = 20, primes = [41, 61], max_rank = 2, expected_rank = 2),
    ]

    for case in cases
        @testset "$(case.name): N=$(case.N)" begin
            classified = ACMG.classify_mtcs_at_conductor(case.N;
                                                         max_rank = case.max_rank,
                                                         primes = case.primes,
                                                         verbose = false)

            @test classified isa Vector{ACMG.ClassifiedMTC}
            @test !isempty(classified)
            @test any(m -> m.rank == case.expected_rank, classified)
            @test all(m -> m.N == case.N, classified)
            @test all(m -> parent(m.T_cyclotomic[1]) == cyclotomic_field(case.N)[1],
                      classified)
        end
    end
end
