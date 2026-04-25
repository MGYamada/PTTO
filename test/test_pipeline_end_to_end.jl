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
end
