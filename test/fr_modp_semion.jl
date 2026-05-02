using Test
using ACMG

@testset "finite-field exact F/R solve: semion" begin
    fr = semion_fr_data_mod_p(17)

    @test fr_metadata(fr)[:solver_status] == :solved
    @test verify_FRData(fr)
    @test simples(fr) == [1, 2]
    @test fusion_rule(fr) == semion_fusion_rules()

    σ1, σ2 = braid_generators_B3(fr, [2, 2, 2]; total_charge = 2)
    @test eltype(σ1) == FpElem
    @test size(σ1) == (1, 1)
    @test σ1 == σ2
    @test !iszero(σ1[1, 1])
    @test σ1[1, 1]^2 == FpElem(-1, 17)

    reference = solve_fr_mod_p(semion_fusion_rules(), 17;
                               solver = :reference,
                               primes = [17],
                               max_solutions = 4)
    @test fr_metadata(reference)[:solver] == :reference
    @test verify_FRData(reference)
    @test_throws ErrorException solve_fr_mod_p(semion_fusion_rules(), 17;
                                               solver = :direct,
                                               primes = [17])
end

@testset "finite-field FRData with no F-symbol variables" begin
    trivial = FusionRule(ones(Int, 1, 1, 1))
    fr = frdata_from_modp_solution(trivial, [1, 1], 5)

    @test isempty(F_values(fr))
    @test verify_pentagon(fr)
    @test verify_hexagon(fr)
    @test verify_FRData(fr)
end
