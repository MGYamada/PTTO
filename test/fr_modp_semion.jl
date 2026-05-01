using Test
using ACMG

@testset "finite-field Phase-4 FR solve: semion" begin
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
end
