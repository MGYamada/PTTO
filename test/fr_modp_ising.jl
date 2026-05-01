using Test
using ACMG

@testset "finite-field Phase-4 FR solve: Ising" begin
    fr = ising_fr_data_mod_p(17)

    @test fr_metadata(fr)[:solver_status] == :solved
    @test verify_FRData(fr)
    @test simples(fr) == [1, 2, 3]
    @test fusion_rule(fr) == ising_fusion_rules()
end
