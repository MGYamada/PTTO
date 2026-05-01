using Test
using ACMG
include("fixtures/frdata.jl")

@testset "finite-field Ising FRData fixture" begin
    fr = test_ising_fr_data_mod_p_17()

    @test fr_metadata(fr)[:solver_status] == :fixture
    @test verify_FRData(fr)
    @test simples(fr) == [1, 2, 3]
    @test fusion_rule(fr) == ising_fusion_rules()
end
