using Test
using ACMG
include("fixtures/frdata.jl")

@testset "finite-field Ising braid representation" begin
    fr = test_ising_fr_data_mod_p_17()
    @test verify_FRData(fr)

    σ1, σ2 = braid_generators_B3(fr, [2, 2, 2]; total_charge = 2)
    @test eltype(σ1) == FpElem
    @test size(σ1) == (2, 2)
    @test σ1 * σ2 * σ1 == σ2 * σ1 * σ2
end
