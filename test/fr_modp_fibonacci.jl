using Test
using ACMG

@testset "finite-field Phase-4 FR solve: Fibonacci" begin
    fr = fibonacci_fr_data_mod_p(101)

    @test fr_metadata(fr)[:solver_status] == :solved
    @test verify_FRData(fr)
    @test simples(fr) == [1, 2]
    @test fusion_rule(fr) == fibonacci_fusion_rules()
end
