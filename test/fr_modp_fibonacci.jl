using Test
using ACMG

@testset "finite-field exact F/R solve: Fibonacci" begin
    fr = fibonacci_fr_data_mod_p(101)

    @test fr_metadata(fr)[:solver_status] == :solved
    @test verify_FRData(fr)
    @test simples(fr) == [1, 2]
    @test fusion_rule(fr) == fibonacci_fusion_rules()

    gf = fr_metadata(fr)[:gauge_fixing]
    @test gf[:gauge_fix_complete]
    @test !isempty(gf[:fixed_f_indices])
    @test gf[:f_variables] < gf[:original_f_variables]
    @test gf[:hexagon_variables] < gf[:original_hexagon_variables]
    @test fr_metadata(fr)[:reduced_solution_count] == fr_metadata(fr)[:solution_count]
end
