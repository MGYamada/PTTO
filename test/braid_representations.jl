using Test
using ACMG

@testset "braid representation fusion bases" begin
    fib = fibonacci_fusion_rules()
    @test fusion_paths(fib, [2, 2, 2], 2) == fusion_paths(fib, [2, 2, 2], 2)
    @test ACMG.dim(fusion_basis(fib, [2, 2, 2], 2)) == 2
    @test ACMG.dim(fusion_basis(fib, [2, 2, 2, 2], 1)) == 2
    @test isempty(fusion_paths(fib, [1, 1], 2))

    semion = semion_fusion_rules()
    @test ACMG.dim(fusion_basis(semion, [2, 2, 2], 2)) == 1
end

@testset "finite-field FRData solves TensorCategories equations" begin
    cases = [
        semion_fr_data_mod_p(17),
        fibonacci_fr_data_mod_p(101),
    ]
    for fr_data in cases
        rules = fr_data.rules
        @test fusion_rule(fr_data) == rules
        @test fr_data isa FRData
        @test isconcretetype(eltype(fr_data.F_values))
        @test eltype(fr_data.R_values) === eltype(fr_data.F_values)
        @test eltype(fr_data.R_inverse_values) === eltype(fr_data.F_values)
        @test !hasproperty(fr_data, :F)
        @test !hasproperty(fr_data, :R)
        @test length(fr_data.F_values) == get_pentagon_system(rules.N, rules.rank)[3]
        @test fr_pentagon_values(fr_data) === F_values(fr_data)
        @test length(fr_hexagon_values(fr_data)) == get_hexagon_fr_system(rules.N, rules.rank)[3]
        @test verify_FRData(fr_data)
    end
end

@testset "braid generator construction" begin
    br = braid_representation(semion_fr_data_mod_p(17), [2, 2, 2], 2)
    @test length(br.generators) == 2
    @test size(br.generators[1]) == (1, 1)

    fib = braid_representation(fibonacci_fr_data_mod_p(101), [2, 2, 2], 2)
    @test fib isa BraidRepresentation
    @test eltype(fib.generators) <: Matrix
    @test isconcretetype(eltype(fib.generators[1]))
    @test length(fib.generators) == 2
    @test size(fib.generators[1]) == (2, 2)

    fib_data = fibonacci_fr_data_mod_p(101)
    fib_tuple = (rules = fib_data.rules,
                 F_values = fib_data.F_values,
                 R_values = fib_data.R_values,
                 R_inverse_values = fib_data.R_inverse_values)
    fib_from_tuple = braid_representation(fib_tuple, [2, 2, 2], 2)
    @test check_braid_relations(fib_from_tuple).ok

    bad = zeros(Int, 2, 2, 2)
    bad[1, 1, 1] = 1
    bad[1, 2, 2] = 1
    bad[2, 1, 2] = 1
    bad[2, 2, 1] = 2
    @test_throws ErrorException fusion_basis(FusionRule(bad), [2, 2], 1)
end
