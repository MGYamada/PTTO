using Test
using ACMG

@testset "Gauge stabilizer infrastructure" begin
    @testset "toy multiplicative action" begin
        p = 5
        group = [ACMG.FpElem(i, p) for i in 1:(p - 1)]
        action = (x, g) -> g * x

        nonzero = stabilizer(ACMG.FpElem(2, p), group; action = action)
        @test stabilizer_order(nonzero) == 1
        @test is_trivial_stabilizer(nonzero)
        @test automorphisms(nonzero) == [ACMG.FpElem(1, p)]
        @test stabilizer_metadata(nonzero).method == :bruteforce

        zero = stabilizer(ACMG.FpElem(0, p), group; action = action)
        @test stabilizer_order(zero) == p - 1
        @test !is_trivial_stabilizer(zero)
        @test ACMG.FpElem(1, p) in automorphisms(zero)

        order_only = stabilizer(ACMG.FpElem(0, p), group;
                                action = action,
                                return_automorphisms = false)
        @test stabilizer_order(order_only) == p - 1
        @test automorphisms(order_only) === nothing
    end

    @testset "problem/result containers" begin
        problem = StabilizerProblem(3, [1, -1]; metadata = (label = :sign_action,))
        @test problem.solution == 3
        @test problem.gauge_group == [1, -1]
        @test problem.metadata.label == :sign_action

        result = StabilizerResult([:id], 1; metadata = (method = :unit_test,))
        @test automorphisms(result) == [:id]
        @test stabilizer_order(result) == 1
        @test is_trivial_stabilizer(result)
        @test stabilizer_metadata(result).method == :unit_test

        @test_throws ErrorException stabilizer(StabilizerProblem(1, nothing))
    end

    @testset "finite-field FRData integration" begin
        fr = semion_fr_data_mod_p(17)
        group = ACMG.finite_field_gauge_group(fr)
        problem = StabilizerProblem(fr, group)

        equations = stabilizer_equations(problem)
        @test equations isa StabilizerEquations
        @test equations.metadata.kind == :toric_character_equations
        @test equations.metadata.field == :F_17
        @test equations.metadata.gauge_convention == :full_channel_scalar
        @test equations.metadata.includes_unit_channels
        @test equations.metadata.includes_ineffective_kernel
        @test length(equations.equations) == length(stabilizer_equations(problem).equations)

        id = identity_gauge(fr)
        fixed = ACMG.apply_gauge(fr, id)
        @test ACMG.F_values(fixed) == ACMG.F_values(fr)
        @test ACMG.R_values(fixed) == ACMG.R_values(fr)

        tiny = StabilizerProblem(fr, [id])
        result = stabilizer(tiny)
        @test stabilizer_order(result) == 1
        @test automorphisms(result)[1] == id

        @test_throws ErrorException stabilizer(problem; max_enumeration = 10)
    end

    @testset "finite-field stabilizer agrees with SNF count for tiny primes" begin
        fr = semion_fr_data_mod_p(5)
        group = ACMG.finite_field_gauge_group(fr)
        @test group.metadata.gauge_convention == :full_channel_scalar
        @test group.metadata.gauge_group_kind == :full_channel_toric_gauge
        @test group.metadata.includes_unit_channels
        @test group.metadata.includes_ineffective_kernel

        brute = stabilizer(fr, group; max_enumeration = 10_000,
                           return_automorphisms = false)
        symbol_data = (
            coordinates = ACMG.fr_symbol_coordinates(fr),
            values = vcat(F_values(fr), R_values(fr)),
            parameters = gauge_parameters(fr),
        )
        snf_order = ACMG.stabilizer_size_mod_p(symbol_data, fusion_rule(fr), 5)
        @test stabilizer_order(brute) == snf_order
    end
end
