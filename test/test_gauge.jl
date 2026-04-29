using Test
using ACMG

@testset "v0.8 gauge infrastructure" begin
    for rules in (semion_fusion_rules(), toric_code_fusion_rules(),
                  fibonacci_fusion_rules(), ising_fusion_rules())
        gvars = gauge_variables(rules)
        @test length(gvars) == length(ACMG.admissible_triples(rules))
        @test length(unique(v.variable.name for v in gvars)) == length(gvars)

        system = fr_equation_system(rules)
        action = gauge_transform(system)
        @test length(action.gauge_variables) == length(gvars)
        @test length(action.F) == length(system.fvars)
        @test length(action.R) == length(system.rvars)

        fixed = gauge_fix(system; strategy = :safe)
        @test validate_fr_system(fixed)
        @test fixed.metadata[:gauge_fix_strategy] == :safe
        @test haskey(fixed.metadata, :fixed_symbols)
        @test haskey(fixed.metadata, :residual_gauge_variables)
    end

    @test_throws ErrorException gauge_fix(fr_equation_system(fibonacci_fusion_rules());
                                          strategy = :canonical)
end
