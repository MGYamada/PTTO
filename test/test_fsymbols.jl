using Test
using ACMG

@testset "v0.8 F-symbol variables" begin
    examples = [
        semion_fusion_rules(),
        toric_code_fusion_rules(),
        fibonacci_fusion_rules(),
        ising_fusion_rules(),
    ]

    for rules in examples
        @test is_multiplicity_free(rules)
        fvars = fsymbol_variables(rules)
        @test !isempty(fvars)
        @test length(unique(v.variable.name for v in fvars)) == length(fvars)
        for v in fvars
            @test is_admissible(rules, v.a, v.b, v.e)
            @test is_admissible(rules, v.e, v.c, v.d)
            @test is_admissible(rules, v.b, v.c, v.f)
            @test is_admissible(rules, v.a, v.f, v.d)
        end
    end

    N = zeros(Int, 2, 2, 2)
    N[1, 1, 1] = 1
    N[1, 2, 2] = 1
    N[2, 1, 2] = 1
    N[2, 2, 1] = 2
    @test !is_multiplicity_free(N)
    @test_throws ErrorException require_multiplicity_free(N)
end
