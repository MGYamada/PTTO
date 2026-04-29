using Test
using ACMG

@testset "v0.8 R-symbol variables" begin
    for rules in (semion_fusion_rules(), toric_code_fusion_rules(),
                  fibonacci_fusion_rules(), ising_fusion_rules())
        rvars = rsymbol_variables(rules)
        @test !isempty(rvars)
        @test length(unique(v.variable.name for v in rvars)) == length(rvars)
        @test length(rvars) == length(ACMG.admissible_triples(rules))
        for v in rvars
            @test is_admissible(rules, v.a, v.b, v.c)
        end
    end
end
