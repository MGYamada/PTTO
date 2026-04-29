using Test
using ACMG

@testset "v0.8 finite-field reduction" begin
    for rules in (semion_fusion_rules(), fibonacci_fusion_rules(), ising_fusion_rules())
        system = gauge_fix(fr_equation_system(rules))
        fp = reduce_mod_p(system, 11)
        @test fp.p == 11
        @test length(fp.variables) == length(system.variables)
        @test length(fp.equations) == length(system.equations)
        @test fp.metadata[:base_field] == :F_11
    end

    @test_throws ErrorException reduce_mod_p(fr_equation_system(semion_fusion_rules()), 12)

    meta = frobenius_metadata(11, 5)
    @test meta[:p_mod_conductor] == 1
    @test meta[:split_prime_hint]
    @test check_modular_data(:candidate, :candidate)
    @test_throws ErrorException cyclotomic_reconstruct(Dict(); conductor = 5)
end
