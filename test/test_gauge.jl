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
        @test isempty(action.F)
        @test isempty(action.R)

        fixed = gauge_fix(system; strategy = :safe)
        @test validate_fr_system(fixed)
        @test fixed.metadata[:gauge_fix_strategy] == :safe
        @test haskey(fixed.metadata, :fixed_symbols)
        @test haskey(fixed.metadata, :residual_gauge_variables)
    end

    @test_throws ErrorException gauge_fix(fr_equation_system(fibonacci_fusion_rules());
                                          strategy = :canonical)
end

@testset "FRData-backed gauge accessors and transforms" begin
    fib = fibonacci_fr_data_mod_p(101)

    @test simples(fib) == [1, 2]
    @test fusion_coeff(fib, 2, 2, 1) == 1
    @test fusion_coeff(fib, "2", "2", "2") == 1
    @test fusion_channels(fib, 2, 2) == [1, 2]
    @test hom_basis(fib, 2, 2, 2) == [1]
    @test gauge_basis_indices(fib, 2, 2, 2) == [1]
    @test validate_frdata_for_gauge(fib)
    @test has_F_symbol(fib, 2, 2, 2, 2; e = 1, f = 2)
    @test has_R_symbol(fib, 2, 2, 1)

    identity = gauge_transform(fib, nothing)
    @test identity isa FRData
    @test identity.F_values == fib.F_values
    @test identity.R_values == fib.R_values

    p = fr_metadata(fib)[:p]
    scalars = Dict((ch[1], ch[2], ch[3]) =>
                   (ch[1] == 1 || ch[2] == 1 ? one(fib.F_values[1]) : FpElem(ch[1] + ch[2] + ch[3], p))
                   for ch in gauge_parameters(fib))
    moved = gauge_transform(fib, GaugeTransform(scalars, Int[], false))
    moved_legacy = gauge_transform(fib.F_values, fib.R_values,
                                   GaugeTransform(scalars, Int[], false);
                                   Nijk = fib.rules.N)
    @test moved.F_values == moved_legacy.F
    @test moved.R_values == moved_legacy.R

    params = GaugeParameters(Dict((a, b, c, 1) => v for ((a, b, c), v) in scalars))
    moved_params = gauge_transform(fib, params)
    @test moved_params.F_values == moved.F_values
    @test moved_params.R_values == moved.R_values

    fixed = canonical_gauge(moved)
    @test fixed isa GaugeFixingResult
    @test fixed.complete
    @test fixed.F == canonical_gauge(moved.F_values, moved.R_values, moved.rules.N).F
    @test gauge_equivalent(fib, moved)
end
