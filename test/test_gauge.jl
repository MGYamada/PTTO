using Test
using ACMG
include("fixtures/frdata.jl")

@testset "gauge infrastructure" begin
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

        fixed = gauge_fix(system; strategy = :toric_snf)
        @test validate_fr_system(fixed)
        @test fixed.metadata[:gauge_fix_strategy] == :toric_snf
        @test haskey(fixed.metadata, :fixed_symbols)
        @test haskey(fixed.metadata, :residual_gauge_variables)
        @test fixed.metadata[:original_f_variables] == system.metadata[:f_variables]
        @test fixed.metadata[:f_variables] ==
              fixed.metadata[:original_f_variables] - length(fixed.metadata[:fixed_f_indices])
        @test fixed.metadata[:hexagon_variables] ==
              fixed.metadata[:original_hexagon_variables] - length(fixed.metadata[:fixed_f_indices])
        @test all(idx -> !(idx in fixed.metadata[:free_hexagon_indices]),
                  fixed.metadata[:fixed_f_indices])
        @test length(fixed.metadata[:fixed_symbols]) == length(fixed.metadata[:fixed_f_indices])
    end

    fib_default = gauge_fix(fr_equation_system(fibonacci_fusion_rules()))
    fib_toric = gauge_fix(fr_equation_system(fibonacci_fusion_rules()); strategy = :toric_snf)
    @test fib_toric.metadata[:fixed_f_indices] == fib_default.metadata[:fixed_f_indices]
    @test fib_toric.metadata[:gauge_fix_complete]
    fib_none = gauge_fix(fr_equation_system(fibonacci_fusion_rules()); strategy = :none)
    @test isempty(fib_none.metadata[:fixed_f_indices])
    @test fib_none.metadata[:gauge_fix_method] == :none
    @test_throws ErrorException gauge_fix(fr_equation_system(fibonacci_fusion_rules());
                                          strategy = :canonical)
    @test_throws ErrorException gauge_fix(fr_equation_system(fibonacci_fusion_rules());
                                          strategy = :safe)
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

@testset "FRData GaugeAction API" begin
    examples = [
        semion_fr_data_mod_p(17),
        fibonacci_fr_data_mod_p(101),
    ]

    for fr in examples
        p = fr_metadata(fr)[:p]
        one_value = fr_value_one(fr)
        identity = identity_gauge(fr)

        @test identity isa GaugeAction
        @test validate_gauge_action(fr, identity)
        @test gauge_degrees_of_freedom(fr) isa Vector{GaugeDegreeOfFreedom}
        @test apply_gauge(fr, identity) == fr
        @test validate_frdata(apply_gauge(fr, identity))

        scalars = Dict((ch[1], ch[2], ch[3]) =>
                       (ch[1] == 1 || ch[2] == 1 ? one_value :
                        FpElem(ch[1] + 2ch[2] + 3ch[3], p))
                       for ch in gauge_parameters(fr))
        action = GaugeAction(scalars; field = Symbol("F_$p"))
        @test validate_gauge_action(fr, action)

        moved = apply_gauge(fr, action)
        @test verify_pentagon(moved)
        @test verify_hexagon(moved)
        @test verify_FRData(moved)

        returned = apply_gauge(moved, inverse_gauge(action))
        @test returned == fr

        composed = compose_gauge(action, inverse_gauge(action))
        @test apply_gauge(fr, composed) == fr

        f_only = apply_gauge(fr, action; target = :F)
        r_only = apply_gauge(fr, action; target = :R)
        @test f_only.R_values == fr.R_values
        @test r_only.F_values == fr.F_values

        constraints = build_gauge_constraints(moved; strategy = :default)
        solved = solve_gauge_constraints(moved, constraints)
        @test solved isa GaugeAction
        normal1 = gauge_normal_form(moved; constraints = constraints)
        normal2 = gauge_normal_form(moved; constraints = constraints)
        @test normal1.F == normal2.F
        @test normal1.R == normal2.R
        @test validate_gauge_fixed(normal1.metadata[:frdata])
    end

    ising = test_ising_fr_data_mod_p_17()
    p = fr_metadata(ising)[:p]
    scalars = Dict((ch[1], ch[2], ch[3]) =>
                   (ch[1] == 1 || ch[2] == 1 ? fr_value_one(ising) :
                    FpElem(ch[1] + 2ch[2] + 3ch[3], p))
                   for ch in gauge_parameters(ising))
    moved_ising = apply_gauge(ising, GaugeAction(scalars; field = :F_17))
    @test verify_pentagon(moved_ising)
    @test verify_hexagon(moved_ising)
    @test apply_gauge(moved_ising, inverse_gauge(GaugeAction(scalars; field = :F_17))) == ising
end
