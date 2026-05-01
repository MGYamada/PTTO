@testset "API stability documentation" begin
    root = dirname(@__DIR__)
    @test isfile(joinpath(root, "docs", "src", "api_stability.md"))

    documented_experimental = [
        :solve_fr_mod_p,
        :lift_higher_central_charge,
        :solve_finite_field,
        :cyclotomic_reconstruct,
        :generated_subgroup,
        :finite_group_diagnostics,
        :generated_matrix_algebra,
        :commutant,
        :zariski_closure_diagnostics,
    ]
    @test all(sym -> isdefined(ACMG, sym), documented_experimental)
    @test :solve_FR_mod_p ∉ names(ACMG)
end

@testset "experimental warning helper" begin
    old = get(ENV, "ACMG_WARN_EXPERIMENTAL", nothing)
    try
        delete!(ENV, "ACMG_WARN_EXPERIMENTAL")
        @test ACMG.warn_experimental("unit-test-unset") == false

        ENV["ACMG_WARN_EXPERIMENTAL"] = "true"
        @test_warn r"unit-test-enabled is an experimental ACMG API" ACMG.warn_experimental("unit-test-enabled")
        @test ACMG.warn_experimental("unit-test-enabled") == false
    finally
        if old === nothing
            delete!(ENV, "ACMG_WARN_EXPERIMENTAL")
        else
            ENV["ACMG_WARN_EXPERIMENTAL"] = old
        end
    end
end
