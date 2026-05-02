@testset "API stability documentation" begin
    root = dirname(@__DIR__)
    @test isfile(joinpath(root, "docs", "src", "api_stability.md"))

    documented_experimental = [
        :solve_fr_mod_p,
        :solve_finite_field,
        :cyclotomic_reconstruct,
        :generated_subgroup,
        :finite_group_diagnostics,
        :generated_matrix_algebra,
        :commutant,
        :zariski_closure_diagnostics,
        :StabilizerProblem,
        :StabilizerEquations,
        :StabilizerResult,
        :stabilizer,
        :stabilizer_equations,
        :stabilizer_order,
        :automorphisms,
        :is_trivial_stabilizer,
        :stabilizer_metadata,
    ]
    @test all(sym -> isdefined(ACMG, sym), documented_experimental)
    @test :solve_FR_mod_p ∉ names(ACMG)
end

@testset "Julia 1.11 public surface" begin
    root = dirname(@__DIR__)
    project = read(joinpath(root, "Project.toml"), String)
    @test occursin("julia = \"1.11\"", project)
    @test occursin("LinearAlgebra = \"1.11\"", project)

    stable_public = [
        :CyclotomicContext,
        :modular_data,
        :validate_exact_modular_data,
        :higher_central_charge,
        :higher_central_charge_result,
        :higher_central_charge_sequence,
        :HCCGeneratingFunction,
        :fr_status,
        :FRStatus,
        :FRSolved,
        :FRSkipped,
        :FRData,
        :F_symbol,
        :braid_representation,
        :classify_mtcs_auto,
    ]
    @test all(sym -> Base.ispublic(ACMG, sym), stable_public)
    @test Base.isexported(ACMG, :higher_central_charge)
    @test !Base.isexported(ACMG, :higher_central_charge_result)
    @test !Base.isexported(ACMG, :HCCGeneratingFunction)
    @test Base.isexported(ACMG, :fr_status)
    @test Base.isexported(ACMG, :solve_fr_mod_p)
    @test Base.isexported(ACMG, :StabilizerProblem)
    @test Base.isexported(ACMG, :stabilizer)

    removed_hcc_prototype_names = [
        :FRSolutionModP,
        :HigherCentralChargeModPResult,
        :lift_higher_central_charge,
    ]
    @test all(sym -> !isdefined(ACMG, sym), removed_hcc_prototype_names)
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
