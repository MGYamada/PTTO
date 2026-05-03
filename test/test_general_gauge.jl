using Test
using ACMG
import ACMG: FpElem

function _gg_semion_tensor()
    N = zeros(Int, 2, 2, 2)
    N[1, 1, 1] = 1
    N[1, 2, 2] = 1
    N[2, 1, 2] = 1
    N[2, 2, 1] = 1
    return N
end

function _gg_multiplicity_two_tensor()
    N = _gg_semion_tensor()
    N[2, 2, 2] = 2
    return N
end

function _gg_matrix_dict(gauge, value)
    return Dict(f.index => fill(value, f.dim, f.dim) for f in gauge.factors)
end

@testset "GeneralGauge foundation" begin
    @testset "toric examples are GL(1) products" begin
        for rules in (semion_fusion_rules(), fibonacci_fusion_rules(), ising_fusion_rules())
            gauge = general_gauge_data(rules)
            @test is_toric(gauge)
            @test is_multiplicity_free(gauge)
            @test all(f -> f.dim == 1, gauge.factors)
            @test gauge_group_dimension(gauge) == length(gauge.factors)
            @test validate_gauge_transformation(identity_gauge_transformation(gauge))
        end
    end

    @testset "multiplicityful factors and dimensions" begin
        N = _gg_multiplicity_two_tensor()
        gauge = general_gauge_data(N)
        idx = FusionSpaceIndex(2, 2, 2)
        factor = only([f for f in gauge.factors if f.index == idx])

        @test !is_toric(gauge)
        @test !is_multiplicity_free(gauge)
        @test factor.dim == 2
        @test gauge_group_dimension(gauge) == sum(n^2 for n in N if n > 0)
        @test validate_gauge_transformation(identity_gauge_transformation(gauge))
        @test gauge_orbit_dimension(nothing, gauge).upper_bound == gauge_group_dimension(gauge)
        @test gauge_stabilizer(nothing, gauge).reason == :not_computed
    end

    @testset "validation rejects malformed transformations" begin
        gauge = general_gauge_data(_gg_multiplicity_two_tensor(); field = :F_5)
        identity = identity_gauge_transformation(gauge)
        idx = FusionSpaceIndex(2, 2, 2)

        missing = copy(identity.matrices)
        delete!(missing, idx)
        @test_throws GeneralGaugeValidationError validate_gauge_transformation(
            GaugeTransformation(gauge, missing))

        extra = copy(identity.matrices)
        extra[FusionSpaceIndex(2, 1, 1)] = [FpElem(1, 5);;]
        @test_throws GeneralGaugeValidationError validate_gauge_transformation(
            GaugeTransformation(gauge, extra))

        wrong_size = copy(identity.matrices)
        wrong_size[idx] = [FpElem(1, 5);;]
        @test_throws GeneralGaugeValidationError validate_gauge_transformation(
            GaugeTransformation(gauge, wrong_size))

        singular = copy(identity.matrices)
        singular[idx] = [FpElem(1, 5) FpElem(2, 5);
                         FpElem(2, 5) FpElem(4, 5)]
        @test !is_invertible_matrix_over_field(singular[idx], :F_5)
        @test_throws GeneralGaugeValidationError validate_gauge_transformation(
            GaugeTransformation(gauge, singular))
    end

    @testset "action hooks are honest about nonabelian work" begin
        gauge = general_gauge_data(_gg_multiplicity_two_tensor(); field = :F_5)
        identity = identity_gauge_transformation(gauge)
        F = [FpElem(3, 5), FpElem(4, 5)]
        R = [FpElem(2, 5)]

        @test apply_gauge_to_F(F, identity) == F
        @test apply_gauge_to_R(R, identity) == R

        nontrivial = copy(identity.matrices)
        nontrivial[FusionSpaceIndex(2, 2, 2)] =
            [FpElem(1, 5) FpElem(1, 5);
             FpElem(0, 5) FpElem(1, 5)]
        @test validate_gauge_transformation(GaugeTransformation(gauge, nontrivial))
        @test_throws GeneralGaugeActionNotImplementedError apply_gauge_to_F(
            F, GaugeTransformation(gauge, nontrivial))
    end

    @testset "toric integration and pipeline metadata" begin
        data = toric_gauge_data(fibonacci_fusion_rules(); field = :F_101)
        @test data.metadata[:general_gauge] isa GeneralGaugeData
        @test is_toric(data.metadata[:general_gauge])

        nonfree = _gg_multiplicity_two_tensor()
        skipped = ACMG._toric_pre_reconstruction_data(nonfree; gauge_fixing = :general)
        @test !skipped.apply
        @test skipped.general isa GeneralGaugeData
        @test skipped.reason == :not_multiplicity_free
        @test_throws ToricGaugeFixingError ACMG._toric_pre_reconstruction_data(
            nonfree; gauge_fixing = :toric)

        p1, p2 = 73, 97
        Nijk = ones(Int, 1, 1, 1)
        group = Dict(
            p1 => ACMG.MTCCandidate(p1, :searched, [1;;], [1], 1, Nijk, [1], 1),
            p2 => ACMG.MTCCandidate(p2, :searched, [1;;], [1], 1, Nijk, [1], 1),
        )
        classified = ACMG.classify_from_group(group, 24, ACMG.Stratum(Dict(1 => 1), 1),
                                              [p1, p2];
                                              gauge_fixing = :general,
                                              skip_FR = true,
                                              verbose = false)
        @test classified.gauge_data isa GeneralGaugeData
        @test is_toric(classified.gauge_data)
    end
end
