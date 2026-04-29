using Test
using ACMG

function _api_semion_fusion()
    N = zeros(Int, 2, 2, 2)
    N[1, 1, 1] = 1
    N[1, 2, 2] = 1
    N[2, 1, 2] = 1
    N[2, 2, 1] = 1
    return N
end

function _api_fibonacci_fusion()
    N = zeros(Int, 2, 2, 2)
    N[1, 1, 1] = 1
    N[1, 2, 2] = 1
    N[2, 1, 2] = 1
    N[2, 2, 1] = 1
    N[2, 2, 2] = 1
    return N
end

@testset "Public API additions" begin
    data = semion_modular_data()
    Nijk = _api_semion_fusion()

    @test check_modular_relations(data.S, data.T).valid
    @test check_unitarity(data.S).valid
    @test check_verlinde_integrality(data.S).valid
    @test check_twist_balance(data.S, data.T, Nijk).valid
    @test check_vafa_constraints(data.T, Nijk).valid
    @test check_galois_symmetry(data.S, data.T).valid
    @test validate_exact_modular_data(data).valid
    @test validate_exact_mtc(nothing, nothing, data.S, data.T, Nijk).valid

    @test fusion_automorphisms(Nijk) == [[1, 2]]
    @test is_fusion_automorphism(Nijk, [1, 2])
    @test modular_data_automorphisms(data) == [[1, 2]]
    @test is_modular_data_automorphism(data, [1, 2])

    st3 = galois_action(data.S, data.T, 3; context = data.context)
    @test st3.S == galois_action(data, 3).S
    @test st3.T == galois_action(data, 3).T
    anyon_action = galois_anyon_action(data, 3)
    @test anyon_action.perm == [1, 2]
    @test galois_anyon_orbits(data) == [[1], [2]]

    sanity = conductor_sanity_table(data)
    @test length(sanity) == 1
    @test sanity[1].N == 8
    @test sanity[1].fusion_automorphisms == 1
    @test sanity[1].modular_data_automorphisms == 1
    @test sanity[1].galois_units == 4
    @test sanity[1].exact_ok
    @test occursin("| N | rank |", conductor_sanity_markdown(data))

    result = compute_FR_from_ST(Nijk;
                                conductor = 8,
                                primes = [17, 41],
                                S = data.S,
                                T = [data.T[i, i] for i in 1:2])
    canon = canonical_gauge(result.F, result.R, Nijk)
    @test gauge_equivalent(result.F, result.R, canon.F, canon.R, Nijk)
    @test gauge_transform(result.F, result.R, canon.gauge).F == result.F

    phase4 = estimate_phase4_complexity(Nijk)
    @test phase4.score >= 1
    @test recommend_skip_FR(Nijk).estimate == phase4
    @test length(recommend_primes(8, 2; min_count = 2, window = 200)) >= 4

    mtc = ACMG.ClassifiedMTC(8, 8, 2, ACMG.Stratum(Dict(1 => 2), 2), Nijk,
                             2, [17, 41], Int[], true, nothing,
                             data.S, data.T, result.F, result.R,
                             result.report, 1)
    @test fr_status(mtc) == FRSolved
end

@testset "Gauge fixing for Fibonacci coordinates" begin
    Nijk = _api_fibonacci_fusion()
    F = Rational{Int}[3, 5, 7, 11, 13]
    R = Rational{Int}[2, 3, 5, 7, 11]

    plan = gauge_fixing_plan(F, Nijk)
    @test [entry.var_idx for entry in plan] == [2]

    canon = canonical_gauge(F, R, Nijk)
    @test canon.gauge.complete
    @test is_gauge_fixed(canon.F, Nijk)
    @test canon.F[2] == 1

    moved = gauge_transform(F, R,
                            GaugeTransform(Dict((1, 1, 1) => 1//1,
                                                (1, 2, 2) => 1//1,
                                                (2, 1, 2) => 1//1,
                                                (2, 2, 1) => 3//1,
                                                (2, 2, 2) => 2//1),
                                           Int[], true);
                            Nijk = Nijk)
    @test !is_gauge_fixed(moved.F, Nijk)
    @test gauge_equivalent(F, R, moved.F, moved.R, Nijk)
    @test canonical_gauge(moved.F, moved.R, Nijk).F == canon.F
end
