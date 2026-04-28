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
