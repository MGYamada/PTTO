using Test
using ACMG
using Oscar

function _semion_fusion()
    Nijk = zeros(Int, 2, 2, 2)
    Nijk[1, 1, 1] = 1
    Nijk[1, 2, 2] = 1
    Nijk[2, 1, 2] = 1
    Nijk[2, 2, 1] = 1
    return Nijk
end

@testset "Phase 4 solution limits are configurable" begin
    R, x = polynomial_ring(QQ, 1, :x)
    eqs = [x[1]^2 - 1]

    limited = solve_pentagon_modular_crt(eqs, 1;
                                         conductor = 8,
                                         primes = Int[],
                                         max_solutions = 1)
    expanded = solve_pentagon_modular_crt(eqs, 1;
                                          conductor = 8,
                                          primes = Int[],
                                          max_solutions = 4)
    crt_lifted = solve_pentagon_modular_crt(eqs, 1;
                                            conductor = 1,
                                            primes = [5, 7],
                                            max_solutions = 4,
                                            reconstruction_bound = 1,
                                            denominator_bound = 1,
                                            exact_fallback = false)
    gb_data = ACMG._compute_modular_groebner_data(eqs, 1, [17])
    fp_points = ACMG._enumerate_modular_triangular_solutions(first(gb_data);
                                                             max_points = 4)

    @test length(limited) == 1
    @test length(expanded) == 2
    K_crt = parent(crt_lifted[1][1])
    @test Set([s[1] for s in crt_lifted]) == Set([K_crt(-1), K_crt(1)])
    @test fp_points.complete
    @test sort(fp_points.points) == [[1], [16]]
    @test_throws ErrorException solve_pentagon_modular_crt(eqs, 1;
                                                           conductor = 8,
                                                           primes = Int[],
                                                           max_solutions = 0)

    half_eqs = [2 * x[1] - 1]
    half = solve_pentagon_modular_crt(half_eqs, 1;
                                      conductor = 1,
                                      primes = [5, 7],
                                      reconstruction_bound = 1,
                                      denominator_bound = 2,
                                      exact_fallback = false)
    @test length(half) == 1
    @test half[1][1] == parent(half[1][1])(1) // parent(half[1][1])(2)
end

function _fibonacci_fusion()
    Nijk = _semion_fusion()
    Nijk[2, 2, 2] = 1
    return Nijk
end

function _ising_fusion()
    Nijk = zeros(Int, 3, 3, 3)
    for a in 1:3
        Nijk[1, a, a] = 1
        Nijk[a, 1, a] = 1
    end
    Nijk[2, 2, 1] = 1
    Nijk[2, 2, 3] = 1
    Nijk[2, 3, 2] = 1
    Nijk[3, 2, 2] = 1
    Nijk[3, 3, 1] = 1
    return Nijk
end

function _twists(data)
    r = length(data.labels)
    return [data.T[i, i] for i in 1:r]
end

function _eval_exact_equation(eq, values)
    K = parent(values[1])
    v = zero(K)
    for (c, m) in zip(coefficients(eq), monomials(eq))
        term = K(c)
        for (i, d) in enumerate(degrees(m))
            d > 0 && (term *= values[i]^d)
        end
        v += term
    end
    return v
end

@testset "Exact cyclotomic Phase 4 roundtrip" begin
    cases = [
        (name = "semion", N = 8, primes = [17, 41],
         fusion = _semion_fusion(), data = semion_modular_data()),
        (name = "Fibonacci", N = 20, primes = [41, 61],
         fusion = _fibonacci_fusion(), data = fibonacci_modular_data()),
        (name = "Ising", N = 16, primes = [17, 97, 113],
         fusion = _ising_fusion(), data = ising_modular_data()),
    ]

    for case in cases
        @testset "$(case.name)" begin
            result = compute_FR_from_ST(case.fusion;
                                        conductor = case.N,
                                        primes = case.primes,
                                        S = case.data.S,
                                        T = _twists(case.data))
            roundtrip = ACMG._modular_data_roundtrip(result.F, result.R,
                                                     case.fusion,
                                                     case.data.S,
                                                     _twists(case.data),
                                                     case.N)

            @test !isempty(result.F)
            @test !isempty(result.R)
            @test !hasproperty(result, :candidates)
            @test is_gauge_fixed(result.F, case.fusion)
            @test canonical_gauge(result.F, result.R, case.fusion).gauge.complete
            @test roundtrip.ok
            @test roundtrip isa ACMG.FRRoundtripReport
            @test iszero(roundtrip.S_max)
            @test iszero(roundtrip.T_max)
            @test iszero(roundtrip.S_error)
            @test iszero(roundtrip.T_error)

            K = parent(result.F[1])
            symmetric_channel_value(ch) =
                case.name == "Ising" ? one(K) :
                (ch[1] == 1 || ch[2] == 1) ? one(K) :
                K(10 * min(ch[1], ch[2]) + max(ch[1], ch[2]) + ch[3])
            gauge = GaugeTransform(Dict(ch => symmetric_channel_value(ch)
                                        for ch in gauge_parameters(case.fusion)),
                                   Int[], true)
            moved = gauge_transform(result.F, result.R, gauge; Nijk = case.fusion)
            _, pentagon_eqs, _ = get_pentagon_system(case.fusion, size(case.fusion, 1))
            @test all(iszero(_eval_exact_equation(eq, moved.F)) for eq in pentagon_eqs)
            if case.name != "Ising"
                _, hexagon_eqs, _ = get_hexagon_system(case.fusion, size(case.fusion, 1), moved.F;
                                                       context = case.data.context)
                @test all(iszero(_eval_exact_equation(eq, moved.R)) for eq in hexagon_eqs)
            end
        end
    end
end

@testset "Phase 4 roundtrip is checked against the requested target" begin
    data = semion_modular_data()
    Nijk = _semion_fusion()
    twists = _twists(data)
    result = compute_FR_from_ST(Nijk;
                                conductor = 8,
                                primes = [17, 41],
                                S = data.S,
                                T = twists,
                                return_all = true)
    @test hasproperty(result, :candidates)

    K = parent(twists[1])
    bad_twists = [one(K), one(K)]
    selection = ACMG._select_fr_for_st(result.candidates, Nijk,
                                       data.S, bad_twists, 8)
    target_roundtrip = ACMG._modular_data_roundtrip(selection.selected.F,
                                                    selection.selected.R,
                                                    Nijk, data.S,
                                                    bad_twists, 8)
    self_roundtrip = ACMG._modular_data_roundtrip(selection.selected.F,
                                                  selection.selected.R,
                                                  Nijk,
                                                  selection.score.S_roundtrip,
                                                  selection.score.T_roundtrip,
                                                  8)

    @test !selection.selected_ok
    @test !target_roundtrip.ok
    @test !ACMG._fr_roundtrip_attachable(target_roundtrip)
    @test self_roundtrip.ok
end

@testset "Reference Ising FRData roundtrips target T" begin
    data = ising_modular_data()
    fr = ising_fr_data()
    twists = _twists(data)
    roundtrip = ACMG._modular_data_roundtrip(fr.F_values, fr.R_values,
                                             fr.rules.N, data.S, twists, 16)

    @test roundtrip.ok
    @test iszero(roundtrip.S_error)
    @test iszero(roundtrip.T_error)
    @test roundtrip.T_roundtrip == twists
end
