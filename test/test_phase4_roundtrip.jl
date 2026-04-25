using Test
using ACMG

function _semion_fusion()
    Nijk = zeros(Int, 2, 2, 2)
    Nijk[1, 1, 1] = 1
    Nijk[1, 2, 2] = 1
    Nijk[2, 1, 2] = 1
    Nijk[2, 2, 1] = 1
    return Nijk
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

@testset "Exact cyclotomic Phase 4 roundtrip" begin
    cases = [
        (name = "semion", N = 8, primes = [17, 41],
         fusion = _semion_fusion(), data = semion_modular_data()),
        (name = "Fibonacci", N = 20, primes = [41, 61],
         fusion = _fibonacci_fusion(), data = fibonacci_modular_data()),
        (name = "Ising", N = 16, primes = [17, 97],
         fusion = _ising_fusion(), data = ising_modular_data()),
    ]

    for case in cases
        @testset "$(case.name)" begin
            result = compute_FR_from_ST(case.fusion;
                                        conductor = case.N,
                                        primes = case.primes)
            roundtrip = ACMG._modular_data_roundtrip(result.F, result.R,
                                                     case.fusion,
                                                     case.data.S,
                                                     _twists(case.data),
                                                     case.N)

            @test !isempty(result.F)
            @test !isempty(result.R)
            @test roundtrip.ok
            @test iszero(roundtrip.S_max)
            @test iszero(roundtrip.T_max)
        end
    end
end
