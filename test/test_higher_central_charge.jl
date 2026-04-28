using Test
using ACMG

@testset "Higher central charge" begin
    @testset "ordinary central charges" begin
        semion = semion_modular_data()
        z8 = ACMG.zeta(semion.context)
        @test gauss_sum_plus(semion) / inv(semion.S[1, 1]) == z8
        @test normalized_gauss_sum(semion) == z8
        @test central_charge(semion) == z8
        @test higher_central_charge(semion).ok
        @test higher_central_charge(semion).value == z8

        ising = ising_modular_data()
        z16 = ACMG.zeta(ising.context)
        @test central_charge(ising) == z16
        @test higher_central_charge(ising).gauss_sum == gauss_sum_plus(ising)

        fib = fibonacci_modular_data()
        z20 = ACMG.zeta(fib.context)
        # This is the convention currently encoded by fibonacci_modular_data().
        @test central_charge(fib) == 2 * z20^7 - z20^3 + z20
    end

    @testset "exact dimensions and normalizations" begin
        semion = semion_modular_data()
        K = ACMG.field(semion)
        @test parent(total_quantum_dimension_squared(semion)) === K
        @test total_quantum_dimension_squared(semion) == K(2)
        @test gauss_sum_minus(semion) == 1 + ACMG.zeta(semion.context)^(-2)
        @test normalized_gauss_sum(semion; normalization = :D2) ==
              gauss_sum_plus(semion) / total_quantum_dimension_squared(semion)
    end

    @testset "higher values are cyclotomic and periodic" begin
        ising = ising_modular_data()
        K = ACMG.field(ising)
        h3 = higher_central_charge(ising; n = 3)
        h19 = higher_central_charge(ising; n = 19)
        @test h3.ok
        @test h19.ok
        @test h3.value == h19.value
        @test parent(h3.value) === K

        semion = semion_modular_data()
        @test gauss_sum_plus(semion; n = 2) == gauss_sum_plus(semion; n = 10)
        blocked = higher_central_charge(semion; n = 2)
        @test !blocked.ok
        @test blocked.status == :not_galois
        @test blocked.value === nothing
        @test normalized_gauss_sum(semion; n = 2, normalization = :D) ==
              gauss_sum_plus(semion; n = 2) / inv(semion.S[1, 1])
    end
end
