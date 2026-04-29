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

        toric = toric_code_modular_data()
        @test central_charge(toric) == one(ACMG.field(toric))
        @test higher_central_charge(toric, 1; normalization = :D).value ==
              one(ACMG.field(toric))
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

    @testset "D-normalized periodicity" begin
        for data in (semion_modular_data(), toric_code_modular_data(),
                     fibonacci_modular_data(), ising_modular_data())
            N = ACMG.conductor(data)
            for n in -2:5
                lhs = higher_central_charge(data, n + N; normalization = :D)
                rhs = higher_central_charge(data, n; normalization = :D)
                @test lhs.ok
                @test rhs.ok
                @test lhs.value == rhs.value
            end
        end
    end

    @testset "finite-field Fibonacci prototype agrees after reduction" begin
        data = fibonacci_modular_data()
        sol = solve_FR_mod_p(:fibonacci, 41)

        @test sol.category == :fibonacci
        @test sol.p == 41
        @test sol.conductor == 20
        @test !isempty(sol.F)
        @test !isempty(sol.R)

        for n in 1:5
            exact = normalized_gauss_sum(data; n = n, normalization = :D)
            exact_mod_p = reduce_mod_p(data.context, exact, sol.p)
            finite_field = higher_central_charge(sol, n)
            @test finite_field.ok
            @test finite_field.value == exact_mod_p
        end
    end

    @testset "finite-field method dispatch" begin
        direct = solve_FR_mod_p(:semion, 17; compute_fr = false)
        via_method = higher_central_charge(:semion, 1; method = :finite_field, p = 17)
        @test via_method.value == higher_central_charge(direct, 1).value
    end
end
