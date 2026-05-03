using Test
using Oscar
using ACMG

@testset "CyclotomicContext exact modular data" begin
    @testset "semion over Q(zeta_8)" begin
        ctx = CyclotomicContext(8)
        data = semion_modular_data(ctx)
        K = ACMG.field(ctx)
        z = ACMG.zeta(ctx)

        @test data.context === ctx
        @test ACMG.field(data) === K
        @test ACMG.conductor(data) == 8
        @test ACMG.cond_S(data) == 8
        @test ACMG.cond_T(data) == 4
        @test ACMG.cond_F(data) === nothing
        @test base_ring(data.S) === K
        @test base_ring(data.T) === K
        @test data.T[2, 2] == z^2
        @test data.S[1, 1] == inv(z + z^(-1))
        @test validate_exact_modular_data(data).valid

        sigma3 = galois_action(data, 3)
        @test sigma3.T[2, 2] == z^6
        @test length(galois_orbit(data)) == 4

        red = reduce_mod_p(data, 17)
        @test size(red.S) == (2, 2)
        @test size(red.T) == (2, 2)
        @test red.p == 17
        @test red.conductor == 8
    end

    @testset "Fibonacci over Q(zeta_20)" begin
        ctx = CyclotomicContext(20)
        data = fibonacci_modular_data(ctx)
        z = ACMG.zeta(ctx)

        @test ACMG.conductor(data) == 20
        @test ACMG.cond_S(data) == 20
        @test ACMG.cond_T(data) == 5
        @test data.labels == [:one, :tau]
        @test data.T[2, 2] == z^8
        @test parent(data.S[1, 1]) === ACMG.field(ctx)
        @test validate_exact_modular_data(data).valid

        sigma3 = galois_action(data, 3)
        @test sigma3.T[2, 2] == z^4

        red = reduce_mod_p(data, 41)
        @test size(red.S) == (2, 2)
        @test size(red.T) == (2, 2)
    end

    @testset "Ising over Q(zeta_16)" begin
        ctx = CyclotomicContext(16)
        data = ising_modular_data(ctx)
        z = ACMG.zeta(ctx)

        @test ACMG.conductor(data) == 16
        @test ACMG.cond_S(data) == 16
        @test ACMG.cond_T(data) == 16
        @test data.labels == [:one, :sigma, :psi]
        @test data.T[2, 2] == z
        @test data.T[3, 3] == -one(ACMG.field(ctx))
        @test data.S[2, 2] == zero(ACMG.field(ctx))
        @test validate_exact_modular_data(data).valid

        sigma3 = galois_action(data, 3)
        @test sigma3.T[2, 2] == z^3

        red = reduce_mod_p(data, 17)
        @test size(red.S) == (3, 3)
        @test size(red.T) == (3, 3)
    end

    @testset "family constructor checks conductor" begin
        @test modular_data(:semion).context.N == 8
        @test modular_data(:Fibonacci).context.N == 20
        @test modular_data(:ising).context.N == 16
        @test_throws Exception semion_modular_data(CyclotomicContext(16))
        @test_throws Exception fibonacci_modular_data(CyclotomicContext(5))
        @test_throws Exception ising_modular_data(CyclotomicContext(8))
    end
end
