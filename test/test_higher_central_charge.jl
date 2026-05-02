using Test
using ACMG
include("test_imports.jl")

@testset "Higher central charge" begin
    @testset "D2-normalized exact moments" begin
        semion = semion_modular_data()
        K = ACMG.field(semion)
        z8 = ACMG.zeta(semion.context)
        im = z8^2

        @test central_charge(semion) == z8
        @test ACMG.higher_central_charge_result(semion).ok
        @test ACMG.higher_central_charge_result(semion).value == z8
        @test higher_central_charge(semion, 0) == one(K)
        @test higher_central_charge_period(semion) == 4
        @test higher_central_charge_sequence(semion) ==
              [(one(K) + im^n) / K(2) for n in 0:3]
        @test higher_central_charge(semion, -1) == higher_central_charge(semion, 3)
        @test_throws ErrorException higher_central_charge(semion, 1; normalization = :D)

        fib = fibonacci_modular_data()
        Kfib = ACMG.field(fib)
        z20 = ACMG.zeta(fib.context)
        phi = quantum_dimensions(fib)[2]
        theta = z20^8
        D2fib = one(Kfib) + phi^2
        @test central_charge(fib) == z20^7
        @test higher_central_charge(fib, 0) == one(Kfib)
        @test higher_central_charge_period(fib) == 5
        @test higher_central_charge_sequence(fib) ==
              [(one(Kfib) + phi^2 * theta^n) / D2fib for n in 0:4]
        @test higher_central_charge(fib, 6) == higher_central_charge(fib, 1)

        ising = ising_modular_data()
        Kising = ACMG.field(ising)
        z16 = ACMG.zeta(ising.context)
        @test central_charge(ising) == z16
        @test higher_central_charge(ising, 0) == one(Kising)
        @test higher_central_charge_period(ising) == 16
        @test higher_central_charge_sequence(ising) ==
              [(one(Kising) + Kising(2) * z16^n + (-one(Kising))^n) / Kising(4)
               for n in 0:15]

        toric = toric_code_modular_data()
        @test central_charge(toric) == one(ACMG.field(toric))
    end

    @testset "generating function object" begin
        fib = fibonacci_modular_data()
        gf = higher_central_charge_generating_function(fib)
        @test gf.period == 5
        @test length(gf) == 2
        for n in -3:8
            @test gf(n) == higher_central_charge(fib, n)
        end
    end

    @testset "finite-field Fibonacci p=11 regression" begin
        fib = fibonacci_modular_data()
        @test higher_central_charge_sequence_modp(fib, 11; embedding = 5 => 3) ==
              [1, 6, 7, 5, 9]
        @test higher_central_charge(fib, 1; method = :finite_field, p = 11,
                                    embedding = 5 => 3) == 6

        factors = hcc_local_factors(fib, 11; embedding = 5 => 3)
        @test [f.coefficients for f in factors] ==
              [[1, 10], [1, 5], [1, 4], [1, 6], [1, 2]]
        @test hcc_local_factor(fib, 3, 11; embedding = 5 => 3).coefficients == [1, 6]

        @test_throws ErrorException higher_central_charge_sequence_modp(fib, 5)
        @test_throws ErrorException higher_central_charge_sequence_modp(fib, 11;
                                                                        embedding = 5 => 1)

        ctx = CyclotomicContext(2)
        K = ACMG.field(ctx)
        fake = ModularData(ctx, [:one, :x],
                           ACMG.matrix(K, 2, 2, [1, 1, 1, 1]),
                           ACMG.matrix(K, 2, 2, [1, 0, 0, 1]),
                           1, 1, nothing)
        @test_throws ErrorException higher_central_charge_modp(fake, 0, 2)
    end

    @testset "finite-field F/R prototype remains compatible" begin
        sol = fibonacci_fr_data_mod_p(41; solver = :reference, primes = [41, 61])

        @test fr_metadata(sol)[:name] == :fibonacci
        @test fr_metadata(sol)[:p] == 41
        @test fr_metadata(sol)[:conductor] == 20
        @test !isempty(F_values(sol))
        @test !isempty(R_values(sol))

        for n in 1:5
            exact_mod_p = higher_central_charge_modp(fibonacci_modular_data(), n, 41)
            finite_field = higher_central_charge(sol, n)
            @test finite_field.ok
            @test finite_field.value == exact_mod_p
        end
    end

    @testset "FRData finite-field higher central charge is toric-gauge invariant" begin
        fr = fibonacci_fr_data_mod_p(101)
        p = fr_metadata(fr)[:p]
        scalars = Dict(ch => (ch[1] == 1 || ch[2] == 1 ? FpElem(1, p) :
                              FpElem(i + 2, p))
                       for (i, ch) in enumerate(gauge_parameters(fr)))
        moved = apply_gauge(fr, GaugeAction(scalars; field = :F_101))

        @test verify_FRData(moved)
        raw = higher_central_charge(moved, 1)
        fixed = higher_central_charge(fr, 1)
        @test raw.ok
        @test fixed.ok
        @test raw.value == fixed.value
        @test higher_central_charge(moved; n = 1, method = :finite_field).value == raw.value
        @test higher_central_charges(fr, 1:3)[1].value == fixed.value
    end
end
