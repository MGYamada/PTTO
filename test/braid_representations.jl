using Test
using ACMG
using Oscar

function _eval_tc_poly_at(eq, values)
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

@testset "braid representation fusion bases" begin
    fib = fibonacci_fusion_rules()
    @test fusion_paths(fib, [:τ, :τ, :τ], :τ) == fusion_paths(fib, [:τ, :τ, :τ], :τ)
    @test ACMG.dim(fusion_basis(fib, [:τ, :τ, :τ], :τ)) == 2
    @test ACMG.dim(fusion_basis(fib, [:τ, :τ, :τ, :τ], :one)) == 2
    @test isempty(fusion_paths(fib, [:one, :one], :τ))

    semion = semion_fusion_rules()
    @test ACMG.dim(fusion_basis(semion, [:s, :s, :s], :s)) == 1
end

@testset "built-in FRData solves TensorCategories equations" begin
    cases = [
        semion_fr_data(),
        fibonacci_fr_data(),
        ising_fr_data(),
    ]
    for fr_data in cases
        rules = fr_data.rules
        pent = pentagon_equations(rules)
        hex = hexagon_equations(rules)
        @test fr_data isa FRData
        @test !hasproperty(fr_data, :F)
        @test !hasproperty(fr_data, :R)
        @test length(fr_data.F_values) == get_pentagon_system(rules.N, rules.rank)[3]
        @test length(ACMG._fr_hexagon_values(fr_data)) == get_hexagon_fr_system(rules.N, rules.rank)[3]
        @test all(iszero(_eval_tc_poly_at(eq, ACMG._fr_pentagon_values(fr_data))) for eq in pent)
        @test all(iszero(_eval_tc_poly_at(eq, ACMG._fr_hexagon_values(fr_data))) for eq in hex)
    end
end

@testset "braid generator construction" begin
    br = braid_representation(semion_fr_data(), [:s, :s, :s], :s)
    @test length(br.generators) == 2
    @test size(br.generators[1]) == (1, 1)

    fib = braid_representation(fibonacci_fr_data(), [:τ, :τ, :τ], :τ)
    @test length(fib.generators) == 2
    @test size(fib.generators[1]) == (2, 2)

    fib_data = fibonacci_fr_data()
    fib_tuple = (rules = fib_data.rules,
                 F_values = fib_data.F_values,
                 R_values = fib_data.R_values,
                 R_inverse_values = fib_data.R_inverse_values)
    fib_from_tuple = braid_representation(fib_tuple, [:τ, :τ, :τ], :τ)
    @test check_braid_relations(fib_from_tuple).ok

    bad = zeros(Int, 2, 2, 2)
    bad[1, 1, 1] = 1
    bad[1, 2, 2] = 1
    bad[2, 1, 2] = 1
    bad[2, 2, 1] = 2
    @test_throws ErrorException fusion_basis(FusionRule(bad), [2, 2], 1)
end
