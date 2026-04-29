using Test
using ACMG

@testset "v0.8 FR equations" begin
    semion = semion_fusion_rules()
    fib = fibonacci_fusion_rules()
    ising = ising_fusion_rules()

    @test simple_objects(fib) == [1, 2]
    @test fusion_channels(fib, 2, 2) == [1, 2]
    @test (2, 2, 1) in ACMG.admissible_triples(fib)
    @test !isempty(ACMG.admissible_quadruples(fib))

    @test_throws UndefVarError fsymbol_variables(semion)
    @test_throws UndefVarError rsymbol_variables(semion)

    pent = pentagon_equations(fib)
    hex = hexagon_equations(fib)
    @test !isempty(pent)
    @test !isempty(hex)
    @test pent == get_pentagon_system(fib.N, fib.rank)[2]
    @test hex == get_hexagon_fr_system(fib.N, fib.rank)[2]

    isys = fr_equation_system(ising)
    @test validate_fr_system(isys)
    @test isys.metadata[:f_variables] == get_pentagon_system(ising.N, ising.rank)[3]
    @test isys.metadata[:r_variables] == sum(ising.N[a, b, c]^2 for a in 1:ising.rank, b in 1:ising.rank, c in 1:ising.rank)
    @test isys.metadata[:include_hexagon]
end
