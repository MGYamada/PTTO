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

    sf = fsymbol_variables(semion)
    sr = rsymbol_variables(semion)
    @test !isempty(pentagon_equations(semion, sf))
    @test !isempty(hexagon_equations(semion, sf, sr))

    ff = fsymbol_variables(fib)
    fr = rsymbol_variables(fib)
    pent = pentagon_equations(fib, ff)
    hex = hexagon_equations(fib, ff, fr)
    @test !isempty(pent)
    @test !isempty(hex)
    @test any(eq -> eq.metadata[:kind] == :pentagon, pent)

    isys = fr_equation_system(ising)
    @test validate_fr_system(isys)
    @test isys.metadata[:f_variables] == length(isys.fvars)
    @test isys.metadata[:r_variables] == length(isys.rvars)
end
