using Test
using ACMG

@testset "finite-field braid reduction" begin
    br = braid_representation(semion_fr_data_mod_p(17), [2, 2, 2], 2)
    brp = reduce_mod_p(br)
    @test typeof(brp).parameters[1] == typeof(br)
    @test eltype(brp.generators) == Matrix{Int}
    @test check_braid_relations(brp).ok
    gdiag = finite_group_diagnostics(brp; max_size = 1000)
    @test gdiag.dimension == 1
    @test gdiag.group_size !== nothing

    fib = braid_representation(fibonacci_fr_data_mod_p(101), [2, 2, 2], 2)
    fibp = reduce_mod_p(fib)
    @test check_braid_relations(fibp).ok
end
