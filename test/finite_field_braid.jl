using Test
using ACMG

@testset "finite-field braid reduction" begin
    br = braid_representation(semion_fr_data(), [:s, :s, :s], :s)
    brp = reduce_mod_p(br, 17; conductor = 8)
    @test check_braid_relations(brp).ok
    gdiag = finite_group_diagnostics(brp; max_size = 1000)
    @test gdiag.dimension == 1
    @test gdiag.group_size !== nothing

    fib = braid_representation(fibonacci_fr_data(), [:τ, :τ, :τ], :τ)
    fibp = reduce_mod_p(fib, 41; conductor = 20)
    @test check_braid_relations(fibp).ok

    @test_throws ErrorException reduce_mod_p(br, 5; conductor = 8)
end
