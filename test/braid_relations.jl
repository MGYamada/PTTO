using Test
using ACMG

@testset "braid relations" begin
    semion = braid_representation(semion_fr_data(), [:s, :s, :s], :s)
    @test check_braid_relations(semion).ok

    fib = braid_representation(fibonacci_fr_data(), [:τ, :τ, :τ], :τ)
    @test check_braid_relations(fib).ok
end
