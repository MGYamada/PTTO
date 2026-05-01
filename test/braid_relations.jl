using Test
using ACMG

@testset "braid relations" begin
    semion = braid_representation(semion_fr_data_mod_p(17), [2, 2, 2], 2)
    @test check_braid_relations(semion).ok

    fib = braid_representation(fibonacci_fr_data_mod_p(101), [2, 2, 2], 2)
    @test check_braid_relations(fib).ok
end
