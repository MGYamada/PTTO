using Test
using ACMG

@testset "Zariski diagnostics over finite fields" begin
    br = braid_representation(semion_fr_data_mod_p(17), [2, 2, 2], 2)
    brp = reduce_mod_p(br)
    alg = generated_matrix_algebra(brp)
    @test alg.dimension == 1
    @test alg.is_full_matrix_algebra

    c = commutant(brp)
    @test c.dimension == 1

    z = zariski_closure_diagnostics(brp; max_words = 20, max_degree = 2)
    @test z.diagnostic_only
    @test !z.exact_closure_computed
    @test z.matrix_algebra_dimension == 1
end
