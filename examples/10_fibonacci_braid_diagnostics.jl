using ACMG

br = braid_representation(fibonacci_fr_data_mod_p(101), [2, 2, 2], 2)
println("basis dimension: ", dim(br.basis))
println("braid relations: ", check_braid_relations(br).ok)

brp = reduce_mod_p(br)
println(generated_matrix_algebra(brp))
println(commutant(brp))
