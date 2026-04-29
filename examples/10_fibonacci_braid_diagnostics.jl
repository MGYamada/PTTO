using ACMG

br = braid_representation(fibonacci_fr_data(), [:τ, :τ, :τ], :τ)
println("basis dimension: ", dim(br.basis))
println("braid relations: ", check_braid_relations(br).ok)

brp = reduce_mod_p(br, 41; conductor = 20)
println(generated_matrix_algebra(brp))
println(commutant(brp))
