using ACMG

br = braid_representation(semion_fr_data_mod_p(17), [2, 2, 2], 2)
println("basis dimension: ", dim(br.basis))
println("braid relations: ", check_braid_relations(br).ok)

brp = reduce_mod_p(br)
println(finite_group_diagnostics(brp; max_size = 1000))
println(zariski_closure_diagnostics(brp; max_words = 50))
