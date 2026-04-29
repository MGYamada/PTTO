using ACMG

br = braid_representation(semion_fr_data(), [:s, :s, :s], :s)
println("basis dimension: ", dim(br.basis))
println("braid relations: ", check_braid_relations(br).ok)

brp = reduce_mod_p(br, 17; conductor = 8)
println(finite_group_diagnostics(brp; max_size = 1000))
println(zariski_closure_diagnostics(brp; max_words = 50))
