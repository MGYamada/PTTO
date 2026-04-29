using ACMG

br = braid_representation(ising_fr_data(), [:σ, :σ, :σ], :σ)
println("basis dimension: ", dim(br.basis))
println("braid relations: ", check_braid_relations(br).ok)

brp = reduce_mod_p(br, 17; conductor = 16)
println(finite_group_diagnostics(brp; max_size = 1000))
