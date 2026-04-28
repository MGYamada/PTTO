using ACMG

result = classify_mtcs_at_conductor(8;
                                    max_rank = 3,
                                    skip_FR = true,
                                    verbose = false)

println("N=8 modular-data results: ", length(result))
for mtc in result
    println(mtc)
end
