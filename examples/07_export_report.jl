using ACMG

result = classify_mtcs_at_conductor(8;
                                    max_rank = 2,
                                    primes = [41, 73],
                                    skip_FR = false,
                                    verbose = false)

save_classification("N8.json", result)
write_report("N8.md", result)

println("Wrote N8.json and N8.md")
