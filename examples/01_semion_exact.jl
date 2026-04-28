using ACMG

data = semion_modular_data()
validation = validate_exact_modular_data(data)

println("Semion conductor: ", conductor(data))
println("Exact modular data valid: ", validation.ok)
