using ACMG

data = semion_modular_data()
ctx = data.context
p = 17

S_fp = reduce_mod_p(ctx, data.S, p)
T_fp = [reduce_mod_p(ctx, data.T[i, i], p) for i in 1:2]
validated = validate_modular_data(S_fp, T_fp, p, conductor(data))

println("S over F_", p, ": ", S_fp)
println("T over F_", p, ": ", T_fp)
println("F_p modular datum valid: ", validated.valid)
