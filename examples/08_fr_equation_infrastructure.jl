using ACMG

function summarize_fr_example(name, rules; p = 11)
    fvars = fsymbol_variables(rules)
    rvars = rsymbol_variables(rules)
    pent = pentagon_equations(rules, fvars)
    hex = hexagon_equations(rules, fvars, rvars)
    system = fr_equation_system(rules)
    fixed = gauge_fix(system; strategy = :safe)
    fp_system = reduce_mod_p(fixed, p)

    println(name)
    println("  F variables: ", length(fvars))
    println("  R variables: ", length(rvars))
    println("  pentagon equations: ", length(pent))
    println("  hexagon equations: ", length(hex))
    println("  fixed equations: ", length(fixed.equations))
    println("  finite field: F_", fp_system.p)
end

summarize_fr_example("Semion", semion_fusion_rules())
summarize_fr_example("Toric code", toric_code_fusion_rules())
summarize_fr_example("Fibonacci", fibonacci_fusion_rules())
summarize_fr_example("Ising", ising_fusion_rules())
