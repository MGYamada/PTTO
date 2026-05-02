# Experimental gauge-stabilizer example.
#
# v0.9.1 computes Aut(x) = Stab_G(x) for already available solution data.
# The resulting order is future input for stacky weights such as 1 / |Aut(x)|.

using ACMG

function main()
    p = 5
    G = [ACMG.FpElem(i, p) for i in 1:(p - 1)]
    action = (x, g) -> g * x

    x = ACMG.FpElem(2, p)
    result = stabilizer(x, G; action = action)
    order = stabilizer_order(result)

    println("Aut(x) order = ", order)
    println("future stacky weight = ", 1 // order)

    zero = ACMG.FpElem(0, p)
    zero_result = stabilizer(zero, G; action = action)
    zero_order = stabilizer_order(zero_result)

    println("Aut(0) order = ", zero_order)
    println("future stacky weight at 0 = ", 1 // zero_order)
end

main()
