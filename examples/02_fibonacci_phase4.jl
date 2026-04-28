using ACMG

function fibonacci_fusion()
    N = zeros(Int, 2, 2, 2)
    N[1, 1, 1] = 1
    N[1, 2, 2] = 1
    N[2, 1, 2] = 1
    N[2, 2, 1] = 1
    N[2, 2, 2] = 1
    return N
end

data = fibonacci_modular_data()
twists = [data.T[i, i] for i in 1:2]
fr = compute_FR_from_ST(fibonacci_fusion();
                        conductor = 20,
                        primes = [41, 61],
                        S = data.S,
                        T = twists)

println("F entries: ", length(fr.F))
println("R entries: ", length(fr.R))
println("Roundtrip: ", fr.report === nothing ? "not scored" : fr.report.ok)
