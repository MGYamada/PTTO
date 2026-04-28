using ACMG

function ising_fusion()
    N = zeros(Int, 3, 3, 3)
    N[1, 1, 1] = 1
    for i in 1:3
        N[1, i, i] = 1
        N[i, 1, i] = 1
    end
    N[2, 2, 1] = 1
    N[2, 2, 3] = 1
    N[2, 3, 2] = 1
    N[3, 2, 2] = 1
    N[3, 3, 1] = 1
    return N
end

data = ising_modular_data()
twists = [data.T[i, i] for i in 1:3]
fr = compute_FR_from_ST(ising_fusion();
                        conductor = 16,
                        primes = [17, 97, 113],
                        S = data.S,
                        T = twists)

report = fr.report
println("Ising roundtrip ok: ", report.ok)
println("Best permutation: ", report.best_perm)
