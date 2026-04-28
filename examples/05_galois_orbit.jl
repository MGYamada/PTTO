using ACMG

data = fibonacci_modular_data()
orbit = galois_orbit(data)

println("Fibonacci Galois orbit size: ", length(orbit))
for (i, item) in enumerate(orbit)
    println("branch ", i, ": T22 = ", item.T[2, 2])
end
