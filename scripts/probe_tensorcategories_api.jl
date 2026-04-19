"""
scripts/probe_tensorcategories_api.jl

Probe the TensorCategories.jl API to identify which morphism helpers
(dual object, evaluation, coevaluation, quantum trace) are available
for building χ³ via categorical composition, in the same style as
`pentagon_equations`.

Output guides the design of `chi3_constraints` written categorically
(not as a numerical matrix).
"""

using TensorCategories
using Oscar

function ising_Nijk()
    r = 3
    N = zeros(Int, r, r, r)
    for j in 1:r
        N[1, j, j] = 1
        N[j, 1, j] = 1
    end
    N[2, 2, 1] = 1
    N[2, 3, 3] = 1
    N[3, 2, 3] = 1
    N[3, 3, 1] = 1
    N[3, 3, 2] = 1
    return N
end

Nijk = ising_Nijk()
one_vec = [1, 0, 0]

# Build the six_j_category with a polynomial base ring, same way
# pentagon_equations does internally.
println("="^68)
println("Building six_j_category with polynomial base")
println("="^68)

# Option A: use QQ
C_qq = TensorCategories.six_j_category(QQ, Nijk)
C_qq.one = one_vec
println("Type of C_qq: ", typeof(C_qq))
println("Fieldnames: ", fieldnames(typeof(C_qq)))

println("\n--- Methods on the category ---")
# Try common categorical helpers
for fn in [:dual, :ev, :coev, :evaluation, :coevaluation,
          :quantum_trace, :partial_trace, :ptrace, :trace,
          :associator, :braiding]
    try
        ms = methods(getfield(TensorCategories, fn))
        println("$fn  →  $(length(ms)) methods")
    catch e
        println("$fn  →  NOT EXPORTED (or unknown)")
    end
end

# Check whether simple objects have a `dual` field or method
try
    obj = simples(C_qq)[2]       # ψ or σ
    println("\nA simple object type: $(typeof(obj))")
    println("Fieldnames: $(fieldnames(typeof(obj)))")
    # Try to get its dual
    try
        d = dual(obj)
        println("dual(obj) = $d  OK")
    catch e
        println("dual(obj) failed: $e")
    end
catch e
    println("simples(C_qq) failed: $e")
end

# Inspect pentagon_equations source to understand the morphism machinery
println("\n--- Look at pentagon_equations method ---")
try
    for m in methods(TensorCategories.pentagon_equations)
        println("  $m")
    end
catch e
    println("  (methods lookup failed: $e)")
end

# Associator access
println("\n--- Associator access ---")
println("C_qq.ass type: ", typeof(C_qq.ass))
println("Size: ", size(C_qq.ass))
println("C_qq.ass[3,3,3,3] type: ", typeof(C_qq.ass[3,3,3,3]))
println("C_qq.ass[3,3,3,3] size: ", size(C_qq.ass[3,3,3,3]))
