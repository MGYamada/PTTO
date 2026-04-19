"""
scripts/probe_6j_solver.jl

Read the actual `pentagon_equations` source and understand:
  1. How the polynomial ring is constructed (how are variables assigned?)
  2. How the associator is turned into a symbolic MatElem
  3. How f = g pentagon equation reduces to polynomial rows
  4. What `ev` / `coev` actually do on SixJCategory

This tells us exactly how to write `chi3_equations` in the same style.
"""

using TensorCategories
using Oscar

# Locate and cat the source file
src_path = "~/.julia/packages/TensorCategories/l1fLP/src/TensorCategoryFramework/6j-Solver.jl"
expanded = expanduser(src_path)
println("="^68)
println("Source: $expanded")
println("="^68)
try
    println(read(expanded, String))
catch e
    println("Read failed: $e")
end

println("\n" * "="^68)
println("ev / coev method signatures")
println("="^68)
for m in methods(TensorCategories.ev)
    println("  $m")
end
println()
for m in methods(TensorCategories.coev)
    println("  $m")
end

println("\n" * "="^68)
println("Do ev / coev work on a simple SixJObject?")
println("="^68)
function ising_Nijk()
    r = 3
    N = zeros(Int, r, r, r)
    for j in 1:r
        N[1, j, j] = 1; N[j, 1, j] = 1
    end
    N[2,2,1]=1; N[2,3,3]=1; N[3,2,3]=1; N[3,3,1]=1; N[3,3,2]=1
    return N
end

C = TensorCategories.six_j_category(QQ, ising_Nijk())
C.one = [1, 0, 0]
σ = simples(C)[3]
println("σ = $σ,  dual(σ) = $(dual(σ))")

try
    e = ev(σ)
    println("\nev(σ) type: $(typeof(e))")
    println("ev(σ): $e")
catch err
    println("\nev(σ) failed: $err")
end

try
    c = coev(σ)
    println("\ncoev(σ) type: $(typeof(c))")
    println("coev(σ): $c")
catch err
    println("\ncoev(σ) failed: $err")
end

# Try a morphism composition expressing the cap: id_X ⊗ ev_σ  on  X ⊗ σ ⊗ σ^*
# This is what we need for χ³.
try
    # Build X ⊗ σ ⊗ dual(σ)  = X ⊗ σ ⊗ σ (since self-dual)
    ψ = simples(C)[2]
    dual_σ = dual(σ)
    # ev_σ : dual(σ) ⊗ σ → 1   (or σ ⊗ dual(σ) → 1 depending on convention)
    # Let's see what TensorCategories' ev gives us.
    e = ev(σ)
    println("\ndomain(ev(σ)) = $(domain(e))")
    println("codomain(ev(σ)) = $(codomain(e))")
catch err
    println("\ndomain/codomain failed: $err")
end
