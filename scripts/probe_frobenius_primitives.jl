"""
scripts/probe_frobenius_primitives.jl

Check what primitives TensorCategories.jl exposes for the canonical
Frobenius algebra construction needed for EGNO Eq. 9.2 implementation:

  A = Hom_{C ⊠ C^op}(1, 1)   (canonical Frobenius algebra of C)
  m : A ⊗ A → A              (multiplication)
  Δ : A → A ⊗ A              (comultiplication)
  e : 1 → A                  (unit)
  u := m ∘ Δ ∘ e : 1 → A     (central element)
  u⁻¹ : 1 → A                (convolution inverse)

  χ(f) = m ∘ (id ⊗ f) ∘ (Δ ∘ u⁻¹ ⊗ id)

Key questions:
  1. Is there a constructor for the canonical Frobenius algebra object?
  2. Are `m`, `Δ`, `e`, `u` available as morphisms?
  3. Alternatively, is there a direct `chi` or `Yetter_cohomology` primitive?
  4. If not, can we build it ourselves from `⊗`, `∘`, `ev`, `coev`?

Also inspect morphism construction from matrix data — to build f ∈ C^n
directly (without ev/coev triggering divexact failure).
"""

using TensorCategories
using Oscar

function ising_Nijk()
    r = 3
    N = zeros(Int, r, r, r)
    for j in 1:r
        N[1, j, j] = 1; N[j, 1, j] = 1
    end
    N[2,2,1]=1; N[2,3,3]=1; N[3,2,3]=1; N[3,3,1]=1; N[3,3,2]=1
    return N
end

println("="^68)
println("Probing TensorCategories for Frobenius primitives")
println("="^68)

C = TensorCategories.six_j_category(QQ, ising_Nijk())
C.one = [1, 0, 0]

# Check for Frobenius-related names
for fn in [:frobenius, :frobenius_algebra, :canonical_frobenius, :frobenius_algebra_of,
          :multiplication, :comultiplication, :unit_morphism, :counit,
          :m, :Δ, :delta, :e, :u, :convolution, :convolution_inverse,
          :yetter_complex, :davydov_yetter, :yetter_cohomology,
          :chi, :contracting_homotopy, :deformation_complex,
          :internal_hom, :End, :endomorphism_ring,
          :tensor_unit, :one, :one_object,
          :tr, :trace, :trace_morphism, :left_dual, :right_dual]
    try
        ms = methods(getfield(TensorCategories, fn))
        if length(ms) > 0
            println("$fn  →  $(length(ms)) methods")
            # Show first 2 signatures
            for (i, m) in enumerate(ms)
                i > 2 && break
                println("    $m")
            end
        end
    catch e
        # silent — means name not exported
    end
end

println("\n--- Constructing morphisms directly from matrix data ---")
println("\nSixJMorphism fields:")
id1 = id(simples(C)[1])
println("  fieldnames(id(1)): $(fieldnames(typeof(id1)))")
println("  id(1): $id1")

println("\nCan we build a custom morphism?")
σ = simples(C)[3]
println("  σ ⊗ σ = $(σ ⊗ σ)")
# Check the type of id(σ ⊗ σ)
idσσ = id(σ ⊗ σ)
println("  id(σ ⊗ σ) matrices: $(idσσ.m)")  # .m is likely the matrix list
