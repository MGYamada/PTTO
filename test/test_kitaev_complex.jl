"""
Tests for `ACMG.Phase4.KitaevComplex`.

Validates Kitaev 2006 App E.6 gauge analysis (linearised Ocneanu rigidity)
in F-coordinate form on:
  - Fibonacci  (rank 2)
  - Ising      (rank 3)

Expected (from Python prototype):
  Fibonacci:  nF = 4 (in σσσσ-block analogue), gauge_orbit = 1, slice = 3, H³ = 0.
  Ising:      nF = 36 (all entries) / 11 deformable, effective gauge = 1,
              slice = 10, H³ = 0.
"""

using Test
using LinearAlgebra
using ACMG
using ACMG.Phase4

# ============================================================================
# Fibonacci: fusion + base F-symbol
# ============================================================================
#
# 1 ⊗ 1 = 1,  1 ⊗ τ = τ,  τ ⊗ 1 = τ,  τ ⊗ τ = 1 ⊕ τ.
# Object indices (1-based):  1 = 1 (unit),  2 = τ.
#
function fibonacci_fusion()
    r = 2
    N = zeros(Int, r, r, r)
    N[1, 1, 1] = 1
    N[1, 2, 2] = 1
    N[2, 1, 2] = 1
    N[2, 2, 1] = 1
    N[2, 2, 2] = 1
    return FusionRule(N)
end

# Standard Fibonacci F-symbol. F^{τττ}_τ is a 2×2 matrix; all other allowed
# entries are forced to ±1 by unit / structural reasons (here = +1 in our
# convention since we have no nontrivial signs in non-σ context).
function fibonacci_F()
    φ = (1 + sqrt(5)) / 2
    sφ = sqrt(φ)
    F_τττ_τ = [1/φ    1/sφ;
               1/sφ  -1/φ]
    function F(a::Int, b::Int, c::Int, d::Int, e::Int, f::Int)
        # Allowed channel check
        fr = fibonacci_fusion()
        (fr.N[a, b, e] == 1 && fr.N[e, c, d] == 1 &&
         fr.N[b, c, f] == 1 && fr.N[a, f, d] == 1) || return 0.0
        # Any vacuum leg → +1
        1 in (a, b, c, d) && return 1.0
        # a = b = c = d = τ = 2:  use the Hadamard-like 2×2.  e, f ∈ {1, 2}
        return F_τττ_τ[e, f]
    end
    return F
end


# ============================================================================
# Ising: fusion + base F-symbol (standard)
# ============================================================================
#
# 1 = unit, ψ = fermion, σ = anyon.   ψ⊗ψ = 1, ψ⊗σ = σ⊗ψ = σ, σ⊗σ = 1 ⊕ ψ.
# Indices:  1 = 1, 2 = ψ, 3 = σ.
#
function ising_fusion()
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
    return FusionRule(N)
end

function ising_F()
    # Standard Ising:  F^{σσσ}_σ = (1/√2) * [[1, 1], [1, -1]]
    # F^{ψσψ}_σ = -1,  F^{σψσ}_ψ = -1,  F^{σψσ}_1 = +1,  others +1.
    F_sss_s = (1.0 / sqrt(2.0)) * [1.0  1.0;
                                    1.0 -1.0]
    fr = ising_fusion()
    function F(a::Int, b::Int, c::Int, d::Int, e::Int, f::Int)
        (fr.N[a, b, e] == 1 && fr.N[e, c, d] == 1 &&
         fr.N[b, c, f] == 1 && fr.N[a, f, d] == 1) || return 0.0
        1 in (a, b, c, d) && return 1.0
        # σ = 3, ψ = 2.  Count σ's:
        nsig = count(==(3), (a, b, c, d))
        nsig == 0 && return 1.0
        if nsig == 2
            (a, b, c, d) == (2, 3, 2, 3) && return -1.0   # F^{ψσψ}_σ
            (a, b, c, d) == (3, 2, 3, 1) && return  1.0   # F^{σψσ}_1
            (a, b, c, d) == (3, 2, 3, 2) && return -1.0   # F^{σψσ}_ψ
            return 1.0
        end
        if nsig == 4
            d != 3 && return 0.0
            # 2×2 indexed by (e, f) ∈ {1, 2}; map to Julia 1-based with 1->1, 2->2
            return F_sss_s[e, f]
        end
        return 1.0
    end
    return F
end


@testset "KitaevComplex" begin

    # --------------------------------------------------------------------
    @testset "Fibonacci F-coordinate gauge analysis" begin
        fr = fibonacci_fusion()
        F  = fibonacci_F()
        fcs = build_F_coord_space(fr, F)

        # Sanity: F-vars include the σσσσ analogue = (τ,τ,τ,τ,e,f) for e,f ∈ {1,2}
        ττττ_keys = filter(k -> k[1:4] == (2, 2, 2, 2), fcs.vars)
        @test length(ττττ_keys) == 4

        ga = analyze_gauge(fcs)

        # Vertices: (τ,τ;1) and (τ,τ;τ) → 2 vertices
        @test length(ga.verts) == 2

        # The deformable (non-vacuum) F-vars are exactly the 4 (τ,τ,τ,τ,*,*) entries
        non_vacuum = [k for k in fcs.vars if !(1 in k[1:4])]
        @test length(non_vacuum) == 4

        # Effective gauge orbit dim (on F-tangent) = 1
        @test gauge_orbit_dim(ga) == 1

        # Ocneanu rigidity: H³ = 0
        @test H3_dimension(ga) == 0
        @test verify_ocneanu_rigidity(ga; tol = 1e-9)

        # Pentagon-gauge consistency:  Δ_pent · Δ_gauge_eff = 0
        @test opnorm(ga.Delta_pent * ga.Delta_gauge_eff) < 1e-10
    end

    # --------------------------------------------------------------------
    @testset "Ising F-coordinate gauge analysis" begin
        fr = ising_fusion()
        F  = ising_F()
        fcs = build_F_coord_space(fr, F)

        # Total F-vars: 36 (matches Python count)
        @test F_var_count(fcs) == 36

        # Deformable (non-vacuum) F-vars: 11
        non_vacuum = [k for k in fcs.vars if !(1 in k[1:4])]
        @test length(non_vacuum) == 11

        ga = analyze_gauge(fcs)

        # Vertices: (ψ,ψ;1), (ψ,σ;σ), (σ,ψ;σ), (σ,σ;1), (σ,σ;ψ) = 5
        @test length(ga.verts) == 5

        # Effective gauge dim after restricting to unit-preserving subspace = 1
        @test gauge_orbit_dim(ga) == 1

        # H³ = 0
        @test H3_dimension(ga) == 0
        @test verify_ocneanu_rigidity(ga; tol = 1e-9)

        # Pentagon-gauge consistency
        @test opnorm(ga.Delta_pent * ga.Delta_gauge_eff) < 1e-9
    end

    # --------------------------------------------------------------------
    @testset "Slice basis sanity" begin
        for (fr_fn, F_fn, expected_slice_dim) in (
                (fibonacci_fusion, fibonacci_F, 3),   # 4 deformable - 1 gauge = 3
                (ising_fusion,    ising_F,    10),    # 11 deformable - 1 gauge = 10
            )
            fr = fr_fn()
            F  = F_fn()
            fcs = build_F_coord_space(fr, F)
            ga = analyze_gauge(fcs)

            S = slice_basis(ga)

            # Slice basis is orthonormal
            @test opnorm(S' * S - I) < 1e-10

            # Slice is gauge-orthogonal
            @test opnorm(S' * ga.Delta_gauge_eff) < 1e-9

            # Slice is unit-axiom compatible (lives in ker U)
            @test opnorm(ga.Unit * S) < 1e-10

            # Slice dimension matches expectation
            @test slice_dim(ga) == expected_slice_dim
        end
    end
end
