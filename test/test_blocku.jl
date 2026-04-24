using Test
using ACMG

"""
Tests for BlockU module — Phase 2 core.

These tests are split into two groups:
1. Pure F_p tests (no Oscar needed): circle enumeration, apply_o2_block,
   verlinde_find_unit, t_eigenspace_decomposition, parameter_dim,
   signed_Fp.
2. Integration tests with Oscar (cyclotomic_to_Fp, find_zeta_in_Fp,
   build_block_diagonal, find_mtcs_at_prime): these require an atomic
   catalog.

Group 1 is implemented here; Group 2 is in test_blocku_oscar.jl (TODO).
"""

@testset "BlockU pure F_p helpers" begin

    @testset "signed_Fp" begin
        @test signed_Fp(0, 7) == 0
        @test signed_Fp(1, 7) == 1
        @test signed_Fp(3, 7) == 3
        @test signed_Fp(4, 7) == -3   # 4 > 7/2 → 4 - 7 = -3
        @test signed_Fp(6, 7) == -1
        @test signed_Fp(36, 73) == 36
        @test signed_Fp(37, 73) == -36
        @test signed_Fp(52, 73) == -21
    end

    @testset "find_zeta_in_Fp" begin
        # N=5, p=41: 40 = 8·5, so ζ_5 exists
        zeta5 = find_zeta_in_Fp(5, 41)
        @test powermod(zeta5, 5, 41) == 1
        @test powermod(zeta5, 1, 41) != 1  # primitive
        # Known value from test_fibonacci: ζ_5 = 10 at p=41
        @test zeta5 == 10

        # N=16, p=17: 16 | 16, so ζ_16 exists
        zeta16 = find_zeta_in_Fp(16, 17)
        @test powermod(zeta16, 16, 17) == 1
        # Known: ζ_16 = 3 at p=17
        @test zeta16 == 3

        # N=24, p=73
        zeta24 = find_zeta_in_Fp(24, 73)
        @test powermod(zeta24, 24, 73) == 1
        # Known: ζ_24 = 52 or another primitive (depends on primitive_root choice)
        # We only check it's a primitive 24th root
        for d in [2, 3, 4, 6, 8, 12]
            @test powermod(zeta24, d, 73) != 1
        end

        # Error cases
        @test_throws ErrorException find_zeta_in_Fp(7, 13)   # 7 ∤ 12
    end

    @testset "t_eigenspace_decomposition" begin
        # Simple case: no degeneracy
        T = [1, 2, 3, 4]
        groups = t_eigenspace_decomposition(T, 7)
        @test groups[1] == [1]
        @test groups[2] == [2]
        @test groups[3] == [3]
        @test groups[4] == [4]
        @test parameter_dim(groups) == 0

        # One 2-fold degeneracy
        T = [1, 1, 2, 3]
        groups = t_eigenspace_decomposition(T, 7)
        @test groups[1] == [1, 2]
        @test groups[2] == [3]
        @test groups[3] == [4]
        @test parameter_dim(groups) == 1  # C(2, 2) = 1

        # Two 2-fold degeneracies
        T = [1, 1, 2, 2, 3]
        groups = t_eigenspace_decomposition(T, 7)
        @test parameter_dim(groups) == 2

        # One 3-fold degeneracy
        T = [1, 1, 1, 2, 3]
        groups = t_eigenspace_decomposition(T, 7)
        @test parameter_dim(groups) == 3  # C(3, 2) = 3
    end

    @testset "o2_circle_points (count)" begin
        # |O(2)(F_p)| = 2(p - ε) where ε = ±1 depending on p mod 4
        # Circle count |{u² + v² = 1}| = p - ε
        p = 13  # p ≡ 1 (mod 4)
        circle = o2_circle_points(p)
        @test length(circle) == 12  # p - 1

        p = 11  # p ≡ 3 (mod 4)
        circle = o2_circle_points(p)
        @test length(circle) == 12  # p + 1

        p = 73  # p ≡ 1 (mod 4)
        circle = o2_circle_points(p)
        @test length(circle) == 72  # p - 1

        # Every point satisfies u² + v² = 1
        for (u, v) in o2_circle_points(17)
            @test mod(u*u + v*v, 17) == 1
        end
    end

    @testset "apply_o2_block (identity check)" begin
        # Applying U = I (det=+1, u=1, v=0) should leave S unchanged
        p = 13
        S = [1 2 3 4; 2 5 6 7; 3 6 8 9; 4 7 9 10]
        S_out = apply_o2_block(S, (1, 2), 1, 0, +1, p)
        @test S_out == S

        # Applying rotation π/2: (u, v) = (0, 1), det = +1
        # U = [[0, -1], [1, 0]] on positions (1, 2)
        S_rot = apply_o2_block(S, (1, 2), 0, 1, +1, p)
        # This rotates cols/rows 1, 2; check result is still symmetric if S was
        @test S_rot[1, 1] == S[2, 2]  # rotation swaps (1,1) and (2,2)
        @test S_rot[2, 2] == S[1, 1]
        # Off-diagonal: S_rot[1, 2] = -S[1, 2]
        @test mod(S_rot[1, 2] + S[1, 2], p) == 0
    end

    @testset "verlinde_find_unit on Fibonacci" begin
        # Hand-built Fibonacci in F_41 from test_fibonacci
        # S = [14 16; 16 27], expected unit = 1 (index 1), d = [1, φ]
        p = 41
        S = [14 16; 16 27]
        result = verlinde_find_unit(S, p; threshold = 3)
        @test result !== nothing
        if result !== nothing
            (u, N) = result
            @test u == 1  # first index should be unit

            # Fibonacci fusion: ψ ⊗ ψ = 1 ⊕ ψ
            # N[1][1][:] should be (1, 0) (unit × unit = unit)
            @test N[1, 1, 1] == 1
            @test N[1, 1, 2] == 0
            # N[2][2][:] should be (1, 1)
            @test N[2, 2, 1] == 1
            @test N[2, 2, 2] == 1
        end
    end

    @testset "passes_unit_axiom" begin
        p = 41
        S = [14 16; 16 27]  # Fibonacci
        @test passes_unit_axiom(S, p, 1)
        @test !passes_unit_axiom(S, p, 2)
    end

    @testset "verlinde_find_unit on Ising" begin
        # Ising in F_17 from test_ising: S = [9 14 9; 14 0 3; 9 3 9]
        # Expected unit = 1, fusion structure Ising
        p = 17
        S = [9 14 9; 14 0 3; 9 3 9]
        result = verlinde_find_unit(S, p; threshold = 3)
        @test result !== nothing
        if result !== nothing
            (u, N) = result
            @test u == 1
            # Ising: σ ⊗ σ = 1 ⊕ ψ
            # obj 1 = 1 (unit), obj 2 = σ (d=√2), obj 3 = ψ (d=1)
            # N[2][2][:] should be (1, 0, 1)
            @test N[2, 2, 1] == 1
            @test N[2, 2, 2] == 0
            @test N[2, 2, 3] == 1
            # ψ ⊗ ψ = 1
            @test N[3, 3, 1] == 1
            @test N[3, 3, 2] == 0
            @test N[3, 3, 3] == 0
        end
    end

    @testset "MTCCandidate show" begin
        # Just check it doesn't throw
        N = zeros(Int, 2, 2, 2)
        # New format: U_params can be Matrix or Tuple
        c1 = MTCCandidate(41, [1 0; 0 1], [14 16; 16 27], [1, 18], 1, N, [1, 7], 9)
        s1 = string(c1)
        @test occursin("p=41", s1)
        @test occursin("unit=1", s1)

        # Legacy Tuple format still works
        c2 = MTCCandidate(41, (1, 0, +1), [14 16; 16 27], [1, 18], 1, N, [1, 7], 9)
        s2 = string(c2)
        @test occursin("p=41", s2)
    end
end
