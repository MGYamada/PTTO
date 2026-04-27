using Test
using ACMG
using Oscar

"""
Tests for general O(n) Cayley parametrisation in BlockU.jl.
"""

# Helper: 2x2 determinant (defined at top-level so @testset blocks see it)
det_2x2(M) = M[1, 1] * M[2, 2] - M[1, 2] * M[2, 1]

@testset "General O(n) Cayley" begin

    @testset "inverse_mod_p" begin
        # 2x2 example
        p = 7
        M = [3 1; 2 4]
        Minv = inverse_mod_p(M, p)
        @test Minv !== nothing
        # Verify: M * Minv ≡ I
        prod_M = [mod(sum(M[i, k] * Minv[k, j] for k in 1:2), p) for i in 1:2, j in 1:2]
        @test prod_M == [1 0; 0 1]

        # 3x3 identity
        I3 = [1 0 0; 0 1 0; 0 0 1]
        @test inverse_mod_p(I3, 11) == I3

        # Singular matrix
        sing = [1 2; 2 4]
        @test inverse_mod_p(sing, 7) === nothing
    end

    @testset "cayley_so_n at n=2 reproduces O(2) circle points" begin
        # Cayley with single param a → SO(2) rotation.
        # A = [[0, a], [-a, 0]], U = (I-A)(I+A)^{-1}
        p = 13
        for a in 0:(p-1)
            U = cayley_so_n([a], 2, p)
            U === nothing && continue  # 1+a² ≡ 0 mod p (singular)
            # U should be a 2x2 SO(2) matrix
            @test size(U) == (2, 2)
            @test mod(det_2x2(U), p) == 1  # det = +1
            # U U^T = I
            UUT = [mod(sum(U[i, k] * U[j, k] for k in 1:2), p) for i in 1:2, j in 1:2]
            @test UUT == [1 0; 0 1]
        end
    end

    @testset "enumerate_so_n_Fp counts" begin
        # SO(1) = {1}
        @test length(enumerate_so_n_Fp(1, 7)) == 1

        # Cayley map on SO(2) covers p - #{a: 1+a²=0} points.
        # For p=13 (≡1 mod 4): 1+a²=0 at a=±5, so 2 excluded → 13 - 2 = 11
        # For p=11 (≡3 mod 4): no solutions → 11
        so2_p13 = enumerate_so_n_Fp(2, 13)
        @test length(so2_p13) == 11

        so2_p11 = enumerate_so_n_Fp(2, 11)
        @test length(so2_p11) == 11

        # All elements should be valid SO(2)
        for U in so2_p13
            @test mod(det_2x2(U), 13) == 1
            UUT = [mod(sum(U[i, k] * U[j, k] for k in 1:2), 13) for i in 1:2, j in 1:2]
            @test UUT == [1 0; 0 1]
        end
    end

    @testset "enumerate_o_n_Fp size" begin
        p = 13
        o2 = enumerate_o_n_Fp(2, p)
        so2 = enumerate_so_n_Fp(2, p)
        @test length(o2) == 2 * length(so2)
        # Half should have det = +1, other half det = -1
        det_pos = count(U -> mod(det_2x2(U), p) == 1, o2)
        det_neg = count(U -> mod(det_2x2(U), p) == p - 1, o2)
        @test det_pos == length(so2)
        @test det_neg == length(so2)
    end

    @testset "search backend mode parity" begin
        p = 11
        for n in (2, 3)
            cands_exhaustive = enumerate_block_candidates(n, p, :exhaustive)
            cands_groebner = enumerate_block_candidates(n, p, :groebner)
            @test !isempty(cands_groebner)
            @test length(cands_groebner) >= length(cands_exhaustive)
            # Determinism check across repeated calls.
            cands_groebner2 = enumerate_block_candidates(n, p, :groebner)
            @test cands_groebner2 == cands_groebner
        end
        @test_throws ErrorException enumerate_block_candidates(2, p, :unknown_mode)
    end

    @testset "validate_search_mode" begin
        @test validate_search_mode(:exhaustive) === nothing
        @test validate_search_mode(:groebner) === nothing
        @test_throws ErrorException validate_search_mode(:bad_mode)
    end

    @testset "verlinde groebner equation builder shape" begin
        # Pure shape/count sanity using small finite-field matrices.
        # n_block=2 => orthogonality equations: 3
        # r=2 => inverse witness equations: 2
        # r=2 => unit axiom equations: 4
        # total expected = 9
        varsU = [1, 0, 0, 1]
        varsW = [1, 1]
        S = [1 0; 0 1]
        eqs = build_verlinde_unit_equations(varsU, varsW, S, [1, 2], 11, 1)
        @test length(eqs) == 9
    end

    @testset "cayley link equation builder shape" begin
        # n=2 => one Cayley parameter and 2x2 U variables.
        varsA = [0]
        varsU = [1, 0, 0, 1]
        eqs = build_cayley_link_equations(varsA, varsU, 2)
        @test length(eqs) == 4
    end

    @testset "solve_cayley_unit_filtered_blocks" begin
        p = 41
        S = [14 16; 16 27]  # Fibonacci S in F_41
        blocks = solve_cayley_unit_filtered_blocks(S, [1, 2], p; max_units = 2)
        @test !isempty(blocks)
        @test all(U -> is_orthogonal_mod_p(U, p), blocks)
    end

    @testset "apply_block_U sanity" begin
        # Note: apply_block_U reduces entries mod p. For identity-block test,
        # keep all S entries already in [0, p) so that S and apply_block_U(S,...,I,...)
        # compare equal without mod reduction.
        p = 13
        S = [1 2 3 4; 2 5 6 7; 3 6 8 9; 4 7 9 10]  # all entries < 13
        U_id = [1 0; 0 1]
        @test apply_block_U(S, [1, 2], U_id, p) == S

        # 5x5 with 3x3 block test — keep entries in [0, p)
        S5 = [1 2 3 4 5;
              2 6 7 8 9;
              3 7 10 11 12;
              4 8 11 0 1;
              5 9 12 1 2]
        U_id3 = [1 0 0; 0 1 0; 0 0 1]
        @test apply_block_U(S5, [1, 3, 4], U_id3, p) == S5

        # Non-trivial: 90° rotation on (1, 2) block should permute rows/cols 1, 2
        # U = [[0, -1], [1, 0]] (mod p: [[0, p-1], [1, 0]])
        U_rot = [0 (p-1); 1 0]
        S_rot = apply_block_U(S, [1, 2], U_rot, p)
        @test size(S_rot) == (4, 4)
        # Entries outside the (1,2) block should be unchanged
        @test S_rot[3, 3] == S[3, 3]
        @test S_rot[3, 4] == S[3, 4]
        @test S_rot[4, 4] == S[4, 4]
    end

    @testset "multiple degenerate eigenspace product helpers" begin
        p = 13
        S = [1 2 3 4;
             2 5 6 7;
             3 6 8 9;
             4 7 9 10]
        U_left = [0 (p - 1); 1 0]
        U_right = [0 1; 1 0]
        choices = [
            (theta = 1, indices = [1, 2], U = U_left),
            (theta = 2, indices = [3, 4], U = U_right),
        ]

        sequential = apply_block_U(apply_block_U(S, [1, 2], U_left, p),
                                   [3, 4], U_right, p)
        @test ACMG._apply_block_U_product(S, choices, p) == sequential
    end

    @testset "find_mtcs_at_prime accepts multiple degenerate eigenspaces" begin
        K, _ = cyclotomic_field(1)
        S = identity_matrix(K, 4)
        T = diagonal_matrix(K, [K(1), K(1), K(2), K(2)])
        atom = AtomicIrrep(4, 1, "two_deg_blocks", S, T, 1, K, 1)
        stratum = Stratum(Dict(1 => 1), 4)

        candidates = find_mtcs_at_prime([atom], stratum, 5;
                                        search_mode = :exhaustive,
                                        max_block_dim = 2,
                                        precheck_unit_axiom = false)
        @test candidates isa Vector{MTCCandidate}
    end

    @testset "is_orthogonal_mod_p" begin
        p = 13
        @test is_orthogonal_mod_p([1 0; 0 1], p)
        @test is_orthogonal_mod_p([0 p-1; 1 0], p)  # 90-degree rotation
        @test !is_orthogonal_mod_p([1 1; 0 1], p)
        @test !is_orthogonal_mod_p([1 0 0; 0 1 0], p)  # non-square
    end
end
