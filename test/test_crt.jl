using Test
using ACMG

"""
Tests for CRT module (Phase 3) — pure F_p parts, no Oscar.
"""

@testset "CRT primitives" begin

    @testset "crt2" begin
        # x ≡ 1 (mod 3), x ≡ 2 (mod 5) → x = 7 (mod 15)
        (x, M) = crt2(1, 3, 2, 5)
        @test x == 7
        @test M == 15

        # x ≡ 5 (mod 7), x ≡ 3 (mod 11) → ?
        (x, M) = crt2(5, 7, 3, 11)
        @test M == 77
        @test mod(x, 7) == 5
        @test mod(x, 11) == 3
    end

    @testset "acmg_crt (multi-modulus)" begin
        residues = [2, 3, 2]
        moduli = [3, 5, 7]
        # x ≡ 2 (mod 3), x ≡ 3 (mod 5), x ≡ 2 (mod 7)
        # classic: x = 23 (mod 105)
        (x, M) = acmg_crt(residues, moduli)
        @test M == 105
        @test mod(x, 3) == 2
        @test mod(x, 5) == 3
        @test mod(x, 7) == 2
        @test x == 23

        # Single residue
        (x, M) = acmg_crt([5], [7])
        @test x == 5
        @test M == 7
    end

    @testset "rational_reconstruct" begin
        # Simple: 3 mod 100 → 3/1
        result = rational_reconstruct(3, 100)
        @test result !== nothing
        @test result == (BigInt(3), BigInt(1))

        # 1/3 mod 100: 1/3 ≡ 67 (since 3·67 = 201 ≡ 1)
        # Reconstruct 67 → (1, 3)
        result = rational_reconstruct(67, 100)
        @test result !== nothing
        (a, b) = result
        # Should satisfy a/b ≡ 67 (mod 100) with small a, b
        @test mod(a * invmod(b, 100), 100) == 67
        @test abs(a) <= 10 && abs(b) <= 10

        # -1 mod 100 = 99. Reconstruct 99 → (-1, 1)
        result = rational_reconstruct(99, 100)
        @test result !== nothing
        (a, b) = result
        @test mod(a - 99 * b, 100) == 0  # a ≡ 99b (mod 100)
    end

    @testset "compute_sqrt_d_mod_p" begin
        @test compute_sqrt_d_mod_p(3, 73) in (21, 52)  # 21² = 441 ≡ 3 (441 = 6·73+3)
        @test mod(compute_sqrt_d_mod_p(3, 73)^2, 73) == 3
        @test mod(compute_sqrt_d_mod_p(2, 17)^2, 17) == 2
        @test compute_sqrt_d_mod_p(3, 5) === nothing  # 3 is not QR mod 5
    end

    @testset "reconstruct_in_Z_sqrt_d (√3)" begin
        # Known: -√3 mod 73 = -21 = 52; -√3 mod 97 = ?
        # sqrt(3) mod 97: let's check
        s3_73 = compute_sqrt_d_mod_p(3, 73)
        s3_97 = compute_sqrt_d_mod_p(3, 97)
        s3_193 = compute_sqrt_d_mod_p(3, 193)

        @test s3_73 !== nothing
        @test s3_97 !== nothing
        @test s3_193 !== nothing

        # Value x = 0 + (-1)·√3 = -√3
        # Construct values by prime
        vals = Dict(73 => mod(-s3_73, 73),
                    97 => mod(-s3_97, 97),
                    193 => mod(-s3_193, 193))
        sqrtd = Dict(73 => s3_73, 97 => s3_97, 193 => s3_193)
        result = reconstruct_in_Z_sqrt_d(vals, 3, sqrtd; bound = 5)
        @test result !== nothing
        @test result == (BigInt(0), BigInt(-1))

        # Value x = 2 + 3√3
        vals2 = Dict(p => mod(2 + 3 * sqrtd[p], p) for p in [73, 97, 193])
        result2 = reconstruct_in_Z_sqrt_d(vals2, 3, sqrtd; bound = 5)
        @test result2 !== nothing
        @test result2 == (BigInt(2), BigInt(3))

        # Non-matching value (not in Z[√3] with small coeffs) — should fail
        vals3 = Dict(73 => 7, 97 => 11, 193 => 13)  # random values
        result3 = reconstruct_in_Z_sqrt_d(vals3, 3, sqrtd; bound = 3)
        @test result3 === nothing
    end

    @testset "reconstruct_matrix_in_Z_sqrt_d" begin
        # Small matrix in Z[√3]: M = [[1, -√3], [√3, 2]]
        # Reduce to primes 73, 97, 193
        test_primes = [73, 97, 193]
        sqrtd = Dict{Int, Int}()
        for p in test_primes
            sqrtd[p] = compute_sqrt_d_mod_p(3, p)
        end

        # Build F_p matrices
        M_true = [(1, 0) (0, -1); (0, 1) (2, 0)]
        matrix_by_prime = Dict{Int, Matrix{Int}}()
        for p in test_primes
            s = sqrtd[p]
            M_p = [mod(a + b*s, p) for (a, b) in M_true]
            matrix_by_prime[p] = M_p
        end

        recon = reconstruct_matrix_in_Z_sqrt_d(matrix_by_prime, 3; bound = 5)
        @test recon == M_true
    end

    @testset "group_mtcs_by_fusion" begin
        # Build two dummy candidates with same N tensor at two primes
        # and a third different one
        N1 = zeros(Int, 2, 2, 2)
        N1[1, 1, 1] = 1; N1[1, 2, 2] = 1; N1[2, 1, 2] = 1; N1[2, 2, 1] = 1; N1[2, 2, 2] = 1

        N2 = zeros(Int, 2, 2, 2)  # different
        N2[1, 1, 1] = 1; N2[2, 2, 1] = 1  # only nonzero

        c1_p1 = MTCCandidate(41, :dummy, [14 16; 16 27], [1, 18], 1, N1, [1, 7], 9)
        c1_p2 = MTCCandidate(61, :dummy, [1 2; 2 3], [1, 2], 1, N1, [1, 5], 8)
        c2_p1 = MTCCandidate(41, :dummy, [1 0; 0 1], [1, 1], 1, N2, [1, 1], 2)

        results = Dict(41 => [c1_p1, c2_p1], 61 => [c1_p2])
        groups = group_mtcs_by_fusion(results)

        # Expect 2 groups: {41→c1, 61→c1'} and {41→c2}
        @test length(groups) == 2
        # Find the group containing N1
        g_N1 = nothing
        for g in groups
            rep = first(values(g))
            if rep.N == N1
                g_N1 = g
                break
            end
        end
        @test g_N1 !== nothing
        @test length(g_N1) == 2
        @test haskey(g_N1, 41) && haskey(g_N1, 61)
    end

    @testset "build_sqrtd_selector + galois grouping branch alignment" begin
        # d=3 uses cyclotomic selector
        cyclo = ACMG.build_sqrtd_selector(3, [73], 73; verbose = false)
        @test cyclo.mode == :cyclotomic
        s3 = cyclo.sqrtd_fn(3, 73)
        @test mod(s3 * s3, 73) == 3

        # d=6 uses anchored selector with per-prime sign cache
        d = 6
        p_anchor = 29
        p_other = 53
        sel = ACMG.build_sqrtd_selector(d, [p_anchor, p_other], p_anchor; verbose = false)
        @test sel.mode == :anchored

        # Build synthetic rank-1 candidate data where p_other requires
        # the opposite sqrt branch to be reconstruction-compatible.
        N1 = zeros(Int, 1, 1, 1)
        N1[1, 1, 1] = 1
        true_x = (a = 5, b = 1)  # x = a + b*√d
        s_anchor = ACMG.compute_sqrt_d_mod_p(d, p_anchor)
        s_other_raw = ACMG.compute_sqrt_d_mod_p(d, p_other)
        s_other_true = mod(-s_other_raw, p_other)  # force opposite branch

        S_anchor = mod((true_x.a + true_x.b * s_anchor) * invmod(2 * s_anchor, p_anchor), p_anchor)
        S_other = mod((true_x.a + true_x.b * s_other_true) * invmod(2 * s_other_true, p_other), p_other)

        c_anchor = MTCCandidate(p_anchor, :dummy, reshape([S_anchor], 1, 1),
                                [1], 1, N1, [1], 1)
        c_other = MTCCandidate(p_other, :dummy, reshape([S_other], 1, 1),
                               [1], 1, N1, [1], 1)
        results = Dict(p_anchor => [c_anchor], p_other => [c_other])

        groups = ACMG.group_mtcs_galois_aware(results, p_anchor;
                                              scale_d = d,
                                              sqrtd_fn = sel.sqrtd_fn,
                                              branch_sign_getter = sel.branch_sign_getter,
                                              branch_sign_setter = sel.branch_sign_setter)

        @test length(groups) == 1
        @test haskey(groups[1], p_anchor)
        @test haskey(groups[1], p_other)
        @test sel.branch_sign_getter(p_other) in (-1, 1)
    end

    @testset "describe_matrix" begin
        # Small 2x2 to test formatting
        M = [(1, 0) (0, -1); (0, 1) (2, 0)]
        s = describe_matrix(M, 3)
        @test occursin("1", s)
        @test occursin("√3", s)
        @test occursin("-", s)  # negative
        @test occursin("2", s)
    end
end
