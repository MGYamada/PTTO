using Test
using ACMG

@testset "Pipeline galois grouping" begin
    @testset "_branch_consistency_precheck resolves anchored sign contradiction" begin
        d = 6
        p_anchor = 29
        p_other = 53
        selector = ACMG.build_sqrtd_selector(d, [p_anchor, p_other], p_anchor; verbose = false)
        selector.branch_sign_setter(p_other, -1)
        orig_sign = selector.branch_sign_getter(p_other)

        N1 = zeros(Int, 1, 1, 1)
        N1[1, 1, 1] = 1
        x = (a = 3, b = 1)
        s_anchor = ACMG.compute_sqrt_d_mod_p(d, p_anchor)
        s_other_raw = ACMG.compute_sqrt_d_mod_p(d, p_other)
        s_other_true = mod(-s_other_raw, p_other)

        S_anchor = mod((x.a + x.b * s_anchor) * invmod(2 * s_anchor, p_anchor), p_anchor)
        S_other = mod((x.a + x.b * s_other_true) * invmod(2 * s_other_true, p_other), p_other)

        c_anchor = ACMG.MTCCandidate(p_anchor, :dummy, reshape([S_anchor], 1, 1),
                                     [1], 1, N1, [1], 1)
        c_other = ACMG.MTCCandidate(p_other, :dummy, reshape([S_other], 1, 1),
                                    [1], 1, N1, [1], 1)
        results = Dict(p_anchor => [c_anchor], p_other => [c_other])

        contradictions = ACMG._branch_consistency_precheck(results, p_anchor, d, selector.sqrtd_fn;
                                                           branch_sign_getter = selector.branch_sign_getter,
                                                           branch_sign_setter = selector.branch_sign_setter,
                                                           verbose = false)
        @test isempty(contradictions)
        @test selector.branch_sign_getter(p_other) == orig_sign
    end

    @testset "_branch_consistency_precheck explores all anchor candidates (rank-2, d=5, p=41/61)" begin
        d = 5
        p_anchor = 41
        p_other = 61

        N2 = zeros(Int, 2, 2, 2)
        for i in 1:2
            N2[i, i, i] = 1
        end

        s_anchor = ACMG.compute_sqrt_d_mod_p(d, p_anchor)
        s_other = ACMG.compute_sqrt_d_mod_p(d, p_other)
        two_s_anchor_inv = invmod(mod(2 * s_anchor, p_anchor), p_anchor)
        two_s_other_inv = invmod(mod(2 * s_other, p_other), p_other)

        # Entries are encoded as x = a + b*√d so that 2*√d*S = x mod p.
        good_x = [
            (3, 1) (0, 1)
            (0, 1) (2, -1)
        ]
        bad_x = [
            (11, 4) (7, -3)
            (5, 2) (9, 1)
        ]

        encode_entry(pair, s, p, two_s_inv) = mod((pair[1] + pair[2] * s) * two_s_inv, p)
        good_anchor_S = [encode_entry(good_x[i, j], s_anchor, p_anchor, two_s_anchor_inv) for i in 1:2, j in 1:2]
        good_other_S = [encode_entry(good_x[i, j], s_other, p_other, two_s_other_inv) for i in 1:2, j in 1:2]
        bad_anchor_S = [encode_entry(bad_x[i, j], s_anchor, p_anchor, two_s_anchor_inv) for i in 1:2, j in 1:2]

        # Put an incompatible anchor candidate first to reproduce the regression.
        c_anchor_bad = ACMG.MTCCandidate(p_anchor, :dummy_bad, bad_anchor_S,
                                         [1, 1], 1, N2, [1, 1], 2)
        c_anchor_good = ACMG.MTCCandidate(p_anchor, :dummy_good, good_anchor_S,
                                          [1, 1], 1, N2, [1, 1], 2)
        c_other = ACMG.MTCCandidate(p_other, :dummy_other, good_other_S,
                                    [1, 1], 1, N2, [1, 1], 2)

        results = Dict(p_anchor => [c_anchor_bad, c_anchor_good], p_other => [c_other])
        contradictions = ACMG._branch_consistency_precheck(results, p_anchor, d, ACMG.compute_sqrt_d_mod_p;
                                                           verbose = false)
        @test isempty(contradictions)
    end

    @testset "cyclotomic-signed selector flips branch in precheck/grouping (d=5)" begin
        d = 5
        p_anchor = 41
        p_other = 61
        selector = ACMG.build_sqrtd_selector(d, [p_anchor, p_other], p_anchor; verbose = false)
        @test selector.mode in (:cyclotomic, :cyclotomic_signed)

        N1 = zeros(Int, 1, 1, 1)
        N1[1, 1, 1] = 1
        x = (a = 4, b = 1)
        s_anchor = selector.sqrtd_fn(d, p_anchor)
        s_other_true = mod(-selector.sqrtd_fn(d, p_other), p_other) # force opposite branch

        S_anchor = mod((x.a + x.b * s_anchor) * invmod(mod(2 * s_anchor, p_anchor), p_anchor), p_anchor)
        S_other = mod((x.a + x.b * s_other_true) * invmod(mod(2 * s_other_true, p_other), p_other), p_other)

        c_anchor = ACMG.MTCCandidate(p_anchor, :dummy, reshape([S_anchor], 1, 1),
                                     [1], 1, N1, [1], 1)
        c_other = ACMG.MTCCandidate(p_other, :dummy, reshape([S_other], 1, 1),
                                    [1], 1, N1, [1], 1)
        results = Dict(p_anchor => [c_anchor], p_other => [c_other])

        selector.branch_sign_setter(p_other, 1)
        sign_trials = Int[]
        wrapped_setter = (p, sgn) -> begin
            push!(sign_trials, sgn >= 0 ? 1 : -1)
            selector.branch_sign_setter(p, sgn)
        end

        contradictions = ACMG._branch_consistency_precheck(results, p_anchor, d, selector.sqrtd_fn;
                                                           branch_sign_getter = selector.branch_sign_getter,
                                                           branch_sign_setter = wrapped_setter,
                                                           verbose = false)
        @test isempty(contradictions)
        @test 1 in sign_trials && -1 in sign_trials
        @test selector.branch_sign_getter(p_other) == 1

        groups = ACMG.group_mtcs_galois_aware(results, p_anchor;
                                              scale_d = d,
                                              sqrtd_fn = selector.sqrtd_fn,
                                              branch_sign_getter = selector.branch_sign_getter,
                                              branch_sign_setter = selector.branch_sign_setter)
        @test length(groups) == 1
        @test haskey(groups[1], p_other)
        @test selector.branch_sign_getter(p_other) == -1
    end

    @testset "regression: N=5, scale_d=5, primes=[41,61], rank-2 is not dropped as single-prime sector" begin
        classified = ACMG.classify_mtcs_at_conductor(5;
                                                     max_rank = 2,
                                                     primes = [41, 61],
                                                     scale_d = 5,
                                                     skip_FR = true,
                                                     verbose = false)

        @test any(c -> c.rank == 2 && length(c.used_primes) >= 2, classified)
    end
end
