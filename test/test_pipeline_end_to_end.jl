using Test
using ACMG

@testset "Pipeline end-to-end" begin
    @testset "classify_mtcs_at_conductor rank-1 smoke (skip_FR)" begin
        classified = ACMG.classify_mtcs_at_conductor(1;
                                                     max_rank = 1,
                                                     primes = [73, 97],
                                                     skip_FR = true,
                                                     verbose = false)
        @test classified isa Vector{ACMG.ClassifiedMTC}
    end

    @testset "classify_mtcs_at_conductor rank-1 auto-prime exact path" begin
        classified = ACMG.classify_mtcs_at_conductor(1;
                                                     max_rank = 1,
                                                     min_primes = 2,
                                                     verbose = false)
        @test length(classified) == 1
        @test classified[1].verify_fresh
        @test classified[1].verify_report !== nothing
        @test classified[1].verify_report.ok
    end

    @testset "N=2 fresh primes validate with cyclotomic CRT" begin
        classified = ACMG.classify_mtcs_at_conductor(2;
                                                     max_rank = 5,
                                                     skip_FR = true,
                                                     verbose = false)
        @test length(classified) == 2
        @test all(m -> m.verify_fresh, classified)
        @test all(m -> m.verify_exact_lift === true, classified)
        @test all(m -> !isempty(m.fresh_primes), classified)
    end

    @testset "N=7 fixed atomic strata lift exactly" begin
        classified = ACMG.classify_mtcs_at_conductor(7;
                                                     max_rank = 3,
                                                     primes = [29, 43],
                                                     skip_FR = true,
                                                     verbose = false)
        @test length(classified) == 3
        @test any(m -> m.rank == 3, classified)
        @test all(m -> m.verify_fresh, classified)
        @test all(m -> m.verify_exact_lift === true, classified)
    end

    @testset "N=8 rejects Ising-like rank-3 sectors with invalid twists" begin
        classified = ACMG.classify_mtcs_at_conductor(8;
                                                     max_rank = 3,
                                                     primes = [41, 73],
                                                     verbose = false)
        @test length(classified) == 4
        @test count(m -> m.rank == 3, classified) == 0
        @test all(m -> m.verify_report !== nothing && m.verify_report.ok,
                  filter(m -> m.rank <= 2, classified))
    end

    @testset "N=8 auto primes keep exact fixed-stratum sectors" begin
        classified = ACMG.classify_mtcs_at_conductor(8;
                                                     max_rank = 3,
                                                     verbose = false)
        @test length(classified) == 4
        @test count(m -> m.rank == 2, classified) == 2
        @test count(m -> m.rank == 3, classified) == 0
        @test all(m -> m.verify_fresh, classified)
        @test all(m -> m.verify_exact_lift === true, classified)
        @test all(m -> m.verify_report !== nothing && m.verify_report.ok,
                  filter(m -> m.rank <= 2, classified))
    end

    @testset "N=20 exact layer is single-pass and roundtrips" begin
        classified = ACMG.classify_mtcs_at_conductor(20;
                                                     max_rank = 2,
                                                     primes = [41, 61],
                                                     skip_FR = false,
                                                     verbose = false)
        @test length(classified) == 6
        @test all(m -> m.verify_report !== nothing && m.verify_report.ok, classified)
        @test all(m -> iszero(m.verify_report.S_max), classified)
        @test all(m -> iszero(m.verify_report.T_max), classified)
        @test !occursin("pent=? hex=?", sprint(show, classified))
    end

    @testset "fallback cyclotomic CRT caps oversized pipeline bound" begin
        p1, p2 = 73, 97
        Nijk = ones(Int, 1, 1, 1)
        group = Dict(
            p1 => ACMG.MTCCandidate(p1, :searched, [1;;], [1], 1, Nijk, [1], 1),
            p2 => ACMG.MTCCandidate(p2, :searched, [1;;], [1], 1, Nijk, [1], 1),
        )

        classified = ACMG.classify_from_group(group, 24, ACMG.Stratum(Dict(1 => 1), 1),
                                              [p1, p2];
                                              reconstruction_bound = 50,
                                              skip_FR = true,
                                              verbose = false)

        @test classified.verify_fresh
        @test classified.verify_exact_lift === nothing
        @test classified.S_cyclotomic[1, 1] == one(parent(classified.S_cyclotomic[1, 1]))
        @test !haskey(ACMG._CRT_MITM_LEFT_CACHE,
                      (4, 50, (1, 1), (p1, p2)))
    end

    @testset "fusion grouping preserves same-rule sign variants" begin
        Nijk = zeros(Int, 2, 2, 2)
        Nijk[1, 1, 1] = 1
        Nijk[1, 2, 2] = 1
        Nijk[2, 1, 2] = 1
        Nijk[2, 2, 1] = 1
        p1, p2 = 17, 41
        cands = Dict(
            p1 => [
                ACMG.MTCCandidate(p1, :a, [1 1; 1 -1], [1, 4], 1, Nijk, [1, 1], 2),
                ACMG.MTCCandidate(p1, :b, [1 -1; -1 -1], [1, 4], 1, Nijk, [1, 1], 2),
            ],
            p2 => [
                ACMG.MTCCandidate(p2, :c, [1 1; 1 -1], [1, 9], 1, Nijk, [1, 1], 2),
                ACMG.MTCCandidate(p2, :d, [1 -1; -1 -1], [1, 9], 1, Nijk, [1, 1], 2),
            ],
        )

        groups = group_mtcs_by_fusion(cands)
        @test length(groups) == 4
        @test all(g -> sort(collect(keys(g))) == [p1, p2], groups)
        @test length(unique([(Tuple(vec(g[p1].S_Fp)), Tuple(vec(g[p2].S_Fp))) for g in groups])) == 4
    end
end
