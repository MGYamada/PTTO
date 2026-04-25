using Test
using ACMG

@testset "Pipeline prime selection and modes" begin
    @testset "N_effective is fixed to N" begin
        @test ACMG.compute_effective_conductor(1) == 1
        @test ACMG.compute_effective_conductor(8) == 8
        @test ACMG.compute_effective_conductor(20) == 20
    end

    @testset "prime validity uses N directly" begin
        err = try
            ACMG.classify_mtcs_at_conductor(8;
                                            max_rank = 1,
                                            primes = [17, 19],
                                            skip_FR = true,
                                            verbose = false)
            nothing
        catch e
            e
        end

        @test err isa ErrorException
        msg = sprint(showerror, err)
        @test occursin("N_effective | p-1", msg)
        @test occursin("input N=8", msg)
        @test occursin("N_effective=8", msg)
    end

    @testset "classify_mtcs_at_conductor auto-selects primes when omitted" begin
        err = try
            ACMG.classify_mtcs_at_conductor(24;
                                            max_rank = 1,
                                            primes = nothing,
                                            min_primes = 2,
                                            prime_start = 29,
                                            prime_window = 10,
                                            skip_FR = true,
                                            verbose = false)
            nothing
        catch e
            e
        end

        @test err isa ErrorException
        msg = sprint(showerror, err)
        @test occursin("insufficient admissible primes", msg)
        @test occursin("N_eff=24", msg)
    end

    @testset "conductor_mode=:T_only is removed" begin
        err = try
            ACMG.classify_mtcs_at_conductor(8;
                                            max_rank = 1,
                                            primes = [17, 41],
                                            conductor_mode = :T_only,
                                            skip_FR = true,
                                            verbose = false)
            nothing
        catch e
            e
        end

        @test err isa ErrorException
        @test occursin("removed", sprint(showerror, err))
    end

    @testset "classify_mtcs_auto returns reproducibility metadata" begin
        auto = ACMG.classify_mtcs_auto(8;
                                       max_rank_candidates = [1],
                                       d_candidates = [1],
                                       min_primes = 2,
                                       prime_start = 11,
                                       prime_max = 50,
                                       max_attempts = 1,
                                       skip_FR = true,
                                       verbose = false)

        @test haskey(auto, :classified)
        @test haskey(auto, :N_effective)
        @test haskey(auto, :primes)
        @test haskey(auto, :max_rank)
        @test auto.N_input == 8
        @test auto.N_effective == 8
        @test auto.conductor_mode == :full_mtc
        @test auto.max_rank == 1
        @test length(auto.primes) == 2
        @test auto.attempts == 1
    end

    @testset "select_admissible_primes picks valid primes" begin
        ps = ACMG.select_admissible_primes(24;
                                           min_count = 3,
                                           start_from = 29,
                                           window = 200)
        @test length(ps) == 3
        @test all(p -> (p - 1) % 24 == 0, ps)
        @test issorted(ps)
    end

    @testset "select_admissible_primes reports searched range on shortage" begin
        err = try
            ACMG.select_admissible_primes(24;
                                          min_count = 2,
                                          start_from = 29,
                                          window = 10)
            nothing
        catch e
            e
        end

        @test err isa ErrorException
        msg = sprint(showerror, err)
        @test occursin("insufficient admissible primes", msg)
        @test occursin("(29, 39]", msg)
        @test occursin("found", msg)
    end
end
