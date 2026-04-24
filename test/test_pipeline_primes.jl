using Test
using ACMG

@testset "Pipeline prime selection and modes" begin
    @testset "conductor_mode default uses auto N_effective in prime validity check" begin
        err = try
            ACMG.classify_mtcs_at_conductor(1;
                                            max_rank = 1,
                                            primes = [5, 13],
                                            scale_d = 2,
                                            skip_FR = true,
                                            verbose = false)
            nothing
        catch e
            e
        end

        @test err isa ErrorException
        msg = sprint(showerror, err)
        @test occursin("N_effective | p-1", msg)
        @test occursin("input N=1", msg)
        @test occursin("N_effective=24", msg)
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

    @testset "conductor_mode=:T_only is removed in v0.5.0" begin
        err = try
            ACMG.classify_mtcs_at_conductor(1;
                                            max_rank = 1,
                                            primes = [73, 97],
                                            scale_d = 2,
                                            conductor_mode = :T_only,
                                            skip_FR = true,
                                            verbose = false)
            nothing
        catch e
            e
        end

        @test err isa ErrorException
        @test occursin("removed in v0.5.0", sprint(showerror, err))
    end

    @testset "classify_mtcs_auto returns reproducibility metadata" begin
        auto = ACMG.classify_mtcs_auto(1;
                                       max_rank_candidates = [1],
                                       scale_d_candidates = [2],
                                       d_candidates = [2],
                                       min_primes = 2,
                                       prime_start = 29,
                                       prime_max = 200,
                                       max_attempts = 1,
                                       skip_FR = true,
                                       verbose = false)

        @test haskey(auto, :classified)
        @test haskey(auto, :N_effective)
        @test haskey(auto, :scale_d)
        @test haskey(auto, :primes)
        @test haskey(auto, :max_rank)
        @test auto.N_input == 1
        @test auto.N_effective == 24
        @test auto.scale_d == 2
        @test auto.conductor_mode == :full_mtc
        @test auto.max_rank == 1
        @test length(auto.primes) == 2
        @test auto.attempts == 1
    end

    @testset "classify_mtcs_auto rejects :T_only in v0.5.0" begin
        err = try
            ACMG.classify_mtcs_auto(1;
                                    max_rank_candidates = [1],
                                    scale_d_candidates = [2],
                                    d_candidates = [1],
                                    conductor_modes = [:T_only],
                                    min_primes = 2,
                                    prime_start = 29,
                                    prime_max = 39,
                                    max_attempts = 1,
                                    skip_FR = true,
                                    verbose = false)
            nothing
        catch e
            e
        end

        @test err isa ErrorException
        @test occursin("removed in v0.5.0", sprint(showerror, err))
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

    @testset "classify_mtcs_auto records prime-search shortage reason" begin
        auto = ACMG.classify_mtcs_auto(1;
                                       max_rank_candidates = [1],
                                       scale_d_candidates = [2],
                                       d_candidates = [1],
                                       min_primes = 2,
                                       prime_start = 29,
                                       prime_max = 39,
                                       max_attempts = 1,
                                       skip_FR = true,
                                       verbose = false)

        @test isempty(auto.classified)
        @test length(auto.history) == 1
        @test auto.history[1].executed
        @test !auto.history[1].success
        @test occursin("insufficient admissible primes", auto.history[1].reason)
        @test occursin("(29, 39]", auto.history[1].reason)
    end
end
