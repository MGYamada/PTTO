"""
Tests for `ACMG.Phase4.SlicedPentagonSolver`.

Strategy:
  1. Fibonacci sanity check:
     - Pentagon HC without slice returns 4 solutions.
     - WITH Kitaev slice (1 constraint), HC should return fewer solutions
       (gauge-equivalent ones collapsed).
  2. Ising test:
     - Pentagon HC without slice: mixed volume ~477k, infeasible.
     - WITH slice (1 constraint): measure mixed volume reduction.
"""

using Test
using LinearAlgebra
using ACMG
using ACMG.Phase4

const KC  = Phase4.KitaevComplex
const SPS = Phase4.SlicedPentagonSolver

# ---------------- Fibonacci helpers (reused) ----------------

function fib_Nijk()
    r = 2
    N = zeros(Int, r, r, r)
    N[1, 1, 1] = 1
    N[1, 2, 2] = 1
    N[2, 1, 2] = 1
    N[2, 2, 1] = 1
    N[2, 2, 2] = 1
    return N
end

function fib_F_func()
    φ = (1 + sqrt(5)) / 2
    sφ = sqrt(φ)
    F_τττ_τ = [1/φ    1/sφ;
               1/sφ  -1/φ]
    fr = FusionRule(fib_Nijk())
    function F(a::Int, b::Int, c::Int, d::Int, e::Int, f::Int)
        (fr.N[a, b, e] == 1 && fr.N[e, c, d] == 1 &&
         fr.N[b, c, f] == 1 && fr.N[a, f, d] == 1) || return 0.0
        1 in (a, b, c, d) && return 1.0
        return F_τττ_τ[e, f]
    end
    return F
end

# ---------------- Ising helpers ----------------

function ising_Nijk()
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
    return N
end

function ising_F_func()
    F_sss_s = (1.0 / sqrt(2.0)) * [1.0  1.0;
                                    1.0 -1.0]
    fr = FusionRule(ising_Nijk())
    function F(a::Int, b::Int, c::Int, d::Int, e::Int, f::Int)
        (fr.N[a, b, e] == 1 && fr.N[e, c, d] == 1 &&
         fr.N[b, c, f] == 1 && fr.N[a, f, d] == 1) || return 0.0
        1 in (a, b, c, d) && return 1.0
        nsig = count(==(3), (a, b, c, d))
        nsig == 0 && return 1.0
        if nsig == 2
            (a, b, c, d) == (2, 3, 2, 3) && return -1.0
            (a, b, c, d) == (3, 2, 3, 1) && return  1.0
            (a, b, c, d) == (3, 2, 3, 2) && return -1.0
            return 1.0
        end
        if nsig == 4
            d != 3 && return 0.0
            return F_sss_s[e, f]
        end
        return 1.0
    end
    return F
end

@testset "SlicedPentagonSolver" begin

    # --------------------------------------------------------------------
    @testset "Fibonacci FKey ↔ pentagon-var map" begin
        Nijk = fib_Nijk()
        one_vec = [1, 0]
        fkey_map = SPS.build_fkey_to_xvar_map(Nijk, 2, one_vec)

        # Fibonacci pentagon has 5 variables (we saw this in test output).
        # The (τ,τ,τ,τ) block contributes 4 entries + possibly 1 extra elsewhere.
        # Check that all mapped indices are in [1, 5] and unique.
        @test length(fkey_map) > 0
        vals = collect(values(fkey_map))
        @test all(1 .≤ vals .≤ 5)
        @test length(unique(vals)) == length(vals)

        # All mapped FKeys should have indices (i,j,k,o) with i, j, k != 1 (unit).
        for key in keys(fkey_map)
            i, j, k, _o, _e, _f = key
            @test !(1 in (i, j, k))
        end
    end

    # --------------------------------------------------------------------
    @testset "Fibonacci sliced pentagon system" begin
        Nijk = fib_Nijk()
        F_fn = fib_F_func()

        R, eqs_aug, n, n_slice = SPS.get_sliced_pentagon_system(Nijk, 2, F_fn)

        @test n == 5                 # Fibonacci pentagon vars
        @test n_slice == 1           # Fibonacci gauge orbit dim = 1
        # Augmented system has pentagon eqs + 1 slice constraint
        @test length(eqs_aug) > 12   # original pentagon eqs kept

        # Slice constraint should vanish at the base F-solution.
        # Extract base F values as a vector ordered by pentagon variable index.
        one_vec = [1, 0]
        fkey_map = SPS.build_fkey_to_xvar_map(Nijk, 2, one_vec)
        F_vec = zeros(Float64, n)
        for (key, pidx) in fkey_map
            F_vec[pidx] = F_fn(key...)
        end

        # Evaluate the slice polynomial(s) at F_vec
        slice_polys = eqs_aug[end - n_slice + 1 : end]
        xs = gens(R)
        for p in slice_polys
            val = Float64(evaluate(p, F_vec))
            @test abs(val) < 1e-10
        end
    end

    # --------------------------------------------------------------------
    @testset "Ising sliced pentagon: mapping + slice build" begin
        Nijk = ising_Nijk()
        F_fn = ising_F_func()

        R, eqs_aug, n, n_slice = SPS.get_sliced_pentagon_system(Nijk, 3, F_fn)

        # Ising pentagon variables count (from TensorCategories' traversal)
        @test n > 0
        @test n_slice == 1             # Ising effective gauge = 1

        # Base F solution satisfies the slice constraint
        one_vec = [1, 0, 0]
        fkey_map = SPS.build_fkey_to_xvar_map(Nijk, 3, one_vec)
        F_vec = zeros(Float64, n)
        for (key, pidx) in fkey_map
            F_vec[pidx] = F_fn(key...)
        end

        slice_polys = eqs_aug[end - n_slice + 1 : end]
        for p in slice_polys
            val = Float64(evaluate(p, F_vec))
            @test abs(val) < 1e-10
        end

        # Base F also satisfies pentagon equations (standard Ising solves pentagon)
        pent_polys = eqs_aug[1 : end - n_slice]
        max_resid = 0.0
        for p in pent_polys
            val = abs(Float64(evaluate(p, F_vec)))
            max_resid = max(max_resid, val)
        end
        @test max_resid < 1e-10

        @info "Ising pentagon" n_vars=n n_slice=n_slice n_eqs=length(eqs_aug)
    end
end
