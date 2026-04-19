"""
Tests for `ACMG.Phase4.SlicedPentagonSolver`.

Categorical χ³ F = 0 slice tests:
  1. Shape check: system builds, pentagon + slice eqs coexist in one ring.
  2. χ³ F_base ≈ 0 check: if the categorical convention is right, the
     known Fibonacci/Ising F-symbols should satisfy χ³ F = 0 to
     numerical precision.
  3. Newton from base F: must converge instantly.
"""

using Test
using LinearAlgebra
using Oscar
using ACMG
using ACMG.Phase4

const KC  = Phase4.KitaevComplex
const SPS = Phase4.SlicedPentagonSolver
const PS  = Phase4.PentagonSolver

function fib_Nijk()
    r = 2
    N = zeros(Int, r, r, r)
    N[1,1,1]=1; N[1,2,2]=1; N[2,1,2]=1; N[2,2,1]=1; N[2,2,2]=1
    return N
end

function fib_F_func()
    φ = (1 + sqrt(5)) / 2
    sφ = sqrt(φ)
    F_τττ_τ = [1/φ    1/sφ;
               1/sφ  -1/φ]
    fr = FusionRule(fib_Nijk())
    function F(a,b,c,d,e,f)
        (fr.N[a,b,e]==1 && fr.N[e,c,d]==1 && fr.N[b,c,f]==1 && fr.N[a,f,d]==1) || return 0.0
        1 in (a,b,c,d) && return 1.0
        return F_τττ_τ[e, f]
    end
    return F
end

function ising_Nijk()
    r = 3
    N = zeros(Int, r, r, r)
    for j in 1:r
        N[1, j, j] = 1; N[j, 1, j] = 1
    end
    N[2,2,1]=1; N[2,3,3]=1; N[3,2,3]=1; N[3,3,1]=1; N[3,3,2]=1
    return N
end

function ising_F_func()
    F_sss_s = (1.0/sqrt(2.0)) * [1.0 1.0; 1.0 -1.0]
    fr = FusionRule(ising_Nijk())
    function F(a,b,c,d,e,f)
        (fr.N[a,b,e]==1 && fr.N[e,c,d]==1 && fr.N[b,c,f]==1 && fr.N[a,f,d]==1) || return 0.0
        1 in (a,b,c,d) && return 1.0
        nsig = count(==(3), (a,b,c,d))
        nsig == 0 && return 1.0
        if nsig == 2
            (a,b,c,d) == (2,3,2,3) && return -1.0
            (a,b,c,d) == (3,2,3,1) && return  1.0
            (a,b,c,d) == (3,2,3,2) && return -1.0
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

# Build base F vector ordered by the same pop! traversal as
# `assign_F_to_associator!`.
function base_F_vector(Nijk, r, F_func)
    one_vec = zeros(Int, r); one_vec[1] = 1
    flat_keys = NTuple{6,Int}[]
    for i in 1:r, j in 1:r, k in 1:r, o in 1:r
        sum(one_vec[[i,j,k]]) > 0 && continue
        rows_e = Int[]; for e in 1:r
            if Nijk[i,j,e] ≥ 1 && Nijk[e,k,o] ≥ 1
                push!(rows_e, e)
            end
        end
        cols_f = Int[]; for f in 1:r
            if Nijk[j,k,f] ≥ 1 && Nijk[i,f,o] ≥ 1
                push!(cols_f, f)
            end
        end
        (isempty(rows_e) || isempty(cols_f)) && continue
        for e in rows_e, f in cols_f
            push!(flat_keys, (i,j,k,o,e,f))
        end
    end
    n = length(flat_keys)
    v = zeros(ComplexF64, n)
    for (pos, key) in enumerate(flat_keys)
        v[n - pos + 1] = ComplexF64(F_func(key...))
    end
    return v
end

@testset "SlicedPentagonSolver" begin

    @testset "get_sliced_pentagon_system: shapes" begin
        R_fib, pent_fib, slice_fib, n_fib = SPS.get_sliced_pentagon_system(fib_Nijk(), 2)
        @test n_fib == 5
        @test length(pent_fib) > 0
        @info "Fibonacci sliced system" n=n_fib n_pent=length(pent_fib) n_slice=length(slice_fib)

        R_is, pent_is, slice_is, n_is = SPS.get_sliced_pentagon_system(ising_Nijk(), 3)
        @test n_is == 14
        @test length(pent_is) > 0
        @info "Ising sliced system" n=n_is n_pent=length(pent_is) n_slice=length(slice_is)
    end

    @testset "χ³ F_base residual" begin
        # Fibonacci
        F_fib = base_F_vector(fib_Nijk(), 2, fib_F_func())
        _R, _pent, slice_fib, _n = SPS.get_sliced_pentagon_system(fib_Nijk(), 2)
        res_fib = [PS.eval_poly_complex(p, F_fib) for p in slice_fib]
        max_fib = isempty(res_fib) ? 0.0 : maximum(abs.(res_fib))
        @info "Fibonacci χ³ F_base residual" n_slice=length(slice_fib) max_resid=max_fib

        # Ising
        F_is = base_F_vector(ising_Nijk(), 3, ising_F_func())
        _R_is, _pent_is, slice_is, _n_is = SPS.get_sliced_pentagon_system(ising_Nijk(), 3)
        res_is = [PS.eval_poly_complex(p, F_is) for p in slice_is]
        max_is = isempty(res_is) ? 0.0 : maximum(abs.(res_is))
        @info "Ising χ³ F_base residual" n_slice=length(slice_is) max_resid=max_is

        # Hard-fail if the residual is too large
        @test max_fib < 1e-6
        @test max_is  < 1e-6
    end

    @testset "Newton from base F (Fibonacci, Ising)" begin
        sols_fib = SPS.solve_pentagon_newton_with_slice(fib_Nijk(), 2, fib_F_func();
            max_trials = 1, max_iter = 50, perturb_scale = 0.0, tol = 1e-12)
        @test length(sols_fib) == 1

        sols_is = SPS.solve_pentagon_newton_with_slice(ising_Nijk(), 3, ising_F_func();
            max_trials = 1, max_iter = 50, perturb_scale = 0.0, tol = 1e-12)
        @test length(sols_is) == 1
    end
end
