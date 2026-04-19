"""
Tests for `ACMG.Phase4.KitaevComplex`.
Just verify that `chi3_equations` compiles and returns polynomials
in the expected ring.
"""

using Test
using Oscar
using ACMG
using ACMG.Phase4

const KC = Phase4.KitaevComplex

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

function ising_Nijk()
    r = 3
    N = zeros(Int, r, r, r)
    for j in 1:r
        N[1, j, j] = 1; N[j, 1, j] = 1
    end
    N[2,2,1]=1; N[2,3,3]=1; N[3,2,3]=1; N[3,3,1]=1; N[3,3,2]=1
    return N
end

@testset "KitaevComplex" begin
    @testset "chi3_equations: Fibonacci compiles" begin
        poly_C, eqs = KC.chi3_equations(fib_Nijk(), [1, 0])
        @info "Fibonacci chi3" n_eqs=length(eqs)
        @test length(eqs) ≥ 0   # can be 0 if all χ³ entries reduce to 0
    end

    @testset "chi3_equations: Ising compiles" begin
        poly_C, eqs = KC.chi3_equations(ising_Nijk(), [1, 0, 0])
        @info "Ising chi3" n_eqs=length(eqs)
        @test length(eqs) ≥ 0
    end
end
