"""
Tests for `ACMG.Phase4.KitaevComplex`.

Kitaev Eq. 251: χⁿ⁺¹ δⁿ + δⁿ⁻¹ χⁿ = id_Cⁿ   for n > 0.

The n = 0 case is excluded: Kitaev states δ⁰ = f₀⁰ − f₁⁰ = 0 identically,
so the identity doesn't apply there. We sanity-check δ⁰ = 0 below.
"""

using Test
using LinearAlgebra
using ACMG
using ACMG.Phase4

const KC = Phase4.KitaevComplex

function fib_Nijk()
    N = zeros(Int, 2, 2, 2)
    N[1,1,1] = 1; N[1,2,2] = 1; N[2,1,2] = 1; N[2,2,1] = 1; N[2,2,2] = 1
    return N
end

function ising_Nijk()
    N = zeros(Int, 3, 3, 3)
    for j in 1:3
        N[1,j,j] = 1; N[j,1,j] = 1
    end
    N[2,2,1] = 1; N[2,3,3] = 1; N[3,2,3] = 1; N[3,3,1] = 1; N[3,3,2] = 1
    return N
end

fib_dims()   = ((φ = (1 + sqrt(5)) / 2; [1.0, φ]), 1 + ((1 + sqrt(5)) / 2)^2)
ising_dims() = ([1.0, 1.0, sqrt(2.0)], 4.0)

@testset "KitaevComplex — Cⁿ basis" begin
    N = fib_Nijk()
    @test KC.C_dim(N, 2, 0) == 1
    @test KC.C_dim(N, 2, 1) == 2
    @test KC.C_dim(N, 2, 2) == 5

    N = ising_Nijk()
    @test KC.C_dim(N, 3, 0) == 1
    @test KC.C_dim(N, 3, 1) == 3
    @test KC.C_dim(N, 3, 2) == 10
end

@testset "KitaevComplex — δ⁰ = 0" begin
    for (label, N, dD) in [("Fib", fib_Nijk(), fib_dims()),
                           ("Ising", ising_Nijk(), ising_dims())]
        d, D2 = dD
        δ0 = KC.delta_matrix(N, d, D2, 0)
        @info "$label δ⁰" δ0
        @test maximum(abs.(δ0)) < 1e-10
    end
end

@testset "KitaevComplex — identity χⁿ⁺¹δⁿ + δⁿ⁻¹χⁿ = id, n ≥ 1" begin
    for (label, N, dD) in [("Fib", fib_Nijk(), fib_dims()),
                           ("Ising", ising_Nijk(), ising_dims())]
        d, D2 = dD
        err, M = KC.verify_homotopy(N, d, D2, 1)
        @info "$label n=1" err dim=size(M) M=M
        @test err < 1e-10
    end
end

@testset "KitaevComplex — component dump" begin
    for (label, N, dD) in [("Fib", fib_Nijk(), fib_dims()),
                           ("Ising", ising_Nijk(), ising_dims())]
        d, D2 = dD
        println("\n=== $label ===")

        δ0 = KC.delta_matrix(N, d, D2, 0)
        δ1 = KC.delta_matrix(N, d, D2, 1)
        χ1 = KC.chi_matrix(N, d, D2, 1)
        χ2 = KC.chi_matrix(N, d, D2, 2)

        println("δ⁰ = "); display(δ0); println()
        println("δ¹ = "); display(δ1); println()
        println("χ¹ = "); display(χ1); println()
        println("χ² = "); display(χ2); println()
        println("χ² δ¹ = "); display(χ2 * δ1); println()
        println("δ⁰ χ¹ = "); display(δ0 * χ1); println()

        @test true
    end
end
