using Test
using LinearAlgebra
using ACMG

@testset "Pipeline modular-data roundtrip logic" begin
    @testset "_modular_data_roundtrip reconstructs S/T from Kitaev (211)(223)" begin
        # Fibonacci fusion ring (rank 2, τ self-dual).
        Nijk = zeros(Int, 2, 2, 2)
        Nijk[1, 1, 1] = 1
        Nijk[1, 2, 2] = 1
        Nijk[2, 1, 2] = 1
        Nijk[2, 2, 1] = 1
        Nijk[2, 2, 2] = 1

        φ = (1 + sqrt(5.0)) / 2
        D = sqrt(1 + φ^2)
        S_target = ComplexF64[
            1 / D   φ / D;
            φ / D  -1 / D
        ]
        T_target = ComplexF64[1.0, exp(4π * im / 5)]

        # One valid Fibonacci braiding branch.
        positions = block_positions_R(Nijk)
        r_var_count = sum(Nijk[i, j, k]^2 for i in 1:2, j in 1:2, k in 1:2 if Nijk[i, j, k] > 0)
        R_fwd = fill(0.0 + 0.0im, r_var_count)
        R_by_triple = Dict(
            (1, 1, 1) => 1.0 + 0.0im,
            (1, 2, 2) => 1.0 + 0.0im,
            (2, 1, 2) => 1.0 + 0.0im,
            (2, 2, 1) => exp(-4π * im / 5),
            (2, 2, 2) => exp(3π * im / 5),
        )
        for (ijk, z) in R_by_triple
            idx = positions[ijk][1]
            R_fwd[idx] = z
        end

        out = ACMG._modular_data_roundtrip(ComplexF64[1.0], R_fwd,
                                           Nijk, S_target, T_target, 20)
        @test out.ok
        @test out.T_max < 1e-8
        @test out.S_max < 1e-8
        @test isapprox(out.T_from[2], T_target[2]; atol = 1e-8)

        R_full = vcat(R_fwd, fill(1.0 + 0.0im, r_var_count))
        out_full = ACMG._modular_data_roundtrip(ComplexF64[1.0], R_full,
                                                Nijk, S_target, T_target, 20)
        @test out_full.ok
        @test out_full.T_max < 1e-8
        @test out_full.S_max < 1e-8
    end

    @testset "multiplicity>1 twist uses trace over the fusion channel" begin
        Nijk = zeros(Int, 2, 2, 2)
        Nijk[1, 1, 1] = 1
        Nijk[1, 2, 2] = 1
        Nijk[2, 1, 2] = 1
        Nijk[2, 2, 1] = 2

        positions = block_positions_R(Nijk)
        r_var_count = sum(Nijk[i, j, k]^2 for i in 1:2, j in 1:2, k in 1:2 if Nijk[i, j, k] > 0)
        R_fwd = fill(0.0 + 0.0im, r_var_count)
        for ijk in ((1, 1, 1), (1, 2, 2), (2, 1, 2))
            R_fwd[positions[ijk][1]] = 1.0 + 0.0im
        end

        θ = exp(2π * im / 20)
        pos = positions[(2, 2, 1)]
        R_fwd[pos[1]] = θ
        R_fwd[pos[2]] = 0.0 + 0.0im
        R_fwd[pos[3]] = 0.0 + 0.0im
        R_fwd[pos[4]] = θ

        @test isapprox(ACMG._monodromy_trace(R_fwd, Nijk, 2, 2, 1),
                       2 * θ^2; atol = 1e-12)

        T = ACMG._twists_from_braiding_trace(R_fwd, Nijk, ComplexF64[1.0, 1.0])
        @test isapprox(T[2], θ; atol = 1e-12)
    end
end

@testset "Pipeline fusion-rule aggregation key" begin
    @testset "fusion_rule_key/canonical_rule are invariant under relabeling and tensor type" begin
        Nijk = zeros(Int, 3, 3, 3)
        for i in 1:3
            Nijk[1, i, i] = 1
            Nijk[i, 1, i] = 1
        end
        # 2 ⊗ 2 = 1 ⊕ 2, 3 ⊗ 3 = 1 ⊕ 3, 2 ⊗ 3 = 3
        Nijk[2, 2, 1] = 1
        Nijk[2, 2, 2] = 1
        Nijk[3, 3, 1] = 1
        Nijk[3, 3, 3] = 1
        Nijk[2, 3, 3] = 1
        Nijk[3, 2, 3] = 1

        perm = [1, 3, 2]
        N_perm = Nijk[perm, perm, perm]
        @test ACMG.fusion_rule_key(Nijk) == ACMG.fusion_rule_key(N_perm)
        @test ACMG.canonical_rule(Nijk) == ACMG.canonical_rule(N_perm)
        @test ACMG.fusion_rule_key(Int32.(Nijk)) == ACMG.fusion_rule_key(Nijk)

        # Full relabeling including unit position should still match.
        full_perm = [2, 1, 3]
        N_full_perm = Nijk[full_perm, full_perm, full_perm]
        @test ACMG.canonical_rule(Nijk) == ACMG.canonical_rule(N_full_perm)
    end
end
