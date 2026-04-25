using Test
using LinearAlgebra
using ACMG

@testset "Pipeline modular-data roundtrip logic" begin
    @testset "_modular_data_roundtrip_up_to_galois aligns S/T under same Galois action" begin
        # Fibonacci fusion ring (rank 2, τ self-dual).
        Nijk = zeros(Int, 2, 2, 2)
        Nijk[1, 1, 1] = 1
        Nijk[1, 2, 2] = 1
        Nijk[2, 1, 2] = 1
        Nijk[2, 2, 1] = 1
        Nijk[2, 2, 2] = 1

        # Same S-from-T reconstruction formula as the pipeline helper.
        function reconstruct_S_from_T_fib(Tvals::Vector{ComplexF64})
            r = size(Nijk, 1)
            A = zeros(Float64, r, r)
            for i in 1:r
                A .+= Float64.(Nijk[i, :, :])
            end
            eig = eigen(A)
            idx = argmax(real(eig.values))
            d = abs.(real(eig.vectors[:, idx]))
            d ./= d[1]
            D = sqrt(sum(d .^ 2))

            S = Matrix{ComplexF64}(undef, r, r)
            for i in 1:r, j in 1:r
                acc = 0.0 + 0.0im
                for k in 1:r
                    Nijk[i, j, k] == 0 && continue
                    acc += Nijk[i, j, k] * (Tvals[k] / (Tvals[i] * Tvals[j])) * d[k]
                end
                S[i, j] = acc / D
            end
            return S
        end

        N = 5
        T_target = ComplexF64[1.0, exp(2π * im / N)]
        S_target = reconstruct_S_from_T_fib(T_target)

        # Build an R vector whose inferred twist is T_target^2.
        T_from = ComplexF64[1.0, exp(2π * im * 2 / N)]
        t2 = 2
        m_by_triple = Dict(
            (1, 1, 1) => 0,
            (1, 2, 2) => 0,
            (2, 1, 2) => 0,
            (2, 2, 1) => mod(2 * t2, N),
            (2, 2, 2) => t2,
        )
        positions = block_positions_R(Nijk)
        r_var_count = sum(Nijk[i, j, k]^2 for i in 1:2, j in 1:2, k in 1:2 if Nijk[i, j, k] > 0)
        R_fwd = Vector{ComplexF64}(undef, r_var_count)
        for (ijk, m) in m_by_triple
            idx = positions[ijk][1]
            # infer_T_candidates_from_R uses z = R^2
            R_fwd[idx] = exp((π * im * m) / N)
        end

        out = ACMG._modular_data_roundtrip_up_to_galois(ComplexF64[1.0], R_fwd,
                                                        Nijk, S_target, T_target, N)
        @test out.ok
        @test out.best_a == 2
        @test out.T_max < 1e-8
        @test out.S_max < 1e-8
        @test isapprox(out.T_best[2], T_from[2]; atol = 1e-8)
    end
end

@testset "Pipeline fusion-rule aggregation key" begin
    @testset "fusion_rule_key is invariant under non-unit relabeling and tensor type" begin
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
        @test ACMG.fusion_rule_key(Int32.(Nijk)) == ACMG.fusion_rule_key(Nijk)
    end
end
