using Test
using ACMG

function _tg_semion_fusion()
    N = zeros(Int, 2, 2, 2)
    N[1, 1, 1] = 1
    N[1, 2, 2] = 1
    N[2, 1, 2] = 1
    N[2, 2, 1] = 1
    return N
end

function _tg_fibonacci_fusion()
    N = _tg_semion_fusion()
    N[2, 2, 2] = 1
    return N
end

function _tg_ising_fusion()
    N = zeros(Int, 3, 3, 3)
    for a in 1:3
        N[1, a, a] = 1
        N[a, 1, a] = 1
    end
    N[2, 2, 1] = 1
    N[2, 2, 3] = 1
    N[2, 3, 2] = 1
    N[3, 2, 2] = 1
    N[3, 3, 1] = 1
    return N
end

function _object_phase_vector(params, obj::Int)
    return [((a == obj) ? 1 : 0) + ((b == obj) ? 1 : 0) - ((c == obj) ? 1 : 0)
            for (a, b, c) in params]
end

@testset "Toric multiplicity-free gauge weights" begin
    cases = [
        (name = :semion, fusion = _tg_semion_fusion(), rank = 2),
        (name = :fibonacci, fusion = _tg_fibonacci_fusion(), rank = 2),
        (name = :ising, fusion = _tg_ising_fusion(), rank = 3),
    ]

    for case in cases
        @testset "$(case.name)" begin
            params = gauge_parameters(case.fusion)
            coords = symbol_coordinates(case.fusion)
            W = gauge_weight_matrix(case.fusion)
            split = smith_gauge_split(W)

            @test size(W) == (length(coords), length(params))
            @test split.effective_rank + split.ineffective_rank == length(params)
            @test ineffective_kernel_rank(W) == split.ineffective_rank
            @test residual_gauge_orders(W) == split.residual_orders
            @test split.ineffective_rank >= case.rank

            for obj in 1:case.rank
                @test W * _object_phase_vector(params, obj) == zeros(Int, size(W, 1))
            end

            p = 17
            symbol_data = (coordinates = coords,
                           values = collect(1:length(coords)),
                           parameters = params)
            gauge = Dict(ch => i + 2 for (i, ch) in enumerate(params))
            moved = apply_gauge_mod_p(symbol_data, gauge, p)
            moved_typed = apply_gauge_mod_p(symbol_data,
                                            GaugeParameters(Dict((a, b, c, 1) => v
                                                                 for ((a, b, c), v) in gauge)),
                                            p)
            moved_action = apply_gauge_mod_p(symbol_data, GaugeAction(gauge; field = :F_17), p)
            @test moved_typed.values == moved.values
            @test moved_action.values == moved.values
            returned = apply_gauge_mod_p(moved,
                                         Dict(ch => invmod(v, p) for (ch, v) in gauge),
                                         p)
            @test returned.values == mod.(symbol_data.values, p)
            @test stabilizer_size_mod_p(symbol_data, case.fusion, p) >= BigInt(p - 1)^case.rank
            @test stacky_weight_mod_p(symbol_data, case.fusion, p) ==
                  1 // stabilizer_size_mod_p(symbol_data, case.fusion, p)
        end
    end

    @test f_symbol_weight(1, 2, 3, 4, 5, 6) ==
          Dict((1, 2, 5) => 1, (5, 3, 4) => 1,
               (2, 3, 6) => -1, (1, 6, 4) => -1)
    @test r_symbol_weight(2, 3, 4) == Dict((2, 3, 4) => 1, (3, 2, 4) => -1)
    @test residual_gauge_orders([2 0; 0 0]) == [2]
end

@testset "SNF toric gauge normal forms feed braid representations" begin
    fr = fibonacci_fr_data_mod_p(101)
    p = fr_metadata(fr)[:p]

    data = toric_gauge_data(fr)
    @test data.split == smith_gauge_split(data.weight_matrix)
    @test length(data.coordinates) == length(F_values(fr)) + length(R_values(fr))

    slice = ACMG.toric_gauge_slice(fusion_rule(fr))
    gf = fr_metadata(fr)[:gauge_fixing]
    @test gf[:gauge_fix_method] == :toric_snf_f_slice
    @test gf[:fixed_f_indices] == slice.fixed_indices
    @test F_values(fr)[gf[:fixed_f_indices]] == fill(FpElem(1, p), length(gf[:fixed_f_indices]))

    normal = toric_gauge_normal_form(fr)
    @test normal isa ToricGaugeNormalFormResult
    @test normal.split.effective_rank > 0
    @test normal.stabilizer_size == stabilizer_size_mod_p(
        (coordinates = fr_symbol_coordinates(normal.frdata),
         values = vcat(F_values(normal.frdata), R_values(normal.frdata)),
         parameters = gauge_parameters(normal.frdata)),
        fusion_rule(normal.frdata),
        p)
    @test normal.stacky_weight == 1 // normal.stabilizer_size
    @test verify_FRData(normal.frdata)
    @test is_gauge_fixed(normal.frdata)

    br = braid_representation(fr, [2, 2, 2], 2)
    @test br.fr_data == fr
    σ1, σ2 = braid_generators(br)
    @test σ1 * σ2 * σ1 == σ2 * σ1 * σ2
end
