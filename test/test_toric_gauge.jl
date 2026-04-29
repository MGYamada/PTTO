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
    @test r_symbol_weight(2, 3, 4) == Dict((3, 2, 4) => 1, (2, 3, 4) => -1)
    @test residual_gauge_orders([2 0; 0 0]) == [2]
end
