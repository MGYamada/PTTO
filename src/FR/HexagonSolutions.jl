"""
Known exact hexagon solutions for small built-in fusion rules.

This file keeps closed-form R-symbol data separate from the hexagon solver's
finite-field preprocessing and public entry point.
"""

function _r_vector_from_channel_values(Nijk::Array{Int,3}, values::AbstractDict)
    r_pos, r_count = _braiding_block_positions(Nijk)
    K = parent(first(values).second)
    forward = [zero(K) for _ in 1:r_count]
    for ((i, j, k), positions) in r_pos
        val = get(values, (i, j, k), one(K))
        for pos in positions
            forward[pos] = val
        end
    end
    reverse = [zero(K) for _ in 1:r_count]
    for ((i, j, k), positions) in r_pos
        val = get(values, (j, i, k), get(values, (i, j, k), one(K)))
        invval = inv(val)
        for pos in positions
            reverse[pos] = invval
        end
    end
    return vcat(forward, reverse)
end

function _hexagon_solution_semion(ctx::CyclotomicContext)
    K, z = field(ctx), zeta(ctx)
    ctx.N % 4 == 0 || error("semion braiding requires conductor divisible by 4")
    r_ss = z^(ctx.N ÷ 4)
    return _r_vector_from_channel_values(
        begin
            N = zeros(Int, 2, 2, 2)
            N[1, 1, 1] = N[1, 2, 2] = N[2, 1, 2] = N[2, 2, 1] = 1
            N
        end,
        Dict((1, 1, 1) => one(K),
             (1, 2, 2) => one(K),
             (2, 1, 2) => one(K),
             (2, 2, 1) => r_ss))
end

function _hexagon_solution_fibonacci(ctx::CyclotomicContext)
    K, z = field(ctx), zeta(ctx)
    N = zeros(Int, 2, 2, 2)
    N[1, 1, 1] = N[1, 2, 2] = N[2, 1, 2] = N[2, 2, 1] = N[2, 2, 2] = 1
    return _r_vector_from_channel_values(
        N,
        Dict((1, 1, 1) => one(K),
             (1, 2, 2) => one(K),
             (2, 1, 2) => one(K),
             (2, 2, 1) => z^12,
             (2, 2, 2) => z^6))
end

function _hexagon_solution_ising(ctx::CyclotomicContext)
    K, z = field(ctx), zeta(ctx)
    return [
        -one(K),
        -K(2) * z^4,
        one(K),
        -z^4 // K(2),
        z^3,
        -z^7,
        one(K),
        one(K),
        one(K),
        one(K),
        -one(K),
        K(2) * z^4,
        one(K),
        z^4 // K(2),
        -z^5,
        z,
        one(K),
        one(K),
        one(K),
        one(K),
    ]
end

function _hexagon_solution_trivial_rank1(ctx::CyclotomicContext)
    K = field(ctx)
    return [one(K), one(K)]
end

function _hexagon_solution_ising(ctx::CyclotomicContext, Nijk::Array{Int,3})
    perm = _ising_label_perm_to_canonical(Nijk)
    perm === nothing && error("exact hexagon reconstruction is not implemented for this fusion rule")
    perm == [1, 2, 3] && return _hexagon_solution_ising(ctx)

    K = field(ctx)
    canonical_Nijk = _canonical_ising_fusion_rule()
    canonical_R = _hexagon_solution_ising(ctx)
    canonical_positions, canonical_total = _braiding_block_positions(canonical_Nijk)
    actual_positions, actual_total = _braiding_block_positions(Nijk)
    actual_R = Vector{typeof(one(K))}(undef, 2 * actual_total)
    for ((i, j, k), positions) in actual_positions
        canonical_key = (perm[i], perm[j], perm[k])
        canonical_positions_for_key = canonical_positions[canonical_key]
        for idx in eachindex(positions)
            actual_R[positions[idx]] = canonical_R[canonical_positions_for_key[idx]]
            actual_R[actual_total + positions[idx]] =
                canonical_R[canonical_total + canonical_positions_for_key[idx]]
        end
    end
    return actual_R
end
