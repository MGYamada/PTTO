"""
Known exact pentagon solutions for small built-in fusion rules.

The public solver still decides when these tables are applicable; this file
only keeps the closed-form F-symbol data out of the solver orchestration.
"""

function _associator_coordinate_slots(Nijk::Array{Int,3}, K)
    r = size(Nijk, 1)
    one_vec = zeros(Int, r)
    one_vec[1] = 1
    C = TensorCategories.six_j_category(K, Nijk)
    C.one = one_vec
    slots = Tuple{Int,Int,Int,Int,Int,Int}[]
    for i in 1:r, j in 1:r, k in 1:r, o in 1:r
        sum(one_vec[[i, j, k]]) > 0 && continue
        rows, cols = size(C.ass[i, j, k, o])
        for a in 1:rows, b in 1:cols
            push!(slots, (i, j, k, o, a, b))
        end
    end
    # TensorCategories' assigner consumes coordinates with pop!, so the
    # external F-vector order is the reverse of the traversal order.
    return reverse(slots)
end

function _default_context_from_kwargs(; context = nothing, conductor = nothing, N = nothing)
    context !== nothing && return context
    n = conductor === nothing ? N : conductor
    n === nothing && error("a CyclotomicContext or conductor N is required for exact Phase 4")
    return CyclotomicContext(n)
end

function _pentagon_solution_semion(ctx::CyclotomicContext)
    K = field(ctx)
    return [-one(K)]
end

function _pentagon_solution_fibonacci(ctx::CyclotomicContext)
    K = field(ctx)
    sqrt5 = _sqrt5(ctx)
    a = (sqrt5 - one(K)) // K(2)
    return [-a, one(K), a, a, one(K)]
end

function _pentagon_solution_ising(ctx::CyclotomicContext)
    K = field(ctx)
    sqrt2 = _sqrt2(ctx)
    h = inv(sqrt2)
    return [
        one(K),
        sqrt2 // K(4),
        -one(K),
        sqrt2,
        K(2),
        h,
        -K(1)//K(2),
        K(1)//K(2),
        sqrt2,
        one(K),
        -sqrt2,
        one(K),
        one(K),
        h,
    ]
end

function _known_pentagon_solution(Nijk::Array{Int,3}, ctx::CyclotomicContext)
    _is_trivial_rank1_fusion(Nijk) && return elem_type(field(ctx))[]
    _is_semion_fusion(Nijk) && return _pentagon_solution_semion(ctx)
    _is_fibonacci_fusion(Nijk) && return _pentagon_solution_fibonacci(ctx)

    perm = _ising_label_perm_to_canonical(Nijk)
    perm === nothing && return nothing
    perm == [1, 2, 3] && return _pentagon_solution_ising(ctx)

    K = field(ctx)
    canonical_Nijk = _canonical_ising_fusion_rule()
    canonical_slots = _associator_coordinate_slots(canonical_Nijk, K)
    canonical_F = _pentagon_solution_ising(ctx)
    slot_values = Dict(canonical_slots[i] => canonical_F[i]
                       for i in eachindex(canonical_slots))
    actual_slots = _associator_coordinate_slots(Nijk, K)
    return [slot_values[(perm[i], perm[j], perm[k], perm[o], a, b)]
            for (i, j, k, o, a, b) in actual_slots]
end
