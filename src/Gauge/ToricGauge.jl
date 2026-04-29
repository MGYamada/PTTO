"""
Smith normal form summaries for toric gauge actions.
"""

function _smith_invariant_factors(W::AbstractMatrix{<:Integer})
    m, n = size(W)
    (m == 0 || n == 0) && return BigInt[]
    entries = [BigInt(W[i, j]) for i in 1:m for j in 1:n]
    S = AbstractAlgebra.snf(matrix(ZZ, m, n, entries))
    invariants = BigInt[]
    for i in 1:min(m, n)
        d = abs(BigInt(S[i, i]))
        d == 0 && continue
        push!(invariants, d)
    end
    return invariants
end

"""
    smith_gauge_split(W)

Summarize the Smith normal form of the integer gauge weight matrix.  The
returned named tuple records the invariant factors, rank of the effective
torus action, free rank of the ineffective kernel, and finite residual
orders from non-unit invariant factors.
"""
function smith_gauge_split(W::AbstractMatrix{<:Integer})
    invariants = _smith_invariant_factors(W)
    rankW = length(invariants)
    nparams = size(W, 2)
    orders = Int[d for d in invariants if d > 1]
    return (invariant_factors = Int.(invariants),
            effective_rank = rankW,
            ineffective_rank = nparams - rankW,
            residual_orders = orders)
end

"""
    ineffective_kernel_rank(W)

Return the dimension of the connected kernel of the toric gauge action.
"""
ineffective_kernel_rank(W::AbstractMatrix{<:Integer}) = smith_gauge_split(W).ineffective_rank

"""
    residual_gauge_orders(W)

Return the finite diagonal kernel orders visible in the Smith normal form.
"""
residual_gauge_orders(W::AbstractMatrix{<:Integer}) = smith_gauge_split(W).residual_orders
