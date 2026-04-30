# Local helpers (small utilities, will move to FpArith if reused broadly)
function I_n(n::Int)
    M = zeros(Int, n, n)
    for i in 1:n
        M[i, i] = 1
    end
    return M
end

function transpose_mod(A::Matrix{Int}, p::Int)
    n, m = size(A)
    B = zeros(Int, m, n)
    for i in 1:m
        for j in 1:n
            B[i, j] = A[j, i]
        end
    end
    return B
end

"""
    inverse_mod_p(M::Matrix{Int}, p::Int) -> Union{Nothing, Matrix{Int}}

Compute the inverse of an n×n matrix mod p via Gauss-Jordan elimination.
Returns `nothing` if M is singular mod p.
"""
function inverse_mod_p(M::Matrix{Int}, p::Int)
    n = size(M, 1)
    n == size(M, 2) || error("M must be square")

    # Augment [M | I]
    aug = zeros(Int, n, 2n)
    for i in 1:n
        for j in 1:n
            aug[i, j] = mod(M[i, j], p)
        end
        aug[i, n + i] = 1
    end

    # Gauss-Jordan
    for col in 1:n
        # Find pivot
        pivot_row = 0
        for r in col:n
            if aug[r, col] != 0
                pivot_row = r
                break
            end
        end
        pivot_row == 0 && return nothing  # singular

        if pivot_row != col
            # Swap rows
            for j in 1:(2n)
                aug[col, j], aug[pivot_row, j] = aug[pivot_row, j], aug[col, j]
            end
        end

        # Scale pivot row
        pivot = aug[col, col]
        pivot_inv = invmod(pivot, p)
        for j in 1:(2n)
            aug[col, j] = mod(aug[col, j] * pivot_inv, p)
        end

        # Eliminate other rows
        for r in 1:n
            r == col && continue
            factor = aug[r, col]
            factor == 0 && continue
            for j in 1:(2n)
                aug[r, j] = mod(aug[r, j] - factor * aug[col, j], p)
            end
        end
    end

    return aug[:, (n+1):(2n)]
end

"""
    signed_Fp(x::Int, p::Int) -> Int

Lift x ∈ [0, p) to the symmetric interval [-⌊p/2⌋, ⌊p/2⌋].
"""
function signed_Fp(x::Int, p::Int)
    return x <= p ÷ 2 ? x : x - p
end
