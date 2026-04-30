"""
    build_orthogonality_equations(vars, n::Int, p::Int)
        -> Vector

Build polynomial equations for `U^T U = I` over `F_p`, where `vars`
encodes an `n×n` matrix in row-major order.
"""
function build_orthogonality_equations(vars, n::Int, p::Int)
    eqs = Any[]
    idx(i, j) = (i - 1) * n + j
    for i in 1:n
        for j in i:n
            acc = zero(vars[1])
            for k in 1:n
                acc += vars[idx(k, i)] * vars[idx(k, j)]
            end
            rhs = (i == j) ? one(vars[1]) : zero(vars[1])
            push!(eqs, acc - rhs)
        end
    end
    return eqs
end

"""
    is_orthogonal_mod_p(U::Matrix{Int}, p::Int) -> Bool

Check `U^T U ≡ I (mod p)`.
"""
function is_orthogonal_mod_p(U::Matrix{Int}, p::Int)
    n, m = size(U)
    n == m || return false
    for i in 1:n
        for j in 1:n
            acc = 0
            for k in 1:n
                acc = mod(acc + U[k, i] * U[k, j], p)
            end
            if i == j
                acc == 1 || return false
            else
                acc == 0 || return false
            end
        end
    end
    return true
end

function _block_var_index(n_block::Int, i::Int, j::Int)
    return (i - 1) * n_block + j
end

function _cayley_param_index(n::Int, i::Int, j::Int)
    i < j || error("require i < j")
    k = 0
    for a in 1:(i-1)
        k += n - a
    end
    k += j - i
    return k
end

function _A_entry_from_params(varsA, n::Int, i::Int, j::Int)
    if i == j
        return zero(varsA[1])
    elseif i < j
        return varsA[_cayley_param_index(n, i, j)]
    else
        return -varsA[_cayley_param_index(n, j, i)]
    end
end

"""
    build_cayley_link_equations(varsA, varsU, n_block::Int)
        -> Vector

Build equations encoding `(I + A)U = I - A`, where `A` is the antisymmetric
matrix parameterized by `varsA` (Cayley parameters) and `U` is an `n_block×n_block`
matrix represented by `varsU` in row-major order.
"""
function build_cayley_link_equations(varsA, varsU, n_block::Int)
    length(varsA) == div(n_block * (n_block - 1), 2) || error("varsA length mismatch")
    length(varsU) == n_block * n_block || error("varsU length mismatch")

    eqs = Any[]
    U(i, j) = varsU[_block_var_index(n_block, i, j)]
    for i in 1:n_block
        for j in 1:n_block
            lhs = zero(varsU[1])
            for k in 1:n_block
                Ipk = (i == k) ? one(varsU[1]) : zero(varsU[1])
                lhs += (Ipk + _A_entry_from_params(varsA, n_block, i, k)) * U(k, j)
            end
            rhs = (i == j ? one(varsU[1]) : zero(varsU[1])) - _A_entry_from_params(varsA, n_block, i, j)
            push!(eqs, lhs - rhs)
        end
    end
    return eqs
end

function _u_full_entry(varsU, indices_deg::Vector{Int}, i::Int, j::Int, p::Int)
    bi = findfirst(==(i), indices_deg)
    bj = findfirst(==(j), indices_deg)
    if bi !== nothing && bj !== nothing
        n_block = length(indices_deg)
        return varsU[_block_var_index(n_block, bi, bj)]
    end
    return i == j ? one(varsU[1]) : zero(varsU[1])
end

"""
    build_verlinde_unit_equations(varsU, varsW, S_atomic, indices_deg, p, unit_idx)
        -> Vector

Build polynomial equations for a fixed candidate unit `unit_idx`:
- orthogonality on the active block U
- invertibility witnesses `w_m * S′[unit_idx,m] - 1 = 0`
- unit axiom `N_{unit_idx,i}^k = δ_{i,k}`

where `S′ = U^T S U` and `U` acts only on `indices_deg`.
"""
function build_verlinde_unit_equations(varsU, varsW,
                                       S_atomic::Matrix{Int},
                                       indices_deg::Vector{Int},
                                       p::Int, unit_idx::Int)
    r = size(S_atomic, 1)
    n_block = length(indices_deg)
    size(S_atomic, 2) == r || error("S_atomic must be square")
    length(varsU) == n_block * n_block || error("varsU length mismatch")
    length(varsW) == r || error("varsW length mismatch")

    eqs = Any[]

    # O(n_block): U^T U = I
    append!(eqs, build_orthogonality_equations(varsU, n_block, p))

    # Build symbolic S' = U^T S U
    Sprime = Matrix{Any}(undef, r, r)
    for i in 1:r
        for j in 1:r
            acc = zero(varsU[1])
            for a in 1:r
                Uai = _u_full_entry(varsU, indices_deg, a, i, p)
                # Split second contraction for better readability:
                for b in 1:r
                    Ubj = _u_full_entry(varsU, indices_deg, b, j, p)
                    acc += Uai * S_atomic[a, b] * Ubj
                end
            end
            Sprime[i, j] = acc
        end
    end

    # Witness inverse entries for S'[u,m]
    for m in 1:r
        push!(eqs, varsW[m] * Sprime[unit_idx, m] - one(varsU[1]))
    end

    # Unit axiom only: N_{u,i}^k = δ_{i,k}
    for i in 1:r
        for k in 1:r
            acc = zero(varsU[1])
            for m in 1:r
                acc += Sprime[unit_idx, m] * Sprime[i, m] * Sprime[k, m] * varsW[m]
            end
            push!(eqs, acc - ((i == k) ? one(varsU[1]) : zero(varsU[1])))
        end
    end

    return eqs
end

"""
    verlinde_find_unit(S_Fp::Matrix{Int}, p::Int; threshold::Int = 5)
        -> Union{Nothing, Tuple{Int, Array{Int,3}}}

Given an F_p modular datum S, try each basis index u as the candidate
"unit object". For unit u to be valid:
1. S[u, m] ≠ 0 for all m (so 1/S[u,m] is defined).
2. N_{u, i}^{k} = δ_{i,k} (unit fusion is identity).
3. All N_{i,j}^k small non-negative integers (|n| ≤ threshold).

If found, returns (u, N) where N is the r×r×r tensor. Otherwise nothing.

The Verlinde formula (assuming S^2 = I, self-dual MTC):
    N_{ij}^k = Σ_m (S[i,m] · S[j,m] · S[k,m]) / S[u,m]
"""
function verlinde_find_unit(S_Fp::Matrix{Int}, p::Int; threshold::Int = 5)
    r = size(S_Fp, 1)

    for u in 1:r
        S_u_row = [S_Fp[u, m] for m in 1:r]
        any(x -> x == 0, S_u_row) && continue

        try
            S_u_inv = [invmod(x, p) for x in S_u_row]

            # First check unit axiom: N[u, i, k] = δ_{i,k}
            is_unit = true
            for i in 1:r
                for k in 1:r
                    val = 0
                    for m in 1:r
                        term = mod(S_Fp[u, m] * S_Fp[i, m], p)
                        term = mod(term * S_Fp[k, m], p)
                        term = mod(term * S_u_inv[m], p)
                        val = mod(val + term, p)
                    end
                    expected = (i == k) ? 1 : 0
                    if val != expected
                        is_unit = false
                        break
                    end
                end
                is_unit || break
            end
            is_unit || continue

            # Compute full tensor, check non-negative integer
            N = zeros(Int, r, r, r)
            ok = true
            for i in 1:r
                for j in 1:r
                    for k in 1:r
                        val = 0
                        for m in 1:r
                            term = mod(S_Fp[i, m] * S_Fp[j, m], p)
                            term = mod(term * S_Fp[k, m], p)
                            term = mod(term * S_u_inv[m], p)
                            val = mod(val + term, p)
                        end
                        s = signed_Fp(val, p)
                        if s < 0 || abs(s) > threshold
                            ok = false
                            break
                        end
                        N[i, j, k] = s
                    end
                    ok || break
                end
                ok || break
            end
            ok && return (u, N)
        catch
            continue
        end
    end
    return nothing
end

"""
    passes_unit_axiom(S_Fp::Matrix{Int}, p::Int, u::Int) -> Bool

Fast check of unit axiom only:
`N_{u,i}^k = δ_{i,k}` using the Verlinde formula over `F_p`.
"""
function passes_unit_axiom(S_Fp::Matrix{Int}, p::Int, u::Int)
    r = size(S_Fp, 1)
    S_u_row = [S_Fp[u, m] for m in 1:r]
    any(x -> x == 0, S_u_row) && return false

    local S_u_inv
    try
        S_u_inv = [invmod(x, p) for x in S_u_row]
    catch
        return false
    end

    for i in 1:r
        for k in 1:r
            val = 0
            for m in 1:r
                term = mod(S_Fp[u, m] * S_Fp[i, m], p)
                term = mod(term * S_Fp[k, m], p)
                term = mod(term * S_u_inv[m], p)
                val = mod(val + term, p)
            end
            expected = (i == k) ? 1 : 0
            val == expected || return false
        end
    end
    return true
end
