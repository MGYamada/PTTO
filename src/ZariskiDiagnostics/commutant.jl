struct CommutantDiagnostics
    p::Union{Int, Nothing}
    dimension::Int
    basis::Vector{Matrix{Int}}
    appears_absolutely_irreducible::Bool
end

function _nullspace_mod_p(A::Matrix{Int}, p::Int)
    m, n = size(A)
    M = mod.(copy(A), p)
    pivots = Int[]
    r = 1
    for c in 1:n
        piv = findfirst(i -> M[i, c] != 0, r:m)
        piv === nothing && continue
        piv = r + piv - 1
        M[r, :], M[piv, :] = M[piv, :], M[r, :]
        invp = invmod(M[r, c], p)
        M[r, :] = mod.(M[r, :] .* invp, p)
        for i in 1:m
            i == r && continue
            f = M[i, c]
            f != 0 && (M[i, :] = mod.(M[i, :] .- f .* M[r, :], p))
        end
        push!(pivots, c)
        r += 1
        r > m && break
    end
    free = [c for c in 1:n if !(c in pivots)]
    basis = Vector{Int}[]
    for f in free
        x = zeros(Int, n)
        x[f] = 1
        for (row, pc) in enumerate(pivots)
            x[pc] = mod(-M[row, f], p)
        end
        push!(basis, x)
    end
    return basis
end

function commutant(sigmas::AbstractVector{<:AbstractMatrix{<:Integer}}, p::Int)
    isempty(sigmas) && error("commutant requires at least one generator")
    n = size(sigmas[1], 1)
    rows = Vector{Int}[]
    for A0 in sigmas
        A = mod.(Int.(A0), p)
        for i in 1:n, j in 1:n
            row = zeros(Int, n*n)
            for k in 1:n
                row[(i - 1) * n + k] = mod(row[(i - 1) * n + k] + A[k, j], p)
                row[(k - 1) * n + j] = mod(row[(k - 1) * n + j] - A[i, k], p)
            end
            push!(rows, row)
        end
    end
    Aeq = reduce(vcat, [reshape(r, 1, n*n) for r in rows])
    ns = _nullspace_mod_p(Aeq, p)
    mats = [reshape(v, n, n) for v in ns]
    return CommutantDiagnostics(p, length(mats), mats, length(mats) == 1)
end

commutant(br::FiniteFieldBraidRepresentation) = commutant(br.generators, br.p)
