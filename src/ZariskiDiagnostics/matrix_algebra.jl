struct MatrixAlgebraDiagnostics
    p::Union{Int, Nothing}
    dimension::Int
    ambient_dimension::Int
    is_full_matrix_algebra::Bool
    basis::Vector{Matrix{Int}}
    truncated::Bool
end

function _rref_basis_add!(basis_vecs::Vector{Vector{Int}}, basis_mats::Vector{Matrix{Int}},
                          M::Matrix{Int}, p::Int)
    v = mod.(vec(M), p)
    rows = copy(basis_vecs)
    push!(rows, v)
    old = _rank_mod_p(basis_vecs, length(v), p)
    new = _rank_mod_p(rows, length(v), p)
    if new > old
        push!(basis_vecs, v)
        push!(basis_mats, mod.(M, p))
        return true
    end
    return false
end

function _rank_mod_p(rows::Vector{Vector{Int}}, ncols::Int, p::Int)
    isempty(rows) && return 0
    M = reduce(vcat, [reshape(mod.(r, p), 1, ncols) for r in rows])
    rank = 0
    r = 1
    for c in 1:ncols
        piv = findfirst(i -> M[i, c] != 0, r:size(M, 1))
        piv === nothing && continue
        piv = r + piv - 1
        M[r, :], M[piv, :] = M[piv, :], M[r, :]
        invp = invmod(M[r, c], p)
        M[r, :] = mod.(M[r, :] .* invp, p)
        for i in 1:size(M, 1)
            i == r && continue
            f = M[i, c]
            f != 0 && (M[i, :] = mod.(M[i, :] .- f .* M[r, :], p))
        end
        rank += 1
        r += 1
        r > size(M, 1) && break
    end
    return rank
end

function generated_matrix_algebra(sigmas::AbstractVector{<:AbstractMatrix{<:Integer}}, p::Int; max_dimension = nothing)
    isempty(sigmas) && error("generated_matrix_algebra requires at least one generator")
    n = size(sigmas[1], 1)
    limit = max_dimension === nothing ? n^2 : Int(max_dimension)
    mats = Matrix{Int}[]
    vecs = Vector{Int}[]
    _rref_basis_add!(vecs, mats, _eye(n, p), p)
    for g in sigmas
        _rref_basis_add!(vecs, mats, mod.(Int.(g), p), p)
    end
    head = 1
    truncated = false
    while head <= length(mats)
        A = mats[head]; head += 1
        current = copy(mats)
        for B in current
            for C in (matmul_mod(A, B, p), matmul_mod(B, A, p))
                _rref_basis_add!(vecs, mats, C, p)
                if length(mats) >= limit
                    truncated = length(mats) < n^2
                    return MatrixAlgebraDiagnostics(p, length(mats), n^2, length(mats) == n^2, mats, truncated)
                end
            end
        end
    end
    return MatrixAlgebraDiagnostics(p, length(mats), n^2, length(mats) == n^2, mats, truncated)
end

generated_matrix_algebra(br::FiniteFieldBraidRepresentation; kwargs...) =
    generated_matrix_algebra(br.generators, br.p; kwargs...)
