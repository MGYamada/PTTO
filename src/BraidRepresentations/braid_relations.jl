struct BraidRelationCheck
    ok::Bool
    failures::Vector{NamedTuple}
    tolerance::Any
end

_mat_close(A, B; atol = 0) = atol == 0 ? A == B :
    size(A) == size(B) && all(abs(A[i] - B[i]) <= atol for i in eachindex(A))

function check_braid_relations(sigmas::AbstractVector; atol = 0)
    failures = NamedTuple[]
    m = length(sigmas)
    for i in 1:m, j in (i+2):m
        _mat_close(sigmas[i] * sigmas[j], sigmas[j] * sigmas[i]; atol = atol) ||
            push!(failures, (kind = :far_commutation, i = i, j = j))
    end
    for i in 1:(m - 1)
        lhs = sigmas[i] * sigmas[i + 1] * sigmas[i]
        rhs = sigmas[i + 1] * sigmas[i] * sigmas[i + 1]
        _mat_close(lhs, rhs; atol = atol) ||
            push!(failures, (kind = :yang_baxter, i = i, j = i + 1))
    end
    return BraidRelationCheck(isempty(failures), failures, atol)
end

check_braid_relations(br::BraidRepresentation; kwargs...) =
    check_braid_relations(br.generators; kwargs...)

function check_braid_relations(sigmas::AbstractVector{<:AbstractMatrix{<:Integer}}, p::Int)
    failures = NamedTuple[]
    m = length(sigmas)
    for i in 1:m, j in (i+2):m
        matmul_mod(sigmas[i], sigmas[j], p) == matmul_mod(sigmas[j], sigmas[i], p) ||
            push!(failures, (kind = :far_commutation, i = i, j = j))
    end
    for i in 1:(m - 1)
        lhs = matmul_mod(matmul_mod(sigmas[i], sigmas[i + 1], p), sigmas[i], p)
        rhs = matmul_mod(matmul_mod(sigmas[i + 1], sigmas[i], p), sigmas[i + 1], p)
        lhs == rhs || push!(failures, (kind = :yang_baxter, i = i, j = i + 1))
    end
    return BraidRelationCheck(isempty(failures), failures, 0)
end
