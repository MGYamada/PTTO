function invariant_subspace_search(sigmas::AbstractVector{<:AbstractMatrix{<:Integer}}, p::Int; max_vectors = 20000)
    isempty(sigmas) && return (found = false, dimension = 0, vector = nothing, searched = 0)
    n = size(sigmas[1], 1)
    searched = 0
    for x in Iterators.product(ntuple(_ -> 0:(p - 1), n)...)
        v = collect(Int, x)
        all(iszero, v) && continue
        searched += 1
        searched > max_vectors && return (found = false, dimension = nothing, vector = nothing, searched = searched)
        ok = true
        for A in sigmas
            Av = mod.(Int.(A) * v, p)
            λ = nothing
            for i in 1:n
                if v[i] != 0
                    λ = mod(Av[i] * invmod(v[i], p), p)
                    break
                end
            end
            if λ === nothing || any(mod(Av[i] - λ * v[i], p) != 0 for i in 1:n)
                ok = false
                break
            end
        end
        ok && return (found = true, dimension = 1, vector = v, searched = searched)
    end
    return (found = false, dimension = 0, vector = nothing, searched = searched)
end
