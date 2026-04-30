function _dedupe_matrices(mats::Vector{Matrix{Int}})
    seen = Set{String}()
    out = Matrix{Int}[]
    for M in mats
        key = join(vec(M), ",")
        if !(key in seen)
            push!(seen, key)
            push!(out, M)
        end
    end
    return out
end

function _sort_matrices_lex(mats::Vector{Matrix{Int}})
    return sort(mats; by = M -> join(vec(M), ","))
end
