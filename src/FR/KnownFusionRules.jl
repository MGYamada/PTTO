"""
Recognition helpers for small built-in fusion rules used by exact F/R solvers.

These predicates are internal glue for the current exact Phase-4 tables; they
do not add new fusion rules or public APIs.
"""

function _is_semion_fusion(Nijk::Array{Int,3})
    size(Nijk) == (2, 2, 2) || return false
    target = zeros(Int, 2, 2, 2)
    target[1, 1, 1] = 1
    target[1, 2, 2] = 1
    target[2, 1, 2] = 1
    target[2, 2, 1] = 1
    return Nijk == target
end

function _is_fibonacci_fusion(Nijk::Array{Int,3})
    size(Nijk) == (2, 2, 2) || return false
    target = zeros(Int, 2, 2, 2)
    target[1, 1, 1] = 1
    target[1, 2, 2] = 1
    target[2, 1, 2] = 1
    target[2, 2, 1] = 1
    target[2, 2, 2] = 1
    return Nijk == target
end

function _is_trivial_rank1_fusion(Nijk::Array{Int,3})
    return size(Nijk) == (1, 1, 1) && Nijk[1, 1, 1] == 1
end

function _is_ising_fusion(Nijk::Array{Int,3})
    size(Nijk) == (3, 3, 3) || return false
    target = zeros(Int, 3, 3, 3)
    for a in 1:3
        target[1, a, a] = 1
        target[a, 1, a] = 1
    end
    target[2, 2, 1] = 1
    target[2, 2, 3] = 1
    target[2, 3, 2] = 1
    target[3, 2, 2] = 1
    target[3, 3, 1] = 1
    return Nijk == target
end

function _ising_label_perm_to_canonical(Nijk::Array{Int,3})
    _is_ising_fusion(Nijk) && return [1, 2, 3]

    # The pipeline can reconstruct Ising with labels (1, ψ, σ), while the
    # closed-form tables above use (1, σ, ψ).
    size(Nijk) == (3, 3, 3) || return nothing
    target = zeros(Int, 3, 3, 3)
    for a in 1:3
        target[1, a, a] = 1
        target[a, 1, a] = 1
    end
    target[2, 2, 1] = 1
    target[2, 3, 3] = 1
    target[3, 2, 3] = 1
    target[3, 3, 1] = 1
    target[3, 3, 2] = 1
    Nijk == target && return [1, 3, 2]
    return nothing
end

function _canonical_ising_fusion_rule()
    N = zeros(Int, 3, 3, 3)
    for a in 1:3
        N[1, a, a] = 1
        N[a, 1, a] = 1
    end
    N[2, 2, 1] = 1
    N[2, 2, 3] = 1
    N[2, 3, 2] = 1
    N[3, 2, 2] = 1
    N[3, 3, 1] = 1
    return N
end
