struct BraidRepresentation{T, M<:AbstractMatrix{T}}
    fr_data::FRData{T}
    objects::Vector{Int}
    total::Int
    basis::FusionTreeBasis
    generators::Vector{M}
    metadata::AbstractDict{Symbol}
end

function _flookup(data::FRData{T}, a,b,c,d,e,f) where {T}
    return F_symbol(data, a, b, c, d; e = e, f = f)
end

function _rlookup(data::FRData{T}, a,b,c) where {T}
    return R_symbol(data, a, b, c)
end

function _matrix_inverse_generic(A::AbstractMatrix{T}) where {T}
    n, m = size(A)
    n == m || error("cannot invert non-square F-move block of size $(size(A))")
    M = T[A[i, j] for i in 1:n, j in 1:n]
    I0 = T[i == j ? one(A[1, 1]) : zero(A[1, 1]) for i in 1:n, j in 1:n]
    for col in 1:n
        piv = findfirst(r -> !iszero(M[r, col]), col:n)
        piv === nothing && error("singular F-move block in braid generator construction")
        piv = col + piv - 1
        if piv != col
            M[col, :], M[piv, :] = M[piv, :], M[col, :]
            I0[col, :], I0[piv, :] = I0[piv, :], I0[col, :]
        end
        invp = inv(M[col, col])
        M[col, :] = M[col, :] .* invp
        I0[col, :] = I0[col, :] .* invp
        for r in 1:n
            r == col && continue
            factor = M[r, col]
            iszero(factor) && continue
            M[r, :] = M[r, :] .- factor .* M[col, :]
            I0[r, :] = I0[r, :] .- factor .* I0[col, :]
        end
    end
    return I0
end

function _local_braid_coeff(data::FRData{T}, y::Int, a::Int, b::Int, w::Int,
                            z_in::Int, z_out::Int)::T where {T}
    fr = fusion_rule(data)
    us = [u for u in fusion_channels(fr, a, b) if is_admissible(fr, y, u, w)]
    zs = [z for z in fusion_channels(fr, y, a) if is_admissible(fr, z, b, w)]
    sort!(us); sort!(zs)
    F = T[_flookup(data, y, a, b, w, z, u) for z in zs, u in us]
    Finv = _matrix_inverse_generic(F)
    iout = findfirst(==(z_out), zs)
    iin = findfirst(==(z_in), zs)
    iout === nothing && return zero(fr_value_one(data))
    iin === nothing && return zero(fr_value_one(data))
    s = zero(F[1,1])
    for (uix, u) in enumerate(us)
        s += F[iout, uix] * _rlookup(data, a, b, u) * Finv[uix, iin]
    end
    return s
end

function _matrix_for_generator(data::FRData{T}, basis::FusionTreeBasis, i::Int)::Matrix{T} where {T}
    nobj = length(basis.objects)
    1 <= i < nobj || error("braid generator index i=$i outside 1:$(nobj - 1)")
    d = dim(basis)
    M = fill(zero(fr_value_one(data)), d, d)
    index = Dict(Tuple(p.channels) => k for (k, p) in enumerate(basis.paths))
    for (col, path) in enumerate(basis.paths)
        ch = path.channels
        if i == 1
            out = copy(ch)
            row = index[Tuple(out)]
            M[row, col] = _rlookup(data, basis.objects[1], basis.objects[2], ch[1])
        else
            y = i == 2 ? basis.objects[1] : ch[i - 2]
            a = basis.objects[i]
            b = basis.objects[i + 1]
            w = ch[i]
            z_in = ch[i - 1]
            for z_out in fusion_channels(fusion_rule(data), y, a)
                is_admissible(fusion_rule(data), z_out, b, w) || continue
                out = copy(ch)
                out[i - 1] = z_out
                row = get(index, Tuple(out), 0)
                row == 0 && continue
                M[row, col] = _local_braid_coeff(data, y, a, b, w, z_in, z_out)
            end
        end
    end
    return M
end

function braid_representation(fr_data::FRData{T}, objects::AbstractVector, total) where {T}
    basis = fusion_basis(fusion_rule(fr_data), objects, total)
    gens = Matrix{T}[_matrix_for_generator(fr_data, basis, i)
                     for i in 1:(length(basis.objects) - 1)]
    metadata = Dict{Symbol, Any}(:multiplicity_free => true)
    return BraidRepresentation{T, Matrix{T}}(fr_data, basis.objects, basis.total, basis,
                                             gens, metadata)
end

braid_generator(br::BraidRepresentation, i::Int) = br.generators[i]
braid_generators(br::BraidRepresentation) = br.generators
braid_generator(fr_data::FRData, objects::AbstractVector, total, i::Int) =
    braid_generator(braid_representation(fr_data, objects, total), i)
braid_generators(fr_data::FRData, objects::AbstractVector, total) =
    braid_generators(braid_representation(fr_data, objects, total))
braid_generator(fr_data::FRData, objects::AbstractVector, i::Int; total_charge) =
    braid_generator(fr_data, objects, total_charge, i)

function _identity_like(A::AbstractMatrix{T}) where {T}
    n, m = size(A)
    n == m || error("identity requested for non-square matrix of size $(size(A))")
    oneval = one(A[1, 1])
    return T[i == j ? oneval : zero(oneval) for i in 1:n, j in 1:n]
end

function _matrix_power_generic(A::AbstractMatrix, n::Integer)
    n == 0 && return _identity_like(A)
    n < 0 && return _matrix_power_generic(_matrix_inverse_generic(A), -n)
    out = _identity_like(A)
    base = A
    k = Int(n)
    while k > 0
        isodd(k) && (out = out * base)
        k ÷= 2
        k > 0 && (base = base * base)
    end
    return out
end

function braid_representation(fr_data::FRData, braid_word::AbstractVector{<:Integer},
                              objects::AbstractVector; total_charge)
    gens = braid_generators(fr_data, objects, total_charge)
    isempty(gens) && return Matrix{fr_scalar_type(fr_data)}(undef, 0, 0)
    out = _identity_like(first(gens))
    for letter in braid_word
        idx = abs(Int(letter))
        1 <= idx <= length(gens) || error("braid generator σ_$idx is outside 1:$(length(gens))")
        out = out * _matrix_power_generic(gens[idx], letter < 0 ? -1 : 1)
    end
    return out
end

braid_generators_B3(fr_data::FRData, objects::AbstractVector; total_charge) =
    Tuple(braid_generators(fr_data, objects, total_charge))

function braid_representation(fr_tuple::NamedTuple, objects::AbstractVector, total)
    return braid_representation(frdata_from_namedtuple(fr_tuple), objects, total)
end

braid_representation(rules, Fvals::AbstractVector, Rvals::AbstractVector,
                     objects::AbstractVector, total) =
    braid_representation(frdata_from_vectors(rules, Fvals, Rvals), objects, total)
