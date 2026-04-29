"""
Galois actions on cyclotomic elements and exact modular data.

The maps are conductor-local: exponents are interpreted using the
CyclotomicContext carried by the input data.
"""

function _galois_elem(ctx::CyclotomicContext, x, a::Int)
    gcd(a, ctx.N) == 1 || error("Galois exponent $a is not a unit modulo $(ctx.N)")
    coeffs = coordinates(x)
    za = ctx.zeta^mod(a, ctx.N)
    y = zero(ctx.field)
    for k in eachindex(coeffs)
        y += ctx.field(Rational{BigInt}(coeffs[k])) * za^(k - 1)
    end
    return y
end

function galois_action(ctx::CyclotomicContext, x, a::Int)
    return _galois_elem(ctx, x, a)
end

function galois_action(ctx::CyclotomicContext, M::MatElem, a::Int)
    K = ctx.field
    out = zero_matrix(K, nrows(M), ncols(M))
    for i in 1:nrows(M), j in 1:ncols(M)
        out[i, j] = galois_action(ctx, M[i, j], a)
    end
    return out
end

function galois_action(ctx::CyclotomicContext, v::AbstractVector, a::Int)
    return [galois_action(ctx, x, a) for x in v]
end

"""
    galois_action(ctx, S, T, a) -> NamedTuple

Apply `σ_a(ζ_N) = ζ_N^a` entrywise to an exact `S,T` pair over the
cyclotomic context `ctx`.
"""
function galois_action(ctx::CyclotomicContext, S, T, a::Int)
    return (S = galois_action(ctx, S, a),
            T = galois_action(ctx, T, a))
end

function galois_action(S, T, a::Int; context = nothing, conductor = nothing, N = nothing)
    ctx = if context !== nothing
        context
    else
        n = conductor === nothing ? N : conductor
        n === nothing && error("a CyclotomicContext or conductor N is required")
        CyclotomicContext(n)
    end
    return galois_action(ctx, S, T, a)
end

function galois_action(data::ModularData, a::Int)
    ctx = data.context
    return ModularData(ctx, copy(data.labels),
                       galois_action(ctx, data.S, a),
                       galois_action(ctx, data.T, a),
                       data.cond_S, data.cond_T, data.cond_F)
end

_galois_units(N::Int) = [a for a in 1:N if gcd(a, N) == 1]

galois_orbit(data::ModularData) =
    [galois_action(data, a) for a in _galois_units(data.context.N)]

function _galois_dim(M)
    return M isa MatElem ? (nrows(M), ncols(M)) : size(M)
end

function _galois_twists(T)
    try
        d = _galois_dim(T)
        length(d) == 2 && d[1] == d[2] && return [T[i, i] for i in 1:d[1]]
    catch
    end
    return collect(T)
end

function _permutes_ST(S, T, perm::AbstractVector{<:Integer})
    r, c = _galois_dim(S)
    r == c == length(perm) || return false
    collect(sort(perm)) == collect(1:r) || return false
    perm[1] == 1 || return false
    twists = _galois_twists(T)
    length(twists) == r || return false
    @inbounds for i in 1:r
        twists[perm[i]] == twists[i] || return false
    end
    @inbounds for i in 1:r, j in 1:r
        S[perm[i], perm[j]] == S[i, j] || return false
    end
    return true
end

"""
    modular_data_automorphisms(data_or_S, T = nothing) -> Vector{Vector{Int}}

Return all unit-fixing anyon relabelings that preserve exact `S` and `T`.
Permutations are 1-indexed and satisfy `S[perm, perm] == S` and
`T[perm] == T`.
"""
function modular_data_automorphisms(S, T)
    r, c = _galois_dim(S)
    r == c || error("S must be square")
    length(_galois_twists(T)) == r || error("T has wrong length")
    r == 1 && return [Int[1]]
    autos = Vector{Vector{Int}}()
    for rest in _permutations(2:r)
        perm = [1; rest...]
        _permutes_ST(S, T, perm) && push!(autos, perm)
    end
    return autos
end

modular_data_automorphisms(data::ModularData) =
    modular_data_automorphisms(data.S, data.T)

is_modular_data_automorphism(S, T, perm::AbstractVector{<:Integer}) =
    _permutes_ST(S, T, perm)

is_modular_data_automorphism(data::ModularData, perm::AbstractVector{<:Integer}) =
    is_modular_data_automorphism(data.S, data.T, perm)

function _galois_anyon_action(data::ModularData, a::Int)
    ctx = data.context
    Sσ = galois_action(ctx, data.S, a)
    r, c = _galois_dim(data.S)
    r == c || error("S must be square")
    perm = zeros(Int, r)
    scalars = Vector{Any}(undef, r)
    used = falses(r)
    for i in 1:r
        found = 0
        for j in 1:r
            used[j] && continue
            scalar = _projective_row_scalar(Sσ, i, data.S, j)
            if scalar !== nothing
                found = j
                scalars[i] = scalar
                break
            end
        end
        found == 0 && error("could not match Galois-transformed anyon row $i for exponent $a")
        perm[i] = found
        used[found] = true
    end
    return (exponent = a, perm = perm, scalars = scalars)
end

function _projective_row_scalar(A, i::Int, B, j::Int)
    _, n = _galois_dim(A)
    scalar = nothing
    for k in 1:n
        b = B[j, k]
        a = A[i, k]
        if iszero(b)
            iszero(a) || return nothing
        else
            cand = a / b
            scalar === nothing ? (scalar = cand) : (scalar == cand || return nothing)
        end
    end
    return scalar === nothing ? 1 : scalar
end

"""
    galois_anyon_action(data, a) -> NamedTuple

Return the projective anyon permutation induced by the Galois exponent `a`.
The result has fields `exponent`, `perm`, and `scalars`; `scalars[i]`
records the projective row factor relating `σ_a(S[i, :])` to
`S[perm[i], :]`.
"""
galois_anyon_action(data::ModularData, a::Int) = _galois_anyon_action(data, a)

"""
    galois_anyon_orbits(data) -> Vector{Vector{Int}}

Return the orbits of simple objects under all Galois-induced anyon
permutations for `data`.
"""
function galois_anyon_orbits(data::ModularData)
    actions = [galois_anyon_action(data, a).perm for a in _galois_units(data.context.N)]
    r = length(data.labels)
    seen = falses(r)
    orbits = Vector{Vector{Int}}()
    for i in 1:r
        seen[i] && continue
        orbit = sort(unique([perm[i] for perm in actions]))
        for j in orbit
            seen[j] = true
        end
        push!(orbits, orbit)
    end
    return orbits
end
