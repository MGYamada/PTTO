"""
    MTCCandidate

Result of a successful block-U sweep at a single prime.

Fields:
- p: the prime used
- U_params: parameters of the block-U used (e.g. (u, v, det) for O(2))
- S_Fp: the S matrix mod p after block-U transformation
- T_Fp: the T eigenvalues mod p (unchanged by block-U)
- unit_index: which basis index is the unit object
- N: fusion tensor N[i][j][k] (signed integers)
- d: quantum dimensions d_i = S[unit, i] / S[unit, unit]
- D2: total quantum dimension squared
"""
struct MTCCandidate
    p::Int
    U_params::Any
    S_Fp::Matrix{Int}
    T_Fp::Vector{Int}
    unit_index::Int
    N::Array{Int, 3}
    d::Vector{Int}
    D2::Int
end

function Base.show(io::IO, c::MTCCandidate)
    d_signed = [signed_Fp(x, c.p) for x in c.d]
    # U_params can be a Matrix (general O(n)) or a Tuple (legacy O(2))
    params_str = if isa(c.U_params, AbstractMatrix)
        n = size(c.U_params, 1)
        "U=$(n)×$(n) block"
    elseif isa(c.U_params, AbstractVector) && all(x -> x isa NamedTuple && haskey(x, :U), c.U_params)
        dims = [size(x.U, 1) for x in c.U_params]
        "U blocks=$(join(dims, ","))"
    else
        "params=$(c.U_params)"
    end
    print(io, "MTCCandidate(p=$(c.p), unit=$(c.unit_index), ",
          "d=$d_signed, D²=$(signed_Fp(c.D2, c.p)), ",
          "$params_str)")
end
