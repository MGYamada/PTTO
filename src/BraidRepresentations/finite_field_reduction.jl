struct FiniteFieldBraidRepresentation
    source::Any
    p::Int
    objects::Vector{Int}
    total::Int
    basis::FusionTreeBasis
    generators::Vector{Matrix{Int}}
    metadata::Dict{Symbol, Any}
end

function _reduce_scalar_mod_p(x, p::Int, zeta_Fp)
    x isa Integer && return mod(Int(x), p)
    x isa Rational && return mod(Int(numerator(x)) * invmod(Int(denominator(x)), p), p)
    zeta_Fp === nothing && error("cannot reduce $(typeof(x)) modulo $p without a conductor/root-of-unity hint")
    return cyclotomic_to_Fp(x, zeta_Fp, p)
end

function _reduce_matrix_mod_p(A, p::Int, zeta_Fp)
    M = zeros(Int, size(A, 1), size(A, 2))
    for i in axes(A, 1), j in axes(A, 2)
        M[i, j] = _reduce_scalar_mod_p(A[i, j], p, zeta_Fp)
    end
    return M
end

function reduce_mod_p(br::BraidRepresentation, p::Integer; conductor = nothing)
    p = Int(p)
    isprime(p) || error("reduce_mod_p for braid representations requires a prime p, got $p")
    zeta_Fp = nothing
    if conductor !== nothing
        (p - 1) % Int(conductor) == 0 ||
            error("F_$p does not contain a primitive $(conductor)-th root of unity; finite-field extensions are not implemented yet")
        zeta_Fp = find_zeta_in_Fp(Int(conductor), p)
    end
    gens = [_reduce_matrix_mod_p(g, p, zeta_Fp) for g in br.generators]
    return FiniteFieldBraidRepresentation(br, p, br.objects, br.total, br.basis, gens,
        Dict(:p => p, :conductor => conductor, :p_mod_conductor => conductor === nothing ? nothing : mod(p, Int(conductor)),
             :extension_field_required => false))
end

check_braid_relations(br::FiniteFieldBraidRepresentation; kwargs...) =
    check_braid_relations(br.generators, br.p)
