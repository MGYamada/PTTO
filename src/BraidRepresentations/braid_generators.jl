struct FRData
    rules::FusionRule
    F_values::Vector
    R_values::Vector
    R_inverse_values::Vector
    metadata::Dict{Symbol, Any}
end

function FRData(rules::FusionRule,
                F_values::AbstractVector,
                R_values::AbstractVector,
                R_inverse_values::AbstractVector,
                metadata::Dict{Symbol, Any} = Dict{Symbol, Any}())
    return _frdata_from_vectors(rules, F_values, vcat(collect(R_values), collect(R_inverse_values));
                                metadata = metadata)
end

struct BraidRepresentation
    fr_data::FRData
    objects::Vector{Int}
    total::Int
    basis::FusionTreeBasis
    generators::Vector{Any}
    metadata::Dict{Symbol, Any}
end

function _fr_value_one(data::FRData)
    !isempty(data.F_values) && return one(data.F_values[1])
    !isempty(data.R_values) && return one(data.R_values[1])
    return 1
end

function _flookup(data::FRData, a,b,c,d,e,f)
    if a == 1 || b == 1 || c == 1
        return e == f ? _fr_value_one(data) : zero(_fr_value_one(data))
    end
    idx = _f_value_index(data.rules, a, b, c, d, e, f)
    idx === nothing &&
        error("missing F-symbol F^($a,$b,$c)_$d[$e,$f] in multiplicity-free braid construction")
    return data.F_values[idx]
end

function _rlookup(data::FRData, a,b,c)
    if a == 1 || b == 1
        return _fr_value_one(data)
    end
    idx = _r_value_index(data.rules, a, b, c)
    idx === nothing &&
        error("missing R-symbol R^($a,$b)_$c in multiplicity-free braid construction")
    return data.R_values[idx]
end

function _matrix_inverse_generic(A)
    n, m = size(A)
    n == m || error("cannot invert non-square F-move block of size $(size(A))")
    M = [A[i, j] for i in 1:n, j in 1:n]
    I0 = [i == j ? one(A[1, 1]) : zero(A[1, 1]) for i in 1:n, j in 1:n]
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

function _local_braid_coeff(data::FRData, y::Int, a::Int, b::Int, w::Int,
                            z_in::Int, z_out::Int)
    fr = data.rules
    us = [u for u in fusion_channels(fr, a, b) if is_admissible(fr, y, u, w)]
    zs = [z for z in fusion_channels(fr, y, a) if is_admissible(fr, z, b, w)]
    sort!(us); sort!(zs)
    F = [_flookup(data, y, a, b, w, z, u) for z in zs, u in us]
    Finv = _matrix_inverse_generic(F)
    iout = findfirst(==(z_out), zs)
    iin = findfirst(==(z_in), zs)
    iout === nothing && return zero(_fr_value_one(data))
    iin === nothing && return zero(_fr_value_one(data))
    s = zero(F[1,1])
    for (uix, u) in enumerate(us)
        s += F[iout, uix] * _rlookup(data, a, b, u) * Finv[uix, iin]
    end
    return s
end

function _matrix_for_generator(data::FRData, basis::FusionTreeBasis, i::Int)
    nobj = length(basis.objects)
    1 <= i < nobj || error("braid generator index i=$i outside 1:$(nobj - 1)")
    d = dim(basis)
    T = typeof(_fr_value_one(data))
    M = [zero(_fr_value_one(data)) for _ in 1:d, _ in 1:d]
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
            for z_out in fusion_channels(data.rules, y, a)
                is_admissible(data.rules, z_out, b, w) || continue
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

function braid_representation(fr_data::FRData, objects::AbstractVector, total)
    basis = fusion_basis(fr_data.rules, objects, total)
    gens = [_matrix_for_generator(fr_data, basis, i) for i in 1:(length(basis.objects) - 1)]
    return BraidRepresentation(fr_data, basis.objects, basis.total, basis, gens,
                               Dict(:multiplicity_free => true))
end

braid_generator(br::BraidRepresentation, i::Int) = br.generators[i]
braid_generators(br::BraidRepresentation) = br.generators
braid_generator(fr_data::FRData, objects::AbstractVector, total, i::Int) =
    braid_generator(braid_representation(fr_data, objects, total), i)
braid_generators(fr_data::FRData, objects::AbstractVector, total) =
    braid_generators(braid_representation(fr_data, objects, total))

function _frdata_from_vectors(rules, Fvals::AbstractVector, Rvals::AbstractVector;
                              metadata = Dict{Symbol, Any}())
    fr = _fusion_rule(rules)
    require_multiplicity_free(fr)
    _, _, nF = get_pentagon_system(fr.N, fr.rank)
    length(Fvals) == nF ||
        error("F_values has length $(length(Fvals)); expected $nF in TensorCategories pentagon order")
    positions, nforward = _braiding_block_positions(fr.N)
    length(Rvals) >= nforward || error("R vector has length $(length(Rvals)); expected at least $nforward")
    length(Rvals) == nforward || length(Rvals) == 2 * nforward ||
        error("R_values has length $(length(Rvals)); expected $nforward or $(2 * nforward)")
    invvals = length(Rvals) >= 2 * nforward ? collect(Rvals[(nforward + 1):(2 * nforward)]) :
              _inverse_R_values_from_forward_values(fr, collect(Rvals[1:nforward]))
    return FRData(fr, collect(Fvals), collect(Rvals[1:nforward]), invvals, metadata)
end

function braid_representation(fr_tuple::NamedTuple, objects::AbstractVector, total)
    if hasproperty(fr_tuple, :rules) && hasproperty(fr_tuple, :F_values) && hasproperty(fr_tuple, :R_values)
        Rinv = hasproperty(fr_tuple, :R_inverse_values) ? fr_tuple.R_inverse_values : Any[]
        return braid_representation(FRData(_fusion_rule(fr_tuple.rules), fr_tuple.F_values,
                                           fr_tuple.R_values, Rinv,
                                           Dict{Symbol, Any}(:source => :namedtuple)),
                                    objects, total)
    end
    error("FR data NamedTuple must contain :rules, :F_values, and :R_values in TensorCategories order")
end

braid_representation(rules, Fvals::AbstractVector, Rvals::AbstractVector,
                     objects::AbstractVector, total) =
    braid_representation(_frdata_from_vectors(rules, Fvals, Rvals), objects, total)

function _inverse_R_values_from_forward_values(fr::FusionRule, R_values::Vector)
    positions, nR = _braiding_block_positions(fr.N)
    length(R_values) == nR || error("R_values has length $(length(R_values)); expected $nR")
    values = Vector{Any}(undef, nR)
    for ((a, b, c), pos) in positions
        opp = positions[(b, a, c)]
        values[first(pos)] = inv(R_values[first(opp)])
    end
    return values
end

function _f_value_index(fr::FusionRule, a::Int, b::Int, c::Int, d::Int, e::Int, f::Int)
    _, _, nF = get_pentagon_system(fr.N, fr.rank)
    nF == 0 && return nothing
    for (idx, m) in enumerate(_pentagon_variable_metadata(fr.N, fr.rank, nF))
        (m.i, m.j, m.k, m.o, m.a, m.b) == (a, b, c, d, e, f) && return idx
    end
    return nothing
end

function _r_value_index(fr::FusionRule, a::Int, b::Int, c::Int)
    positions, _ = _braiding_block_positions(fr.N)
    pos = get(positions, (a, b, c), nothing)
    pos === nothing && return nothing
    length(pos) == 1 || error("multiplicity > 1 R-blocks are not supported by braid representations")
    return first(pos)
end

_fr_hexagon_values(data::FRData) = vcat(data.F_values, data.R_values, data.R_inverse_values)
_fr_pentagon_values(data::FRData) = data.F_values

function semion_fr_data(ctx::CyclotomicContext = CyclotomicContext(8))
    data = semion_modular_data(ctx)
    K, z = field(ctx), zeta(ctx)
    F = Any[-one(K)]
    R = Any[z^2, one(K), one(K), one(K)]
    Rinv = Any[-z^2, one(K), one(K), one(K)]
    return _frdata_from_vectors(semion_fusion_rules(), F, vcat(R, Rinv);
        metadata = Dict{Symbol, Any}(:name => :semion, :conductor => 8,
                                     :modular_data => data,
                                     :format => :tensorcategories_variable_order,
                                     :source => :hardcoded_verified_solution))
end

function fibonacci_fr_data(ctx::CyclotomicContext = CyclotomicContext(20))
    data = fibonacci_modular_data(ctx)
    K, z = field(ctx), zeta(ctx)
    F = Any[z^6 - z^4,
            one(K),
            -z^6 + z^4,
            -z^6 + z^4,
            one(K)]
    R = Any[-z^4,
            z^6 - z^4 + z^2 - one(K),
            one(K),
            one(K),
            one(K)]
    Rinv = Any[z^6,
               -z^2,
               one(K),
               one(K),
               one(K)]
    return _frdata_from_vectors(fibonacci_fusion_rules(), F, vcat(R, Rinv);
        metadata = Dict{Symbol, Any}(:name => :fibonacci, :conductor => 20,
                                     :modular_data => data,
                                     :format => :tensorcategories_variable_order,
                                     :source => :hardcoded_verified_solution))
end

function ising_fr_data(ctx::CyclotomicContext = CyclotomicContext(16))
    data = ising_modular_data(ctx)
    K, z = field(ctx), zeta(ctx)
    s = (-z^6 + z^2) // K(2)
    F = Any[one(K),
            one(K),
            -one(K),
            one(K),
            one(K),
            one(K),
            one(K),
            -one(K),
            -one(K),
            -one(K),
            s,
            -s,
            s,
            s]
    R = Any[-one(K),
            z^4,
            one(K),
            z^4,
            -z^3,
            z^7,
            one(K),
            one(K),
            one(K),
            one(K)]
    Rinv = Any[-one(K),
               -z^4,
               one(K),
               -z^4,
               z^5,
               -z,
               one(K),
               one(K),
               one(K),
               one(K)]
    return _frdata_from_vectors(ising_fusion_rules(), F, vcat(R, Rinv);
        metadata = Dict{Symbol, Any}(:name => :ising, :conductor => 16,
                                     :modular_data => data,
                                     :format => :tensorcategories_variable_order,
                                     :source => :hardcoded_verified_solution))
end
