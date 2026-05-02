"""
Gauge helpers for exact multiplicity-free `(F, R)` data.

For a fusion channel basis vector `V_ab^c`, a gauge scalar `u[a,b,c]`
acts by

    F'^{ijk}_{o;a,b} = u[i,j,a] u[a,k,o] / (u[j,k,b] u[i,b,o]) F^{ijk}_{o;a,b}
    R'^{ij}_k       = u[i,j,k] / u[j,i,k] R^{ij}_k

The current implementation uses the full channel-scalar toric gauge group:
unit channels such as `(1,a,a)`, `(a,1,a)`, and `(1,1,1)` are included as
gauge parameters.  Consequently the action usually has an ineffective kernel.
This convention is natural for stacky stabilizer counts, while unit-normalized
or effective toric quotients should be treated as separate conventions.
The implementation is intentionally exact and multiplicity-free; higher
multiplicity needs full change-of-basis matrices rather than channel scalars.
"""

function _identity_gauge(Nijk::Array{Int,3})
    r = size(Nijk, 1)
    return Dict{Tuple{Int, Int, Int}, Int}(
        (a, b, c) => 1 for a in 1:r, b in 1:r, c in 1:r if Nijk[a, b, c] != 0
    )
end

_identity_gauge(data::FRData) =
    Dict((a, b, c) => 1 for (a, b, c) in gauge_parameters(data))

"""
    validate_frdata_for_gauge(frdata)

Check the FRData invariants required by ACMG's scalar gauge layer.  This
belongs to Gauge, not FRData, because multiplicityful FRData is valid even
though scalar channel gauges are currently multiplicity-free only.
"""
function validate_frdata_for_gauge(data::FRData)
    validate(fusion_rule(data))
    is_multiplicity_free(fusion_rule(data)) ||
        error("ACMG scalar gauge transforms currently require multiplicity-free FRData")
    for a in simples(data), b in simples(data), c in fusion_channels(data, a, b)
        isempty(hom_basis(data, a, b, c)) &&
            error("admissible fusion channel ($a,$b,$c) has empty Hom basis")
    end
    return true
end

function _one_like_from_FR(F, R)
    !isempty(F) && return one(F[1])
    !isempty(R) && return one(R[1])
    return 1
end

function _one_gauge_scalars(fusion_rule, one_value)
    fr = _gauge_rules(fusion_rule)
    scalars = Dict{Tuple{Int,Int,Int}, typeof(one_value)}()
    for a in 1:fr.rank, b in 1:fr.rank, c in 1:fr.rank
        fusion_coeff(fr, a, b, c) == 0 && continue
        scalars[(a, b, c)] = one_value
    end
    return scalars
end

function _require_multiplicity_free_gauge(fusion_rule)
    is_multiplicity_free(_gauge_rules(fusion_rule)) ||
        error("channel-scalar gauge transforms require multiplicity-free fusion rules")
    return true
end

function _gauge_parameters_as_scalars(gauge::GaugeParameters, default_one)
    out = Dict{Tuple{Int,Int,Int}, typeof(default_one)}()
    for ((a, b, c, μ), value) in gauge.values
        μ == 1 || error("channel-scalar gauge transforms currently require μ = 1")
        out[(a, b, c)] = value
    end
    return out
end

function _gauge_scalar(gauge, ch::Tuple{Int,Int,Int}, default_one)
    if gauge isa GaugeTransform
        return get(gauge.scalars, ch, default_one)
    elseif gauge isa GaugeParameters
        return get(gauge.values, (ch[1], ch[2], ch[3], 1), default_one)
    elseif gauge isa GaugeAction
        return get(gauge.parameters, (ch[1], ch[2], ch[3], 1), default_one)
    elseif gauge isa AbstractDict
        return get(gauge, ch, default_one)
    elseif gauge === nothing
        return default_one
    end
    error("unsupported gauge object: $(typeof(gauge))")
end

function _pow_int_exact(x, n::Int)
    n == 0 && return one(x)
    n > 0 && return x^n
    return inv(x)^(-n)
end

function _factor_from_weight(weight::Vector{Int},
                             channels::Vector{Tuple{Int,Int,Int}},
                             scalars,
                             one_value)
    factor = one_value
    for (idx, exponent) in enumerate(weight)
        exponent == 0 && continue
        factor *= _pow_int_exact(_gauge_scalar(scalars, channels[idx], one_value), exponent)
    end
    return factor
end

function _ordered_gauge_channels(fusion_rule, r::Int)
    fr = _gauge_rules(fusion_rule)
    indexed = collect(_gauge_channel_basis(fr.N, r))
    sort!(indexed; by = kv -> kv.second)
    return [kv.first for kv in indexed]
end

function _r_var_count(fusion_rule)
    fr = _gauge_rules(fusion_rule)
    return sum(fusion_coeff(fr, i, j, k)^2 for i in 1:fr.rank, j in 1:fr.rank, k in 1:fr.rank)
end

function _r_forward_transform_factor(Nijk::Array{Int,3},
                                     i::Int, j::Int, k::Int,
                                     gauge,
                                     one_value)
    numerator = _gauge_scalar(gauge, (i, j, k), one_value)
    denominator = _gauge_scalar(gauge, (j, i, k), one_value)
    return numerator / denominator
end

function _r_inverse_transform_factor(Nijk::Array{Int,3},
                                     i::Int, j::Int, k::Int,
                                     gauge,
                                     one_value)
    return _r_forward_transform_factor(Nijk, i, j, k, gauge, one_value)
end

function _transform_F_values(F, fusion_rule, gauge, one_value)
    fr = _gauge_rules(fusion_rule)
    r = fr.rank
    metadata = _pentagon_variable_metadata(fr.N, r, length(F))
    channel_index = _gauge_channel_basis(fr.N, r)
    channels = _ordered_gauge_channels(fr, r)
    out = copy(F)
    scalar_gauge = gauge isa GaugeParameters ? _gauge_parameters_as_scalars(gauge, one_value) :
                   gauge isa GaugeTransform ? gauge.scalars : gauge
    for idx in eachindex(F)
        weight = _gauge_weight(metadata[idx], channel_index)
        out[idx] = F[idx] * _factor_from_weight(weight, channels, scalar_gauge, one_value)
    end
    return out
end

function _transform_R_values(R, fusion_rule, gauge, one_value)
    fr = _gauge_rules(fusion_rule)
    r = fr.rank
    n_forward = _r_var_count(fr)
    (length(R) == n_forward || length(R) == 2n_forward) ||
        error("R has length $(length(R)); expected $n_forward or $(2n_forward)")
    positions, _ = _braiding_block_positions(fr.N)
    out = copy(R)
    for i in 1:r, j in 1:r, k in 1:r
        fusion_coeff(fr, i, j, k) == 0 && continue
        fwd = _r_forward_transform_factor(fr.N, i, j, k, gauge, one_value)
        for pos in positions[(i, j, k)]
            out[pos] = R[pos] * fwd
        end
        if length(R) == 2n_forward
            inv_factor = _r_inverse_transform_factor(fr.N, i, j, k, gauge, one_value)
            for pos in positions[(i, j, k)]
                out[n_forward + pos] = R[n_forward + pos] * inv_factor
            end
        end
    end
    return out
end

function _gauge_fixing_plan(F, fusion_rule)
    fr = _gauge_rules(fusion_rule)
    r = fr.rank
    _require_multiplicity_free_gauge(fr)
    metadata = _pentagon_variable_metadata(fr.N, r, length(F))
    channel_index = _gauge_channel_basis(fr.N, r)
    channels = _ordered_gauge_channels(fr, r)
    isempty(channels) && return NamedTuple[]

    used_pivots = Int[]
    plan = NamedTuple[]
    for var_idx in eachindex(F)
        iszero(F[var_idx]) && continue
        weight = _gauge_weight(metadata[var_idx], channel_index)
        any(!iszero, weight) || continue
        all(iszero(weight[p]) for p in used_pivots) || continue
        pivot = findfirst(p -> !(p in used_pivots) && abs(weight[p]) == 1,
                          eachindex(weight))
        pivot === nothing && continue
        push!(plan, (var_idx = var_idx,
                     weight = weight,
                     pivot = pivot,
                     pivot_channel = channels[pivot],
                     pivot_exponent = weight[pivot]))
        push!(used_pivots, pivot)
    end
    return plan
end

"""
    gauge_fixing_plan(F, Nijk)

Return the independent F-symbol coordinates used by `canonical_gauge` to
fix channel-scalar gauge freedom.  Each entry records the F index and the
pivot fusion channel whose scalar is solved exactly.
"""
gauge_fixing_plan(F, Nijk::Array{Int,3}) = _gauge_fixing_plan(F, FusionRule(Nijk))
gauge_fixing_plan(data::FRData) = _gauge_fixing_plan(F_values(data), data)

function _canonical_gauge_transform(F, R, fusion_rule)
    fr = _gauge_rules(fusion_rule)
    _require_multiplicity_free_gauge(fr)
    one_value = _one_like_from_FR(F, R)
    scalars = _one_gauge_scalars(fr, one_value)
    channels = _ordered_gauge_channels(fr, fr.rank)
    plan = _gauge_fixing_plan(F, fr)

    for entry in reverse(plan)
        rest = one_value
        for (idx, exponent) in enumerate(entry.weight)
            (idx == entry.pivot || exponent == 0) && continue
            rest *= _pow_int_exact(get(scalars, channels[idx], one_value), exponent)
        end
        target = inv(F[entry.var_idx] * rest)
        scalars[entry.pivot_channel] =
            entry.pivot_exponent == 1 ? target : inv(target)
    end

    fixed = Int[entry.var_idx for entry in plan]
    complete = length(fixed) == length(_select_pentagon_gauge_fixed_indices(fr.N, fr.rank, length(F)))
    return GaugeTransform(scalars, fixed, complete)
end

"""
    canonical_gauge(F, R, Nijk)

Return a canonical channel-scalar gauge representative for multiplicity-free
exact `(F, R)` data.  Independent nonzero F-symbol coordinates are fixed to
`1` whenever this can be done with exact multiplication and inversion.
"""
function canonical_gauge(F, R, Nijk::Array{Int,3})
    fr = FusionRule(Nijk)
    gauge = _canonical_gauge_transform(F, R, fr)
    transformed = gauge_transform(F, R, gauge; Nijk = Nijk)
    return (F = transformed.F, R = transformed.R, gauge = gauge)
end

function canonical_gauge(data::FRData)
    validate_frdata_for_gauge(data)
    gauge = _canonical_gauge_transform(F_values(data), R_values(data), data)
    transformed = gauge_transform(data, gauge)
    return GaugeFixingResult(F_values(transformed), R_values(transformed), gauge,
                             gauge.fixed_indices, gauge.complete,
                             Dict{Symbol, Any}(:frdata => transformed))
end

"""
    gauge_transform(F, R, gauge; Nijk)

Apply a channel-scalar gauge transform to `(F, R)`.  `gauge` may be a
`GaugeTransform`, a dictionary keyed by `(a,b,c)`, or `nothing`.
"""
function gauge_transform(F, R, gauge; Nijk::Union{Array{Int,3}, Nothing} = nothing)
    if Nijk === nothing
        if gauge isa AbstractDict
            all(v -> v == 1, values(gauge)) ||
                error("Nijk is required for nontrivial gauge transforms")
            return (F = copy(F), R = copy(R))
        elseif gauge isa GaugeTransform
            all(v -> v == 1, values(gauge.scalars)) ||
                error("Nijk is required for nontrivial gauge transforms")
            return (F = copy(F), R = copy(R))
        elseif gauge === nothing
            return (F = copy(F), R = copy(R))
        end
        error("unsupported gauge object: $(typeof(gauge))")
    end
    fr = FusionRule(Nijk)
    _require_multiplicity_free_gauge(fr)
    one_value = _one_like_from_FR(F, R)
    return (F = _transform_F_values(F, fr, gauge, one_value),
            R = _transform_R_values(R, fr, gauge, one_value))
end

function gauge_transform(data::FRData, gauge)
    validate_frdata_for_gauge(data)
    one_value = fr_value_one(data)
    R_all = vcat(R_values(data), R_inverse_values(data))
    Fout = _transform_F_values(F_values(data), data, gauge, one_value)
    Rout = _transform_R_values(R_all, data, gauge, one_value)
    metadata = copy(fr_metadata(data))
    metadata[:gauge_transform] = gauge
    return frdata_from_vectors(fusion_rule(data), Fout, Rout; metadata = metadata)
end

"""
    is_gauge_fixed(F, Nijk)

Return `true` when all coordinates selected by `gauge_fixing_plan(F, Nijk)`
are exactly `1`.
"""
function is_gauge_fixed(F, Nijk::Array{Int,3})
    one_value = _one_like_from_FR(F, [])
    return all(idx -> F[idx] == one_value, [entry.var_idx for entry in _gauge_fixing_plan(F, FusionRule(Nijk))])
end

is_gauge_fixed(data::FRData) =
    all(idx -> F_values(data)[idx] == fr_value_one(data),
        [entry.var_idx for entry in _gauge_fixing_plan(F_values(data), data)])

"""
    gauge_equivalent(F1, R1, F2, R2, Nijk)

Test gauge equivalence by comparing canonical channel-scalar representatives.
This is a deterministic F-slice based comparison, not a complete quotient
algorithm in all cases.  It is reliable when the selected F-slice is complete
or when the residual gauge action is known to be trivial on R-symbols.
"""
function gauge_equivalent(F1, R1, F2, R2, Nijk::Array{Int,3})
    c1 = canonical_gauge(F1, R1, Nijk)
    c2 = canonical_gauge(F2, R2, Nijk)
    return c1.F == c2.F && c1.R == c2.R
end

function gauge_equivalent(data1::FRData, data2::FRData)
    fusion_rule(data1).N == fusion_rule(data2).N || return false
    c1 = canonical_gauge(data1)
    c2 = canonical_gauge(data2)
    return c1.F == c2.F && c1.R == c2.R
end
