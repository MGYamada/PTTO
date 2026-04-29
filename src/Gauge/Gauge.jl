"""
Gauge helpers for exact multiplicity-free `(F, R)` data.

For a fusion channel basis vector `V_ab^c`, a gauge scalar `u[a,b,c]`
acts by

    F'^{ijk}_{o;a,b} = u[i,j,a] u[a,k,o] / (u[j,k,b] u[i,b,o]) F^{ijk}_{o;a,b}
    R'^{ij}_k       = u[j,i,k] / u[i,j,k] R^{ij}_k

with unit channels normalized to `1`.  The current implementation is
intentionally exact and multiplicity-free; higher multiplicity needs full
change-of-basis matrices rather than channel scalars.
"""

struct GaugeTransform{T}
    scalars::Dict{Tuple{Int,Int,Int}, T}
    fixed_indices::Vector{Int}
    complete::Bool
end

function _identity_gauge(Nijk::Array{Int,3})
    r = size(Nijk, 1)
    return Dict{Tuple{Int, Int, Int}, Int}(
        (a, b, c) => 1 for a in 1:r, b in 1:r, c in 1:r if Nijk[a, b, c] != 0
    )
end

function _one_like_from_FR(F, R)
    !isempty(F) && return one(F[1])
    !isempty(R) && return one(R[1])
    return 1
end

function _one_gauge_scalars(Nijk::Array{Int,3}, one_value)
    r = size(Nijk, 1)
    scalars = Dict{Tuple{Int,Int,Int}, typeof(one_value)}()
    for a in 1:r, b in 1:r, c in 1:r
        Nijk[a, b, c] == 0 && continue
        scalars[(a, b, c)] = one_value
    end
    return scalars
end

function _require_multiplicity_free_gauge(Nijk::Array{Int,3})
    _multiplicity_free_fusion(Nijk) ||
        error("channel-scalar gauge transforms require multiplicity-free fusion rules")
    return true
end

function _gauge_scalar(gauge, ch::Tuple{Int,Int,Int}, default_one)
    if gauge isa GaugeTransform
        return get(gauge.scalars, ch, default_one)
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

function _ordered_gauge_channels(Nijk::Array{Int,3}, r::Int)
    indexed = collect(_gauge_channel_basis(Nijk, r))
    sort!(indexed; by = kv -> kv.second)
    return [kv.first for kv in indexed]
end

function _r_var_count(Nijk::Array{Int,3})
    r = size(Nijk, 1)
    return sum(Nijk[i, j, k]^2 for i in 1:r, j in 1:r, k in 1:r)
end

function _r_forward_transform_factor(Nijk::Array{Int,3},
                                     i::Int, j::Int, k::Int,
                                     gauge,
                                     one_value)
    numerator = _gauge_scalar(gauge, (j, i, k), one_value)
    denominator = _gauge_scalar(gauge, (i, j, k), one_value)
    return numerator / denominator
end

function _r_inverse_transform_factor(Nijk::Array{Int,3},
                                     i::Int, j::Int, k::Int,
                                     gauge,
                                     one_value)
    numerator = _gauge_scalar(gauge, (i, j, k), one_value)
    denominator = _gauge_scalar(gauge, (j, i, k), one_value)
    return numerator / denominator
end

function _transform_F_values(F, Nijk::Array{Int,3}, gauge, one_value)
    r = size(Nijk, 1)
    metadata = _pentagon_variable_metadata(Nijk, r, length(F))
    channel_index = _gauge_channel_basis(Nijk, r)
    channels = _ordered_gauge_channels(Nijk, r)
    out = copy(F)
    for idx in eachindex(F)
        weight = _gauge_weight(metadata[idx], channel_index)
        out[idx] = F[idx] * _factor_from_weight(weight, channels,
                                                gauge isa GaugeTransform ? gauge.scalars : gauge,
                                                one_value)
    end
    return out
end

function _transform_R_values(R, Nijk::Array{Int,3}, gauge, one_value)
    r = size(Nijk, 1)
    n_forward = _r_var_count(Nijk)
    (length(R) == n_forward || length(R) == 2n_forward) ||
        error("R has length $(length(R)); expected $n_forward or $(2n_forward)")
    positions, _ = _braiding_block_positions(Nijk)
    out = copy(R)
    for i in 1:r, j in 1:r, k in 1:r
        Nijk[i, j, k] == 0 && continue
        fwd = _r_forward_transform_factor(Nijk, i, j, k, gauge, one_value)
        for pos in positions[(i, j, k)]
            out[pos] = R[pos] * fwd
        end
        if length(R) == 2n_forward
            inv_factor = _r_inverse_transform_factor(Nijk, i, j, k, gauge, one_value)
            for pos in positions[(i, j, k)]
                out[n_forward + pos] = R[n_forward + pos] * inv_factor
            end
        end
    end
    return out
end

function _gauge_fixing_plan(F, Nijk::Array{Int,3})
    r = size(Nijk, 1)
    _require_multiplicity_free_gauge(Nijk)
    metadata = _pentagon_variable_metadata(Nijk, r, length(F))
    channel_index = _gauge_channel_basis(Nijk, r)
    channels = _ordered_gauge_channels(Nijk, r)
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
gauge_fixing_plan(F, Nijk::Array{Int,3}) = _gauge_fixing_plan(F, Nijk)

function _canonical_gauge_transform(F, R, Nijk::Array{Int,3})
    _require_multiplicity_free_gauge(Nijk)
    one_value = _one_like_from_FR(F, R)
    scalars = _one_gauge_scalars(Nijk, one_value)
    channels = _ordered_gauge_channels(Nijk, size(Nijk, 1))
    plan = _gauge_fixing_plan(F, Nijk)

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
    complete = length(fixed) == length(_select_pentagon_gauge_fixed_indices(Nijk, size(Nijk, 1), length(F)))
    return GaugeTransform(scalars, fixed, complete)
end

"""
    canonical_gauge(F, R, Nijk)

Return a canonical channel-scalar gauge representative for multiplicity-free
exact `(F, R)` data.  Independent nonzero F-symbol coordinates are fixed to
`1` whenever this can be done with exact multiplication and inversion.
"""
function canonical_gauge(F, R, Nijk::Array{Int,3})
    gauge = _canonical_gauge_transform(F, R, Nijk)
    transformed = gauge_transform(F, R, gauge; Nijk = Nijk)
    return (F = transformed.F, R = transformed.R, gauge = gauge)
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
    _require_multiplicity_free_gauge(Nijk)
    one_value = _one_like_from_FR(F, R)
    return (F = _transform_F_values(F, Nijk, gauge, one_value),
            R = _transform_R_values(R, Nijk, gauge, one_value))
end

"""
    is_gauge_fixed(F, Nijk)

Return `true` when all coordinates selected by `gauge_fixing_plan(F, Nijk)`
are exactly `1`.
"""
function is_gauge_fixed(F, Nijk::Array{Int,3})
    one_value = _one_like_from_FR(F, [])
    return all(idx -> F[idx] == one_value, [entry.var_idx for entry in _gauge_fixing_plan(F, Nijk)])
end

"""
    gauge_equivalent(F1, R1, F2, R2, Nijk)

Test gauge equivalence by comparing canonical channel-scalar representatives.
"""
function gauge_equivalent(F1, R1, F2, R2, Nijk::Array{Int,3})
    c1 = canonical_gauge(F1, R1, Nijk)
    c2 = canonical_gauge(F2, R2, Nijk)
    return c1.F == c2.F && c1.R == c2.R
end
