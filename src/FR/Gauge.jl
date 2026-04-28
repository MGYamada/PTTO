"""
Conservative gauge helpers for exact `(F, R)` data.

These public functions currently expose the identity gauge orbit machinery
used by the pipeline.  They deliberately do not claim to solve the full
monoidal gauge-equivalence problem yet; instead they provide a stable API
surface and exact equality semantics that future canonicalization can refine.
"""

function _identity_gauge(Nijk::Array{Int,3})
    r = size(Nijk, 1)
    return Dict{Tuple{Int, Int, Int}, Int}(
        (a, b, c) => 1 for a in 1:r, b in 1:r, c in 1:r if Nijk[a, b, c] != 0
    )
end

"""
    canonical_gauge(F, R, Nijk)

Return a conservative canonical representative for exact `(F, R)` data.
At present this is the identity-gauge representative.
"""
function canonical_gauge(F, R, Nijk::Array{Int,3})
    return (F = copy(F), R = copy(R), gauge = _identity_gauge(Nijk))
end

"""
    gauge_transform(F, R, gauge)

Apply a gauge transform to `(F, R)`.  The current public implementation accepts
identity gauges and returns copied data.
"""
function gauge_transform(F, R, gauge)
    if gauge isa AbstractDict
        all(v -> v == 1, values(gauge)) ||
            error("nontrivial gauge transforms are not implemented yet")
    elseif gauge !== nothing
        error("unsupported gauge object: $(typeof(gauge))")
    end
    return (F = copy(F), R = copy(R))
end

"""
    gauge_equivalent(F1, R1, F2, R2, Nijk)

Conservatively test gauge equivalence by comparing identity-gauge canonical
representatives.  This is exact equality, not a full gauge-orbit search.
"""
function gauge_equivalent(F1, R1, F2, R2, Nijk::Array{Int,3})
    c1 = canonical_gauge(F1, R1, Nijk)
    c2 = canonical_gauge(F2, R2, Nijk)
    return c1.F == c2.F && c1.R == c2.R
end
