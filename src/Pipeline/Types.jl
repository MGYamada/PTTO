"""
Pipeline result types and small result helpers.

This file owns the public ClassifiedMTC record and related presentation
helpers used by the conductor-first pipeline.
"""

# ============================================================
#  ClassifiedMTC: the final pipeline output
# ============================================================

"""
    ClassifiedMTC

The full output of `classify_mtcs_at_conductor` for a single MTC.

Arithmetic / F_p layer (Phase 0–3 output):
- `N`:              conductor used internally by the pipeline
- `N_input`:        user-requested conductor (for provenance)
- `rank`:           rank
- `stratum`:        the SL(2, ℤ/N) irrep decomposition `(m_λ)` that gave
                    rise to this MTC
- `Nijk`:           fusion tensor (Galois-invariant integer array)
- `S_Zsqrtd`:       S-matrix as `(a, b)` = `a + b·√d` in ℤ[√d]
                    (after Phase 3 CRT; entries before dividing by
                    `scale · √d`)
- `quadratic_d`:    the inferred `d` such that `S_Zsqrtd` lives in ℤ[√d]
- `scale_factor`:   the scalar multiplying `S` before reconstruction
                    (so `S_cyclotomic = S_Zsqrtd / (scale_factor · √d)`)
- `used_primes`:    primes used for CRT reconstruction
- `fresh_primes`:   primes used for cross-validation (may be empty)
- `verify_fresh`:   `true` iff all fresh primes cross-check
- `verify_exact_lift`: `true` iff an exact fixed-stratum lift was
                    finite-field checked at all group primes; `nothing`
                    when no fixed exact lift was used

Exact modular data (Phase 4 input):
- `S_cyclotomic`:   S-matrix over `Q(ζ_N)`
- `T_cyclotomic`:   T-eigenvalues over `Q(ζ_N)`

Exact `(F, R)` layer:
- `F_values`:       reserved for exact pentagon data; currently `nothing`.
- `R_values`:       reserved for exact braiding data; currently `nothing`.
- `verify_report`:  reserved for exact verification data; currently `nothing`.

Provenance:
- `galois_sector`:  integer index of the Galois orbit element
                    (1, 2, ... for groups returned by
                    `group_mtcs_galois_aware`)
"""
struct ClassifiedMTC
    N::Int
    N_input::Int
    rank::Int
    stratum::Stratum
    Nijk::Array{Int, 3}
    S_Zsqrtd::Matrix{Tuple{Int, Int}}
    quadratic_d::Int
    scale_factor::Int
    used_primes::Vector{Int}
    fresh_primes::Vector{Int}
    verify_fresh::Bool
    verify_exact_lift::Union{Bool, Nothing}
    S_cyclotomic::Any
    T_cyclotomic::Any
    F_values::Union{Vector, Nothing}
    R_values::Union{Vector, Nothing}
    verify_report::Any
    galois_sector::Int
end

function _permute_fusion_tensor(Nijk::Array{Int, 3}, perm::Vector{Int})
    return Nijk[perm, perm, perm]
end

function _all_permutations(v::Vector{Int})
    if length(v) <= 1
        return [copy(v)]
    end
    out = Vector{Vector{Int}}()
    for i in eachindex(v)
        head = v[i]
        tail = Vector{Int}(undef, length(v) - 1)
        ti = 1
        for j in eachindex(v)
            if j != i
                tail[ti] = v[j]
                ti += 1
            end
        end
        for rest in _all_permutations(tail)
            push!(out, vcat(head, rest))
        end
    end
    return out
end

"""
    fusion_rule_key(Nijk) -> String

Backward-compatible alias for `canonical_rule(Nijk)`.
"""
function fusion_rule_key(Nijk::AbstractArray{<:Integer, 3})
    return canonical_rule(Nijk)
end

function _classify_modular_data_by_fusion_rule(classified::Vector{ClassifiedMTC})
    grouped = Dict{String, Vector{Int}}()
    for (i, cmtc) in enumerate(classified)
        key = fusion_rule_key(cmtc.Nijk)
        push!(get!(grouped, key, Int[]), i)
    end
    return grouped
end

function _with_fr_result(c::ClassifiedMTC,
                         F::Union{Vector, Nothing},
                         R::Union{Vector, Nothing},
                         report;
                         S = c.S_cyclotomic,
                         T = c.T_cyclotomic)
    return ClassifiedMTC(c.N, c.N_input, c.rank, c.stratum, c.Nijk,
                         c.S_Zsqrtd, c.quadratic_d, c.scale_factor,
                         c.used_primes, c.fresh_primes, c.verify_fresh,
                         c.verify_exact_lift,
                         S, T,
                         F, R, report, c.galois_sector)
end
function Base.show(io::IO, m::ClassifiedMTC)
    FR_status = if m.F_values === nothing
        "(F,R)=none"
    else
        rep = m.verify_report
        if rep !== nothing && hasproperty(rep, :ok)
            "(F,R) roundtrip=$(rep.ok ? "✓" : "✗")"
        else
            "(F,R)=attached"
        end
    end
    fresh_str = isempty(m.fresh_primes) ? "no fresh" :
        (m.verify_fresh ? "fresh✓" : "fresh✗")
    exact_str = m.verify_exact_lift === nothing ? "" :
        (m.verify_exact_lift ? ", exact✓" : ", exact✗")
    n_desc = m.N == m.N_input ? string(m.N) : "$(m.N) [input=$(m.N_input)]"
    print(io, "ClassifiedMTC(N=$(n_desc), rank=$(m.rank), ",
          "sector=$(m.galois_sector), $(length(m.used_primes)) primes, ",
          "$fresh_str$exact_str, $FR_status)")
end
