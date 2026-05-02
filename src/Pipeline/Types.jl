"""
Pipeline result types and small result helpers.

This file owns the public ClassifiedMTC record and related presentation
helpers used by the conductor-first pipeline.
"""

# ============================================================
#  Exact F/R reports
# ============================================================

@enum FRStatus begin
    FRSkipped
    FRSolved
    FRNoSolutionFound
    FRTimeoutLikeFailure
    FRReconstructionFailed
    FRVerificationFailed
end

"""
    FRRoundtripReport

Exact `(F, R)` roundtrip verification report.

Fields:
- `ok`:              `true` iff the reconstructed modular data matches the target
- `S_error`:         first nonzero exact S-entry error, or exact zero
- `T_error`:         first nonzero exact T-entry error, or exact zero
- `best_perm`:       fusion-rule automorphism giving the best match
- `S_roundtrip`:     S-matrix reconstructed from the selected R-data
- `T_roundtrip`:     T-eigenvalues reconstructed from the selected R-data
- `candidate_index`: selected F/R candidate index, when known
- `galois_exponent`: selected Galois branch exponent, when known

For compatibility with earlier `NamedTuple` reports, `S_max` and `T_max`
remain available as aliases for `S_error` and `T_error`.
"""
struct FRRoundtripReport
    ok::Bool
    S_error::Any
    T_error::Any
    best_perm::Vector{Int}
    S_roundtrip::Any
    T_roundtrip::Any
    candidate_index::Int
    galois_exponent::Int
end

function FRRoundtripReport(; ok::Bool,
                           S_error,
                           T_error,
                           best_perm,
                           S_roundtrip,
                           T_roundtrip,
                           candidate_index::Int = 0,
                           galois_exponent::Int = 1)
    return FRRoundtripReport(ok, S_error, T_error, collect(best_perm),
                             S_roundtrip, T_roundtrip,
                             candidate_index, galois_exponent)
end

function Base.getproperty(r::FRRoundtripReport, name::Symbol)
    name === :S_max && return getfield(r, :S_error)
    name === :T_max && return getfield(r, :T_error)
    return getfield(r, name)
end

function Base.propertynames(r::FRRoundtripReport, private::Bool = false)
    names = (:ok, :S_error, :T_error, :best_perm, :S_roundtrip, :T_roundtrip,
             :candidate_index, :galois_exponent, :S_max, :T_max)
    return private ? (names..., fieldnames(FRRoundtripReport)...) : names
end

function _with_fr_roundtrip_metadata(r::FRRoundtripReport;
                                     candidate_index::Int = r.candidate_index,
                                     galois_exponent::Int = r.galois_exponent)
    return FRRoundtripReport(r.ok, r.S_error, r.T_error, r.best_perm,
                             r.S_roundtrip, r.T_roundtrip,
                             candidate_index, galois_exponent)
end

# ============================================================
#  ClassifiedMTC: the final pipeline output
# ============================================================

"""
    ClassifiedMTC

The full output of the conductor-first classification pipeline for one
Galois sector.

A `ClassifiedMTC` is best read as a verified arithmetic candidate: it contains
an exact cyclotomic modular-data lift, its finite-field reconstruction
provenance, and optional exact `(F, R)` data.  When F/R reconstruction is
skipped or does not solve, the modular-data layer may still be present; inspect
`fr_status(m)` and `m.verify_report` before treating the result as a fully
realized braided fusion category.

Arithmetic / F_p and modular-data layer:
- `N`:              conductor used internally by the pipeline
- `N_input`:        user-requested conductor (for provenance)
- `rank`:           rank
- `stratum`:        the SL(2, ℤ/N) irrep decomposition `(m_λ)` that gave
                    rise to this MTC
- `Nijk`:           fusion tensor (Galois-invariant integer array)
- `scale_factor`:   the scalar multiplying `S` before rational CRT
                    reconstruction in the fallback finite-field path
- `used_primes`:    primes used for CRT reconstruction
- `fresh_primes`:   primes used for cross-validation (may be empty)
- `verify_fresh`:   `true` iff all fresh primes cross-check
- `verify_exact_lift`: `true` iff an exact fixed-stratum lift was
                    finite-field checked at all group primes; `nothing`
                    when no fixed exact lift was used

Exact modular data:
- `S_cyclotomic`:   S-matrix over `Q(ζ_N)`
- `T_cyclotomic`:   T-eigenvalues over `Q(ζ_N)`

Exact `(F, R)` layer:
- `F_values`:       exact associator coordinates over `Q(ζ_N)`, or `nothing`
                    if F/R reconstruction was skipped or inconclusive
- `R_values`:       exact braiding coordinates over `Q(ζ_N)`, or `nothing`
                    if F/R reconstruction was skipped or inconclusive
- `verify_report`:  exact roundtrip report comparing `S,T` reconstructed
                    from `F,R` against the target modular data, or `nothing`
- `fr_status`:      explicit F/R reconstruction status

Provenance:
- `galois_sector`:  integer sector index retained for provenance
"""
struct ClassifiedMTC
    N::Int
    N_input::Int
    rank::Int
    stratum::Stratum
    Nijk::Array{Int, 3}
    scale_factor::Int
    used_primes::Vector{Int}
    fresh_primes::Vector{Int}
    verify_fresh::Bool
    verify_exact_lift::Union{Bool, Nothing}
    S_cyclotomic::Any
    T_cyclotomic::Any
    F_values::Union{Vector, Nothing}
    R_values::Union{Vector, Nothing}
    verify_report::Union{FRRoundtripReport, Nothing}
    galois_sector::Int
    fr_status::FRStatus
end

function ClassifiedMTC(N::Int, N_input::Int, rank::Int, stratum::Stratum,
                       Nijk::Array{Int, 3}, scale_factor::Int,
                       used_primes::Vector{Int}, fresh_primes::Vector{Int},
                       verify_fresh::Bool, verify_exact_lift::Union{Bool, Nothing},
                       S_cyclotomic, T_cyclotomic,
                       F_values::Union{Vector, Nothing},
                       R_values::Union{Vector, Nothing},
                       verify_report::Union{FRRoundtripReport, Nothing},
                       galois_sector::Int)
    status = F_values === nothing || R_values === nothing ? FRSkipped :
             verify_report !== nothing && !verify_report.ok ? FRVerificationFailed :
             FRSolved
    return ClassifiedMTC(N, N_input, rank, stratum, Nijk, scale_factor,
                         used_primes, fresh_primes, verify_fresh,
                         verify_exact_lift, S_cyclotomic, T_cyclotomic,
                         F_values, R_values, verify_report, galois_sector,
                         status)
end

fr_status(m::ClassifiedMTC) = m.fr_status

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
                         T = c.T_cyclotomic,
                         status = (F === nothing || R === nothing ? c.fr_status : FRSolved))
    return ClassifiedMTC(c.N, c.N_input, c.rank, c.stratum, c.Nijk,
                         c.scale_factor, c.used_primes, c.fresh_primes, c.verify_fresh,
                         c.verify_exact_lift,
                         S, T,
                         F, R, report, c.galois_sector, status)
end

function _with_fr_status(c::ClassifiedMTC, status::FRStatus)
    return ClassifiedMTC(c.N, c.N_input, c.rank, c.stratum, c.Nijk,
                         c.scale_factor, c.used_primes, c.fresh_primes, c.verify_fresh,
                         c.verify_exact_lift,
                         c.S_cyclotomic, c.T_cyclotomic,
                         c.F_values, c.R_values, c.verify_report,
                         c.galois_sector, status)
end
function Base.show(io::IO, m::ClassifiedMTC)
    FR_status = if m.F_values === nothing
        "(F,R)=none"
    else
        rep = m.verify_report
        if rep !== nothing && hasproperty(rep, :ok)
            rep.ok ? "(F,R) roundtrip=✓" :
            (iszero(rep.S_error) ? "(F,R) roundtrip=S✓/T✗" :
             "(F,R) roundtrip=✗")
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
