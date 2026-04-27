"""
Stratum enumeration.

Given an atomic catalog of SL(2, ℤ/N) irreps and a target rank r, enumerate
all possible "strata" — i.e., assignments of multiplicity m_λ to each atomic
irrep such that Σ m_λ · dim(λ) = r.

In the NRWW scheme-theoretic framework, each stratum corresponds to a locally
closed subscheme of the representation variety Hom(SL(2, ℤ/N), GL_r), indexed
by the isomorphism type of the representation as a direct sum of irreducibles.

This is a pure combinatorial step (dim-based partition); T-spectrum matching
and block-U enumeration come later.

Key type: `Stratum`
- `multiplicities`: Dict{Int, Int} mapping catalog index → multiplicity m_λ
- `total_dim`: sum of m_λ · dim(λ), equals r by construction
"""

"""
    Stratum

A stratum of the representation variety = an assignment of multiplicities
to atomic irreps.

Fields:
- multiplicities: Dict mapping catalog index (1-based) to multiplicity m_λ ≥ 1
- total_dim: sum of m_λ · dim(λ) (redundant but cached)
"""
struct Stratum
    multiplicities::Dict{Int, Int}
    total_dim::Int
end

function Base.show(io::IO, s::Stratum)
    parts = sort(collect(s.multiplicities); by = first)
    desc = join(["[$idx]×$m" for (idx, m) in parts], " + ")
    print(io, "Stratum($desc, total_dim=$(s.total_dim))")
end

"""
    describe_stratum(s::Stratum, catalog) -> String

Human-readable description of a stratum using catalog labels.
E.g., "3d_8 + 2×1d_3 + 1d_1" for a stratum with a 3-dimensional irrep at level 8
and two copies of the non-trivial 1-dim irrep at level 3, plus the trivial irrep.
"""
function describe_stratum(s::Stratum, catalog)
    parts = String[]
    for (idx, m) in sort(collect(s.multiplicities); by = first)
        label = catalog[idx].label
        parity_tag = catalog[idx].parity == +1 ? "+" : "-"
        entry = m == 1 ? "$(label)$(parity_tag)" : "$m×$(label)$(parity_tag)"
        push!(parts, entry)
    end
    return join(parts, " ⊕ ")
end

"""
    find_unit_indices(catalog) -> Vector{Int}

Return indices of the trivial irrep (dim 1, level 1, parity +1).
The unit object of the MTC must transform under this irrep.
There should be exactly one such index in a well-formed catalog.
"""
function find_unit_indices(catalog)
    return [i for (i, a) in enumerate(catalog)
            if a.dim == 1 && a.level == 1 && a.parity == +1]
end

"""
    enumerate_strata(catalog, r::Int; require_unit_summand::Bool = false,
                     max_multiplicity::Int = typemax(Int)) -> Vector{Stratum}

Enumerate all strata {m_λ} with Σ m_λ · dim(λ) = r using atomic irreps from
the catalog.

NOTE on unit objects: In general, the unit object of an MTC sits as a
specific basis vector WITHIN an irreducible block (whichever irrep contains
T-eigenvalue 1), not as a separate 1d_1 summand. For example, SU(2)_4 at N=24
has rank 5 decomposition (ρ_8^(3) ⊠ 1_3) ⊕ (1_8 ⊠ ρ_3^(2)) with NO explicit
1d_1 summand, yet two objects have spin 0 (one in each block).

Therefore `require_unit_summand=false` is the correct default. The unit
constraint is enforced later, at the T-spectrum / S-matrix assembly step.

If `require_unit_summand=true`, strata without at least one copy of the
trivial irrep (1d_1, parity=+1) are excluded. Use only for abelian/pointed
categories or debugging.

This is a pure combinatorial enumeration by dimension. T-spectrum matching
and block-U enumeration are subsequent filters.
"""
function enumerate_strata(catalog, r::Int;
                          require_unit_summand::Bool = false,
                          max_multiplicity::Int = typemax(Int))
    results = Stratum[]
    n = length(catalog)
    dims = [a.dim for a in catalog]

    unit_indices = find_unit_indices(catalog)
    if require_unit_summand && isempty(unit_indices)
        error("No trivial irrep (1d_1, parity=+1) in catalog; cannot enforce unit summand")
    end

    function recurse!(results::Vector{Stratum},
                      current::Dict{Int, Int},
                      idx::Int,
                      remaining::Int)
        if remaining == 0
            if require_unit_summand
                has_unit = any(haskey(current, u) && current[u] >= 1
                              for u in unit_indices)
                has_unit || return
            end
            push!(results, Stratum(copy(current), r))
            return
        end
        if idx > n
            return
        end

        d = dims[idx]
        recurse!(results, current, idx + 1, remaining)

        max_m = min(remaining ÷ d, max_multiplicity)
        for m in 1:max_m
            current[idx] = m
            recurse!(results, current, idx + 1, remaining - m * d)
            delete!(current, idx)
        end
    end

    recurse!(results, Dict{Int, Int}(), 1, r)
    return results
end

"""
    count_strata(catalog, r::Int; kwargs...) -> Int

Count strata without materializing the full list (useful for estimating
enumeration size before committing resources).
"""
function count_strata(catalog, r::Int; kwargs...)
    # Currently a wrapper; could be optimized with DP
    return length(enumerate_strata(catalog, r; kwargs...))
end
