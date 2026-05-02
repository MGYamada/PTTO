"""
Stratum enumeration for conductor-level modular-data search.

Fix a conductor `N`, write `G_N = SL(2, ℤ/N)`, and let the atomic catalog be a
chosen list of irreducible `G_N`-representations `ρ_λ` realized over
`Q(ζ_N)`.  For a target rank `r`, a stratum is a multiplicity function

```math
m : λ \\mapsto m_λ \\in ℤ_{\\ge 0},
\\qquad
\\sum_λ m_λ \\dim ρ_λ = r.
```

It indexes the semisimple representation isomorphism type

```math
ρ_m \\cong \\bigoplus_λ ρ_λ^{\\oplus m_λ}.
```

Equivalently, in the representation variety of `G_N` in dimension `r`, the
stratum is the locally closed locus with this fixed Jordan-Hölder
decomposition.  In ACMG this step is purely combinatorial.  The unit object,
T-spectrum constraints, Verlinde integrality, and Block-U equations are
imposed later.

Key type: `Stratum`
- `multiplicities`: `Dict{Int, Int}` mapping catalog index to `m_λ`
- `total_dim`: cached value of `Σ m_λ dim(ρ_λ)`, equal to the searched rank
"""

"""
    Stratum

A semisimple representation stratum for `SL(2, ℤ/N)` at fixed rank.

If the atomic catalog contains irreducibles `ρ_λ`, a `Stratum` stores the
finite-support multiplicity vector `(m_λ)`.  The corresponding representation
type is `⊕_λ ρ_λ^{⊕m_λ}`.  This is not yet a modular tensor category: it is the
ambient representation-theoretic layer in which candidate modular data are
searched.

Fields:
- `multiplicities`: dictionary mapping catalog index (1-based) to positive
  multiplicity `m_λ`
- `total_dim`: sum of `m_λ · dim(ρ_λ)` (redundant but cached)
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

Human-readable description of a stratum using catalog labels.  The result is a
compact direct-sum expression such as `"3d_8+ ⊕ 2×1d_3+"`, where each term
records multiplicity, atomic label, and parity sign.
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
This is a useful diagnostic, but a valid MTC need not contain the unit object
as a separate trivial direct summand of the `SL(2, ℤ/N)` representation.  The
unit is ultimately a basis vector in the modular-data basis and may lie inside
a higher-dimensional irreducible summand.
"""
function find_unit_indices(catalog)
    return [i for (i, a) in enumerate(catalog)
            if a.dim == 1 && a.level == 1 && a.parity == +1]
end

"""
    enumerate_strata(catalog, r::Int; require_unit_summand::Bool = false,
                     max_multiplicity::Int = typemax(Int)) -> Vector{Stratum}

Enumerate all strata `{m_λ}` with `Σ m_λ · dim(ρ_λ) = r` using atomic irreps
from the catalog.

Mathematical note on the unit object: in general, the tensor unit is a
distinguished simple object, hence a distinguished basis vector after modular
data are assembled.  It need not appear as a separate trivial summand of the
`SL(2, ℤ/N)` representation.  It may lie inside a nontrivial irreducible block
that contains a `T`-eigenvalue equal to `1`.

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
