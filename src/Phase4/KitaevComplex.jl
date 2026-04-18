"""
    KitaevComplex

Linearised gauge-orbit analysis at a fixed F-symbol solution, following
Kitaev 2006 Appendix E.6 (Crane–Yetter / Davydov tangent cohomology,
"Yetter cohomology"). Provides the canonical Hodge slice for gauge fixing.

# Approach

Operates in F-symbol coordinates. Given a fusion rule and a base
F-symbol function, linearises pentagon, gauge action, and unit axiom at
the base point and returns:

  * `Δ_gauge`    : gauge tangent  →  F tangent
  * `Δ_pent`     : F tangent      →  pentagon equation tangent
  * `Unit`       : unit-axiom constraints (δF = 0 on vacuum-containing entries)
  * `δ¹`         : object-phase → vertex-gauge map  (Kitaev Eq. 244)
  * Effective gauge orbit (after restricting to unit-preserving subspace)
  * Orthonormal slice basis of the F-tangent gauge fixing
  * Numerical check of Ocneanu rigidity (H³ = 0)

# Use case

For HC pentagon solving at large rank, the bottleneck is the gauge orbit:
the polynomial system has continuous gauge directions that inflate the
mixed volume. The linear constraints from `slice_basis(ga)` can be added
to the pentagon system to kill the gauge directions while preserving
algebraicity (the constraints are linear in F).

# References

  Kitaev, A. "Anyons in an exactly solved model and beyond."
    Annals Phys. 321 (2006), App. E.6.

"""
module KitaevComplex

using LinearAlgebra
using ACMG: FusionRule

export FCoordSpace, build_F_coord_space, F_var_count, F_value
export build_Delta_gauge, build_Delta_pent
export build_unit_constraints, build_delta1, enumerate_vertices
export GaugeAnalysis, analyze_gauge
export gauge_orbit_dim, slice_dim, H3_dimension
export effective_gauge, slice_basis
export verify_ocneanu_rigidity

# ============================================================================
# F-symbol entry indexing
# ============================================================================
#
# An F-symbol entry is identified by a 6-tuple (a,b,c,d,e,f):
#   - a, b, c, d are external simple-object labels (1-indexed)
#   - e is the (a⊗b → e) channel
#   - f is the (b⊗c → f) channel
# Consistency: e⊗c → d  AND  a⊗f → d
#
# F^{abc}_{d; e, f}  is the value at this entry.
# ============================================================================

const FKey = NTuple{6, Int}

"""
    is_F_allowed(fr, a, b, c, d, e, f) -> Bool

Whether the F-entry F^{abc}_{d; e, f} is structurally allowed by fusion.
Requires all four channels to have nonzero multiplicity.
"""
function is_F_allowed(fr::FusionRule, a::Int, b::Int, c::Int,
                                       d::Int, e::Int, f::Int)
    return fr.N[a, b, e] ≥ 1 &&
           fr.N[e, c, d] ≥ 1 &&
           fr.N[b, c, f] ≥ 1 &&
           fr.N[a, f, d] ≥ 1
end

# ============================================================================
# F-coordinate space
# ============================================================================

"""
    FCoordSpace

Tangent space of F-symbols at a base solution.

Fields:
- `fr`           : fusion rule
- `vars`         : Vector{FKey} of all allowed F-entries with nonzero base value
- `var_idx`      : Dict mapping FKey → 1-based index in `vars`
- `F_value`      : function (a,b,c,d,e,f) → Float64 (base F value, 0 if disallowed)
- `unit_indices` : indices of variables forced to 1 by unit axiom
                   (any of (a,b,c,d) equals the unit object 1)
"""
struct FCoordSpace
    fr::FusionRule
    vars::Vector{FKey}
    var_idx::Dict{FKey, Int}
    F_value::Function
    unit_indices::Vector{Int}
end

"""
    build_F_coord_space(fr, F_value::Function) -> FCoordSpace

Enumerate all F-entries allowed by fusion with `F_value(a,b,c,d,e,f) ≠ 0`.
"""
function build_F_coord_space(fr::FusionRule, F_value::Function)
    r = fr.rank
    vars = FKey[]
    for a in 1:r, b in 1:r, c in 1:r, d in 1:r, e in 1:r, f in 1:r
        if is_F_allowed(fr, a, b, c, d, e, f)
            v = F_value(a, b, c, d, e, f)
            if v != 0.0
                push!(vars, (a, b, c, d, e, f))
            end
        end
    end
    var_idx = Dict{FKey, Int}(v => k for (k, v) in enumerate(vars))
    unit_idx = Int[]
    for (k, key) in enumerate(vars)
        a, b, c, d, _, _ = key
        if 1 in (a, b, c, d)
            push!(unit_idx, k)
        end
    end
    return FCoordSpace(fr, vars, var_idx, F_value, unit_idx)
end

F_var_count(fcs::FCoordSpace) = length(fcs.vars)
F_value(fcs::FCoordSpace, key::FKey) = fcs.F_value(key...)

# ============================================================================
# Vertex gauge enumeration
# ============================================================================

const VKey = NTuple{3, Int}

"""
    enumerate_vertices(fr) -> (Vector{VKey}, Dict{VKey,Int})

Non-trivial fusion vertices (a, b, c) with c ∈ a⊗b, a ≠ 1, b ≠ 1.
Vacuum vertices are normalized to u = 1 and excluded.
"""
function enumerate_vertices(fr::FusionRule)
    r = fr.rank
    verts = VKey[]
    for a in 2:r, b in 2:r, c in 1:r
        fr.N[a, b, c] ≥ 1 && push!(verts, (a, b, c))
    end
    vert_idx = Dict{VKey, Int}(v => k for (k, v) in enumerate(verts))
    return verts, vert_idx
end

# ============================================================================
# Δ_gauge : gauge tangent → F tangent
# ============================================================================
#
# Linearised gauge action (Kitaev Eq. 239 with Φ from Eq. 242, restricted to
# vertex gauge u^{ab}_c = 1 + δu):
#
#   F^{abc}_{d;e,f}  ↦  (u^{ab}_e u^{ec}_d / (u^{bc}_f u^{af}_d)) · F
#   δF / F  =  δu^{ab}_e + δu^{ec}_d - δu^{bc}_f - δu^{af}_d
#
# So at row (a,b,c,d,e,f) and column = vertex index, the matrix entry is
# ±F-value if the vertex matches one of the four, else 0. Vacuum vertices
# (a=1 or b=1) are u=1 fixed → no contribution.
# ============================================================================

"""
    build_Delta_gauge(fcs) -> (Matrix{Float64}, Vector{VKey})

Linearised gauge action matrix (size nF × nG) and the vertex list.
"""
function build_Delta_gauge(fcs::FCoordSpace)
    fr = fcs.fr
    verts, vert_idx = enumerate_vertices(fr)
    nF = F_var_count(fcs)
    nG = length(verts)
    M = zeros(Float64, nF, nG)

    for (row, key) in enumerate(fcs.vars)
        a, b, c, d, e, f = key
        Fv = F_value(fcs, key)
        Fv == 0.0 && continue
        for (sign, vert) in ((+1, (a, b, e)), (+1, (e, c, d)),
                             (-1, (b, c, f)), (-1, (a, f, d)))
            haskey(vert_idx, vert) && (M[row, vert_idx[vert]] += sign * Fv)
        end
    end
    return M, verts
end

# ============================================================================
# Δ_pent : F tangent → pentagon tangent
# ============================================================================
#
# Pentagon (standard convention):
#   [F^{fcd}_e]_{g,l} · [F^{abl}_e]_{f,k}
#       = Σ_h [F^{abc}_g]_{f,h} · [F^{ahd}_e]_{g,k} · [F^{bcd}_k]_{h,l}
#
# Each (a,b,c,d,e,f,g,k,l) gives one polynomial equation in F-vars.
# Linearise at the base F-solution.
# ============================================================================

"""
    build_Delta_pent(fcs) -> Matrix{Float64}

Linearised pentagon Jacobian, of size (npent_rows × nF). Only nonzero
rows are kept.
"""
function build_Delta_pent(fcs::FCoordSpace)
    fr = fcs.fr
    r = fr.rank
    nF = F_var_count(fcs)
    rows = Vector{Vector{Float64}}()

    for a in 1:r, b in 1:r, c in 1:r, d in 1:r, e in 1:r,
        f in 1:r, g in 1:r, k in 1:r, l in 1:r

        jac = zeros(Float64, nF)
        has_any = false

        # LHS (sign +1): F^{fcd}_{e, g, l} · F^{abl}_{e, f, k}
        lhs_factors = ((f, c, d, e, g, l), (a, b, l, e, f, k))
        lhs_vals = [F_value(fcs, t) for t in lhs_factors]
        if all(!iszero, lhs_vals)
            has_any = true
            for (i, t) in enumerate(lhs_factors)
                if haskey(fcs.var_idx, t)
                    partial = prod(lhs_vals[j] for j in eachindex(lhs_factors) if j != i)
                    jac[fcs.var_idx[t]] += +partial
                end
            end
        end

        # RHS (sign -1, summed over h): F^{abc}_{g,f,h} · F^{ahd}_{e,g,k} · F^{bcd}_{k,h,l}
        for h in 1:r
            rhs_factors = ((a, b, c, g, f, h), (a, h, d, e, g, k), (b, c, d, k, h, l))
            rhs_vals = [F_value(fcs, t) for t in rhs_factors]
            if all(!iszero, rhs_vals)
                has_any = true
                for (i, t) in enumerate(rhs_factors)
                    if haskey(fcs.var_idx, t)
                        partial = prod(rhs_vals[j] for j in eachindex(rhs_factors) if j != i)
                        jac[fcs.var_idx[t]] += -partial
                    end
                end
            end
        end

        if has_any && any(x -> abs(x) > 1e-14, jac)
            push!(rows, jac)
        end
    end

    isempty(rows) && return zeros(Float64, 0, nF)
    M = zeros(Float64, length(rows), nF)
    for (i, r_) in enumerate(rows)
        M[i, :] = r_
    end
    return M
end

# ============================================================================
# Unit-axiom constraints
# ============================================================================

"""
    build_unit_constraints(fcs) -> Matrix{Float64}

Linear constraints δF = 0 for every F-variable forced by unit axioms
(any of a, b, c, d equals 1).
"""
function build_unit_constraints(fcs::FCoordSpace)
    nF = F_var_count(fcs)
    n_unit = length(fcs.unit_indices)
    M = zeros(Float64, n_unit, nF)
    for (row, idx) in enumerate(fcs.unit_indices)
        M[row, idx] = 1.0
    end
    return M
end

# ============================================================================
# δ¹ : object phases → vertex gauge   (Kitaev Eq. 244)
# ============================================================================

"""
    build_delta1(fr, verts) -> Matrix{Float64}

Map δ¹ : R^r → R^nG sending (X_1, …, X_r) to the vertex-gauge deformation
δu^{ab}_c = X_a + X_b - X_c at each vertex.
"""
function build_delta1(fr::FusionRule, verts::Vector{VKey})
    r = fr.rank
    nG = length(verts)
    M = zeros(Float64, nG, r)
    for (row, (a, b, c)) in enumerate(verts)
        M[row, a] += 1.0
        M[row, b] += 1.0
        M[row, c] -= 1.0
    end
    return M
end

# ============================================================================
# Gauge analysis: full pipeline
# ============================================================================

"""
    GaugeAnalysis

Output of `analyze_gauge`. See field docs in source.
"""
struct GaugeAnalysis
    fcs::FCoordSpace
    verts::Vector{VKey}
    Delta_gauge::Matrix{Float64}            # nF × nG
    Delta_pent::Matrix{Float64}             # npent × nF
    Unit::Matrix{Float64}                   # nUnit × nF
    delta1::Matrix{Float64}                 # nG × r
    unit_pres_gauge_basis::Matrix{Float64}  # nG × kU; columns = ker(U·Δ_g)
    Delta_gauge_eff::Matrix{Float64}        # nF × kU
    gauge_orbit_dim::Int
    slice_basis_::Matrix{Float64}           # nF × s; orthonormal slice
    H3_dim::Int
end

"""
    analyze_gauge(fcs; tol = 1e-10) -> GaugeAnalysis

Full Kitaev-Hodge analysis: pentagon, unit, gauge linearizations + slice.
"""
function analyze_gauge(fcs::FCoordSpace; tol::Real = 1e-10)
    Dg, verts = build_Delta_gauge(fcs)
    Dp = build_Delta_pent(fcs)
    U = build_unit_constraints(fcs)
    d1 = build_delta1(fcs.fr, verts)
    nF = F_var_count(fcs)
    nG = length(verts)

    # Unit-preserving gauge subspace = ker(U · Dg)
    unit_pres_basis = if nG == 0 || size(U, 1) == 0
        Matrix{Float64}(I, nG, nG)
    else
        UD = U * Dg
        rank_UD = _matrix_rank(UD; tol = tol)
        if rank_UD == nG
            zeros(Float64, nG, 0)
        elseif rank_UD == 0
            Matrix{Float64}(I, nG, nG)
        else
            F = svd(UD; full = true)
            F.V[:, rank_UD+1:end]
        end
    end

    Dg_eff = Dg * unit_pres_basis
    rank_Dg_eff = _matrix_rank(Dg_eff; tol = tol)

    # Permissible deformation = ker [Dp; U]
    Combined = isempty(U) ? Dp : vcat(Dp, U)
    rank_Combined = _matrix_rank(Combined; tol = tol)
    ker_Combined_dim = nF - rank_Combined
    H3_dim = ker_Combined_dim - rank_Dg_eff

    # Slice basis: ker(U) ∩ (im Δ_gauge_eff)^⊥
    unit_compat = if size(U, 1) == 0
        Matrix{Float64}(I, nF, nF)
    else
        F_U = svd(U; full = true)
        rank_U = _matrix_rank(U; tol = tol)
        F_U.V[:, rank_U+1:end]
    end

    Dg_in_unit = unit_compat' * Dg_eff
    slice_in_compat = if iszero(Dg_in_unit) || isempty(Dg_in_unit)
        Matrix{Float64}(I, size(unit_compat, 2), size(unit_compat, 2))
    else
        F_g = svd(Dg_in_unit; full = true)
        rank_g = _matrix_rank(Dg_in_unit; tol = tol)
        F_g.U[:, rank_g+1:end]
    end
    slice_basis_mat = unit_compat * slice_in_compat

    return GaugeAnalysis(fcs, verts, Dg, Dp, U, d1,
                         unit_pres_basis, Dg_eff,
                         rank_Dg_eff, slice_basis_mat, H3_dim)
end

# ----- Accessors -----
gauge_orbit_dim(ga::GaugeAnalysis) = ga.gauge_orbit_dim
slice_dim(ga::GaugeAnalysis) = size(ga.slice_basis_, 2)
H3_dimension(ga::GaugeAnalysis) = ga.H3_dim
effective_gauge(ga::GaugeAnalysis) = ga.Delta_gauge_eff
slice_basis(ga::GaugeAnalysis) = ga.slice_basis_

"""
    verify_ocneanu_rigidity(ga; tol = 1e-9) -> Bool

Verify ker(pentagon ∩ unit) = im(effective gauge), i.e., H³ = 0
within tolerance.
"""
function verify_ocneanu_rigidity(ga::GaugeAnalysis; tol::Real = 1e-9)
    ga.H3_dim != 0 && return false

    Combined = isempty(ga.Unit) ? ga.Delta_pent : vcat(ga.Delta_pent, ga.Unit)
    F_C = svd(Combined; full = true)
    rank_C = _matrix_rank(Combined; tol = tol)
    nF = size(Combined, 2)
    ker_basis = F_C.V[:, rank_C+1:end]
    Dg_eff = ga.Delta_gauge_eff

    isempty(ker_basis) && return true
    isempty(Dg_eff) && return iszero(ker_basis)

    for j in 1:size(ker_basis, 2)
        rhs = ker_basis[:, j]
        y = Dg_eff \ rhs
        resid = norm(Dg_eff * y - rhs)
        resid > tol && return false
    end
    return true
end

# ============================================================================
# Helpers
# ============================================================================

function _matrix_rank(M::AbstractMatrix; tol::Real = 1e-10)
    isempty(M) && return 0
    s = svdvals(M)
    isempty(s) && return 0
    thresh = tol * maximum(size(M)) * s[1]
    return count(>(thresh), s)
end

end # module KitaevComplex
