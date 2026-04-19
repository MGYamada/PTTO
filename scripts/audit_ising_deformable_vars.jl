"""
scripts/audit_ising_deformable_vars.jl

List the 11 deformable F-entries that KitaevComplex identifies for Ising,
and check for each one whether it's actually a free variable or forced by
fusion-path uniqueness.

A 1×1 block F^{abc}_{d;e,f} has only one allowed (e,f) pair, so the entry
is structurally 1-dimensional. Its value is determined by pentagon alone;
there is no further gauge freedom within it. Such entries are still F-vars
in the tangent space formally, but KitaevComplex's Δ_gauge should show
them as gauge-fixed (column-zero) if the unit-preserving subspace is
correctly constructed.

If Δ_gauge produces NON-zero rows for entries that are structurally fixed
to ±1 by the pentagon equations alone, then analyze_gauge may be
under-counting gauge directions (the "gauge" it detects is only the
partial Hadamard-block rotation, missing other structural freedoms).
"""

using LinearAlgebra
using ACMG
using ACMG.Phase4

const KC  = Phase4.KitaevComplex
const SPS = Phase4.SlicedPentagonSolver

function ising_Nijk()
    r = 3
    N = zeros(Int, r, r, r)
    for j in 1:r
        N[1, j, j] = 1
        N[j, 1, j] = 1
    end
    N[2, 2, 1] = 1
    N[2, 3, 3] = 1
    N[3, 2, 3] = 1
    N[3, 3, 1] = 1
    N[3, 3, 2] = 1
    return N
end

function ising_F_func()
    F_sss_s = (1.0 / sqrt(2.0)) * [1.0  1.0;
                                    1.0 -1.0]
    fr = FusionRule(ising_Nijk())
    function F(a::Int, b::Int, c::Int, d::Int, e::Int, f::Int)
        (fr.N[a, b, e] == 1 && fr.N[e, c, d] == 1 &&
         fr.N[b, c, f] == 1 && fr.N[a, f, d] == 1) || return 0.0
        1 in (a, b, c, d) && return 1.0
        nsig = count(==(3), (a, b, c, d))
        nsig == 0 && return 1.0
        if nsig == 2
            (a, b, c, d) == (2, 3, 2, 3) && return -1.0
            (a, b, c, d) == (3, 2, 3, 1) && return  1.0
            (a, b, c, d) == (3, 2, 3, 2) && return -1.0
            return 1.0
        end
        if nsig == 4
            d != 3 && return 0.0
            return F_sss_s[e, f]
        end
        return 1.0
    end
    return F
end

function main()
    Nijk = ising_Nijk()
    F_fn = ising_F_func()
    r = 3
    fr = FusionRule(Nijk)

    fcs = KC.build_F_coord_space(fr, F_fn)
    ga  = KC.analyze_gauge(fcs)

    println("KitaevComplex analysis:")
    println("  Total F-vars: $(KC.F_var_count(fcs))")
    println("  Unit-fixed (any of abcd = 1): $(length(fcs.unit_indices))")
    println("  Gauge orbit dim: $(ga.gauge_orbit_dim)")

    # List the deformable (non-unit-fixed) entries with block structure
    println("\nDeformable (non-unit-fixed) entries, with block size:")
    println(rpad("idx", 4), rpad("FKey", 28), rpad("value", 22), "block_shape")
    deformable = [(k, i) for (i, k) in enumerate(fcs.vars) if !(1 in k[1:4])]

    # Compute block shape: for each (a,b,c,d), count allowed (e,f) pairs
    function block_shape(a, b, c, d, fr)
        rows = 0
        for e in 1:fr.rank
            if fr.N[a, b, e] ≥ 1 && fr.N[e, c, d] ≥ 1
                rows += 1
            end
        end
        cols = 0
        for f in 1:fr.rank
            if fr.N[b, c, f] ≥ 1 && fr.N[a, f, d] ≥ 1
                cols += 1
            end
        end
        return (rows, cols)
    end

    for (key, idx) in deformable
        a, b, c, d, e, f = key
        shp = block_shape(a, b, c, d, fr)
        val = F_fn(key...)
        println(rpad(idx, 4), rpad(string(key), 28), rpad(round(val, digits=4), 22), shp)
    end

    # Δ_gauge rows: which deformable entries are actually in the gauge orbit image?
    println("\nΔ_gauge_eff rows (rows = F-entries):")
    Dg_eff = ga.Delta_gauge_eff
    println("  shape: $(size(Dg_eff))")
    # Normalize each deformable row's contribution
    println("  Norms of Δ_gauge_eff rows per deformable entry:")
    for (key, idx) in deformable
        a, b, c, d, e, f = key
        row_norm = norm(Dg_eff[idx, :])
        marker = row_norm > 1e-10 ? "(in gauge orbit)" : "(not in gauge)"
        shp = block_shape(a, b, c, d, fr)
        print(rpad(idx, 4), rpad(string(key), 28), rpad(shp, 10))
        println(rpad(round(row_norm, digits=6), 15), "  ", marker)
    end

    # The key question: are 1×1-block deformable entries really free?
    # For each 1×1-block entry, the pentagon alone fixes the value. If we
    # solve pentagon at this entry set to a DIFFERENT value, we should get
    # residual ≠ 0. Let's check by perturbing.
    println("\n" * "="^68)
    println("1×1-block entries: perturb each and measure pentagon residual")
    println("="^68)

    n = KC.F_var_count(fcs)
    one_vec = [1, 0, 0]
    fkey_map = SPS.build_fkey_to_xvar_map(Nijk, r, one_vec)
    F_base = zeros(ComplexF64, 14)
    for (key, pidx) in fkey_map
        F_base[pidx] = ComplexF64(F_fn(key...))
    end

    # Hand-coded pentagon residual for comparison
    function F_of(Fvec, a, b, c, d, e, f)
        (fr.N[a, b, e] ≥ 1 && fr.N[e, c, d] ≥ 1 &&
         fr.N[b, c, f] ≥ 1 && fr.N[a, f, d] ≥ 1) || return ComplexF64(0)
        1 in (a, b, c, d) && return ComplexF64(1)
        key = (a, b, c, d, e, f)
        haskey(fkey_map, key) || return ComplexF64(0)
        return Fvec[fkey_map[key]]
    end

    function pent_max_resid(Fvec)
        m = 0.0
        for a in 1:r, b in 1:r, c in 1:r, d in 1:r, e in 1:r,
            f in 1:r, g in 1:r, k in 1:r, l in 1:r
            lhs1 = F_of(Fvec, f, c, d, e, g, l);  lhs1 == 0 && continue
            lhs2 = F_of(Fvec, a, b, l, e, f, k);  lhs2 == 0 && continue
            lhs = lhs1 * lhs2
            rhs = ComplexF64(0)
            for h in 1:r
                x1 = F_of(Fvec, a, b, c, g, f, h); x1 == 0 && continue
                x2 = F_of(Fvec, a, h, d, e, g, k); x2 == 0 && continue
                x3 = F_of(Fvec, b, c, d, k, h, l); x3 == 0 && continue
                rhs += x1 * x2 * x3
            end
            d_res = abs(lhs - rhs)
            d_res > m && (m = d_res)
        end
        return m
    end

    println("\nBase F pentagon residual (sanity): $(pent_max_resid(F_base))")

    for (key, _) in deformable
        a, b, c, d, e, f = key
        shp = block_shape(a, b, c, d, fr)
        if shp == (1, 1)
            pidx = fkey_map[key]
            # Perturb just this coordinate
            Fpert = copy(F_base)
            Fpert[pidx] += 0.1
            res = pent_max_resid(Fpert)
            println("  $(key)  [1×1]  perturb +0.1:  pent residual = $(round(res, digits=5))")
        end
    end
end

main()
