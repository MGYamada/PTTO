"""
    KitaevComplex

Kitaev 2006 Eq. 250 contracting homotopy `χ^n : C^n → C^{n-1}` for n = 3,
implemented by accessing the polynomial associator entries directly.

# χ³ formula (hand-derived from Kitaev Eq. 250, multiplicity-free)

    (χ³ F)^{a1, a2}_y = (1/D²) Σ_c d_c · F^{c, a1, a2}_{y; a1, y}

valid when the associator entry `F^{c, a1, a2}_{y; a1, y}` exists, i.e.
all four fusion rules are satisfied:

    N^{c a1}_{a1} ≥ 1        (row index e = a1 ∈ c ⊗ a1)
    N^{a1 a2}_{y} ≥ 1        (col index f = y ∈ a1 ⊗ a2)
    N^{a1 a2}_{y} ≥ 1        (rows exit into overall charge y via a2)
    N^{c y}_{y} ≥ 1          (cols exit into overall charge y via y ⊗ c = y)

This form avoids the use of `ev` / `coev` on polynomial-coefficient
SixJCategory objects, which trigger `divexact` by polynomial and fail.

# Output

Same style as `pentagon_equations`: returns `(poly_C, eqs)`. Eqs are
linear polynomials in the pentagon variables.
"""
module KitaevComplex

using Oscar
using TensorCategories
using ACMG: FusionRule

export chi3_equations

function chi3_equations(mult::Array{Int,3}, one_vec::Vector{Int})
    _C = six_j_category(QQ, mult)
    _C.one = one_vec
    m = _C.simples

    var_count = TensorCategories._number_of_variables_in_pentagon_equations(_C)

    R, x = polynomial_ring(QQ, var_count)
    y_stack = deepcopy(x)

    poly_C = six_j_category(R, mult)
    poly_C.one = one_vec

    # Fill associator blocks with polynomial variables using the same
    # traversal as pentagon_equations / assign_F_to_associator!.
    for i in 1:m, j in 1:m, k in 1:m, o in 1:m
        sum(poly_C.one[[i, j, k]]) > 0 && continue
        (rr, tt) = size(poly_C.ass[i, j, k, o])
        poly_C.ass[i, j, k, o] = matrix(R, rr, tt, [pop!(y_stack) for _ in 1:(rr*tt)])
    end

    # Quantum dims via fpdim (QQBar), convert to rational.
    fpdims_qqbar = fpdim.(simples(_C))
    fpdims_f64 = [Float64(d) for d in fpdims_qqbar]
    D2_f64 = sum(fpdims_f64 .^ 2)

    # Helper: look up the F-entry F^{i,j,k}_{o; e, f} as a polynomial in R.
    # Based on how TensorCategories indexes its `ass` blocks.
    function F_entry(i, j, k, o, e, f)
        # Is block size nonzero?
        # Basis of rows: e with N^{ij}_e ≥ 1 and N^{ek}_o ≥ 1.
        # Basis of cols: f with N^{jk}_f ≥ 1 and N^{if}_o ≥ 1.
        (mult[i, j, e] ≥ 1 && mult[e, k, o] ≥ 1) || return nothing
        (mult[j, k, f] ≥ 1 && mult[i, f, o] ≥ 1) || return nothing

        # Find row index of e in the block
        rows_e = Int[]
        for ee in 1:m
            if mult[i, j, ee] ≥ 1 && mult[ee, k, o] ≥ 1
                push!(rows_e, ee)
            end
        end
        cols_f = Int[]
        for ff in 1:m
            if mult[j, k, ff] ≥ 1 && mult[i, ff, o] ≥ 1
                push!(cols_f, ff)
            end
        end

        row_idx = findfirst(==(e), rows_e)
        col_idx = findfirst(==(f), cols_f)
        (row_idx === nothing || col_idx === nothing) && return nothing

        # Is the block unit-axiom-fixed (i, j, or k is the unit)?
        # In that case TensorCategories sets it to identity matrix (1 on diagonal).
        if sum(one_vec[[i, j, k]]) > 0
            # Unit-fixed block: F is identity → F[row,col] = 1 if row==col else 0
            return row_idx == col_idx ? one(R) : zero(R)
        end

        return poly_C.ass[i, j, k, o][row_idx, col_idx]
    end

    eqs = elem_type(R)[]

    # For each (a1, a2) ∈ {1..m}² with y in a1 ⊗ a2, build slice eq.
    for a1 in 1:m, a2 in 1:m, y in 1:m
        mult[a1, a2, y] ≥ 1 || continue
        # χ³ should vanish to give the slice constraint.
        chi_expr = zero(R)
        for c in 1:m
            Fe = F_entry(c, a1, a2, y, a1, y)
            Fe === nothing && continue
            d_c_over_D2 = fpdims_f64[c] / D2_f64
            rat_coef = Rational(rationalize(BigInt, d_c_over_D2; tol = 1e-12))
            chi_expr += rat_coef * Fe
        end
        if !iszero(chi_expr)
            push!(eqs, chi_expr)
        end
    end

    return poly_C, unique(eqs)
end

end # module KitaevComplex
