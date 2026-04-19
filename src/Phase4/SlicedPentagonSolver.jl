"""
    SlicedPentagonSolver

Combine pentagon equations with Kitaev χ³ F = 0 slice constraints.
Both live in the same Oscar polynomial ring, built via a single
polynomial-ring + polynomial-associator SixJCategory.

# Design

Since `ev` / `coev` on a SixJCategory with polynomial coefficients fails
(tries to invert a polynomial via `divexact`), we use the formula from
Kitaev Eq. 250 specialised to n = 3 in multiplicity-free form:

    (χ³ F)^{a1,a2}_y = (1/D²) Σ_c d_c · F^{c,a1,a2}_{y; a1, y}

and extract the polynomial by looking up `poly_C.ass[c,a1,a2,y][row,col]`
with `row = index_of(a1)` in rows-of-block and `col = index_of(y)` in
cols-of-block.

# Main API

- `get_sliced_pentagon_system(Nijk, r)`
     → (R, pent_eqs, slice_polys, n)
- `solve_pentagon_newton_with_slice(Nijk, r, base_F_func; …)`
     → Vector{Vector{ComplexF64}}
"""
module SlicedPentagonSolver

using LinearAlgebra
using SparseArrays
using KrylovKit
using Oscar
using TensorCategories
using ACMG: FusionRule
using ..KitaevComplex
using ..PentagonEquations
using ..PentagonSolver

export get_sliced_pentagon_system
export solve_pentagon_newton_with_slice

"""
    get_sliced_pentagon_system(Nijk::Array{Int,3}, r::Int)
        -> (R, pent_eqs, slice_polys, n)

Pentagon equations and Kitaev χ³ slice in ONE polynomial ring. Slice
is computed using the hand-derived formula

    (χ³ F)^{a1,a2}_y = (1/D²) Σ_c d_c · F^{c, a1, a2}_{y; a1, y}

avoiding ev/coev which can't be used on a polynomial-coefficient
SixJCategory.
"""
function get_sliced_pentagon_system(Nijk::Array{Int,3}, r::Int)
    one_vec = zeros(Int, r); one_vec[1] = 1

    _C = six_j_category(QQ, Nijk)
    _C.one = one_vec
    m = _C.simples
    var_count = TensorCategories._number_of_variables_in_pentagon_equations(_C)

    R, x = polynomial_ring(QQ, var_count)
    y_stack = deepcopy(x)

    poly_C = six_j_category(R, Nijk)
    poly_C.one = one_vec

    for i in 1:m, j in 1:m, k in 1:m, o in 1:m
        sum(poly_C.one[[i, j, k]]) > 0 && continue
        (rr, tt) = size(poly_C.ass[i, j, k, o])
        poly_C.ass[i, j, k, o] = matrix(R, rr, tt, [pop!(y_stack) for _ in 1:(rr*tt)])
    end

    # Pentagon equations (TensorCategories style)
    pent_eqs_raw = elem_type(R)[]
    for X in simples(poly_C), Y in simples(poly_C),
        Z in simples(poly_C), W in simples(poly_C)
        f = (id(X) ⊗ associator(Y, Z, W)) ∘
            associator(X, Y ⊗ Z, W) ∘
            (associator(X, Y, Z) ⊗ id(W))
        g = associator(X, Y, Z ⊗ W) ∘ associator(X ⊗ Y, Z, W)
        append!(pent_eqs_raw, collect(matrix(f - g))[:])
    end
    pent_eqs = filter(e -> e != 0 && !iszero(e), unique(pent_eqs_raw))

    # χ³ slice equations via direct associator entry access
    fpdims_qqbar = fpdim.(simples(_C))
    fpdims_f64 = [Float64(d) for d in fpdims_qqbar]
    D2_f64 = sum(fpdims_f64 .^ 2)

    # Helper to get F^{i,j,k}_{o; e, f} as polynomial in R
    function F_entry(i, j, k, o, e, f)
        (Nijk[i, j, e] ≥ 1 && Nijk[e, k, o] ≥ 1) || return nothing
        (Nijk[j, k, f] ≥ 1 && Nijk[i, f, o] ≥ 1) || return nothing

        rows_e = Int[]
        for ee in 1:m
            if Nijk[i, j, ee] ≥ 1 && Nijk[ee, k, o] ≥ 1
                push!(rows_e, ee)
            end
        end
        cols_f = Int[]
        for ff in 1:m
            if Nijk[j, k, ff] ≥ 1 && Nijk[i, ff, o] ≥ 1
                push!(cols_f, ff)
            end
        end
        row_idx = findfirst(==(e), rows_e)
        col_idx = findfirst(==(f), cols_f)
        (row_idx === nothing || col_idx === nothing) && return nothing

        if sum(one_vec[[i, j, k]]) > 0
            # Unit-fixed block: identity matrix
            return row_idx == col_idx ? one(R) : zero(R)
        end
        return poly_C.ass[i, j, k, o][row_idx, col_idx]
    end

    slice_polys = elem_type(R)[]
    for a1 in 1:m, a2 in 1:m, yv in 1:m
        Nijk[a1, a2, yv] ≥ 1 || continue
        chi = zero(R)
        for c in 1:m
            Fe = F_entry(c, a1, a2, yv, a1, yv)
            Fe === nothing && continue
            d_over_D2 = fpdims_f64[c] / D2_f64
            rat = Rational(rationalize(BigInt, d_over_D2; tol = 1e-12))
            chi += rat * Fe
        end
        if !iszero(chi)
            push!(slice_polys, chi)
        end
    end
    slice_polys = unique(slice_polys)

    return R, pent_eqs, slice_polys, var_count
end

"""
    solve_pentagon_newton_with_slice(Nijk, r, base_F_func; …)
        -> Vector{Vector{ComplexF64}}
"""
function solve_pentagon_newton_with_slice(Nijk::Array{Int,3}, r::Int,
        base_F_func::Function;
        initial_points::Union{Nothing, Vector{Vector{ComplexF64}}} = nothing,
        max_trials::Int = 5,
        max_iter::Int = 200,
        tol::Real = 1e-12,
        perturb_scale::Real = 0.05,
        verbose::Bool = false)

    R, pent_eqs, slice_polys, n = get_sliced_pentagon_system(Nijk, r)
    all_eqs = vcat(pent_eqs, slice_polys)
    derivs = [[derivative(eq, j) for j in 1:n] for eq in all_eqs]

    # Base F vector ordered via pop! traversal
    one_vec = zeros(Int, r); one_vec[1] = 1
    flat_keys = NTuple{6,Int}[]
    for i in 1:r, j in 1:r, k in 1:r, o in 1:r
        sum(one_vec[[i, j, k]]) > 0 && continue
        rows_e = Int[]; cols_f = Int[]
        for e in 1:r; if Nijk[i,j,e]≥1 && Nijk[e,k,o]≥1; push!(rows_e, e); end; end
        for f in 1:r; if Nijk[j,k,f]≥1 && Nijk[i,f,o]≥1; push!(cols_f, f); end; end
        (isempty(rows_e) || isempty(cols_f)) && continue
        for e in rows_e, f in cols_f
            push!(flat_keys, (i, j, k, o, e, f))
        end
    end
    F_base_vec = zeros(ComplexF64, n)
    for (pos, key) in enumerate(flat_keys)
        pidx = n - pos + 1
        F_base_vec[pidx] = ComplexF64(base_F_func(key...))
    end

    if initial_points === nothing
        initial_points = Vector{Vector{ComplexF64}}()
        push!(initial_points, copy(F_base_vec))
        for _ in 2:max_trials
            push!(initial_points,
                  F_base_vec + perturb_scale * randn(ComplexF64, n))
        end
    end

    solutions = Vector{Vector{ComplexF64}}()
    for (trial, x0) in enumerate(initial_points)
        x = copy(x0)
        for iter in 1:max_iter
            F_val = ComplexF64[PentagonSolver.eval_poly_complex(eq, x) for eq in all_eqs]
            res = maximum(abs.(F_val))
            if res < tol
                is_dup = any(norm(x - s) < 1e-8 for s in solutions)
                is_dup || push!(solutions, copy(x))
                verbose && println("  trial $trial: converged at iter $iter, res=$res")
                break
            end
            J = PentagonSolver.sparse_jacobian(all_eqs, derivs, x, n)
            delta, _ = linsolve(v -> J' * (J * v), J' * F_val;
                                ishermitian = true, isposdef = true, verbosity = 0)
            alpha = 1.0
            for _ in 1:20
                x_new = x - alpha * delta
                F_new = ComplexF64[PentagonSolver.eval_poly_complex(eq, x_new) for eq in all_eqs]
                if maximum(abs.(F_new)) < res
                    x = x_new
                    break
                end
                alpha *= 0.5
            end
        end
    end

    verbose && println("  Newton+slice: $(length(solutions)) unique solutions across $(length(initial_points)) trials")
    return solutions
end

end # module SlicedPentagonSolver
