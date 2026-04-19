"""
    Verify

Independent verification of Phase 4 solutions.

A Phase 4 solve returns numerical `(F, R)` that was found to satisfy
pentagon and hexagon equations up to some tolerance *by the solver*.
This module re-checks those conditions independently and exposes
additional MTC-level consistency checks:

- `pentagon_residuals(F, Nijk)`:   pentagon equation residuals given F
- `hexagon_residuals(F, R, Nijk)`: hexagon equation residuals given (F, R)
- `extract_R_block(R_values, mult, i, j, k)`: pull the R^{ij}_k block out
                                              of a flat R-variable vector
- `ribbon_residuals(R, T, Nijk, dual)`: check θ_i consistency with R
- `verify_mtc`: run all of the above and return a summary struct

All residuals are reported in the ∞-norm (`maximum(abs.(…))`).

This module is **purely downstream**: it does not depend on
HomotopyContinuation or KrylovKit. Only Oscar (for polynomial evaluation
using `eval_poly_complex`) and TensorCategories (for rebuilding a
category with numerical F and R to extract block structure).
"""
module Verify

using Oscar
using TensorCategories

# Reach into sibling modules
using ..PentagonEquations: get_pentagon_system
using ..PentagonSolver: eval_poly_complex
using ..HexagonEquations: get_hexagon_system, _number_of_variables_in_hexagon_equations

export pentagon_residuals, hexagon_residuals
export extract_R_block, block_positions_R
export ribbon_residuals, VerifyReport, verify_mtc

# ============================================================
#  Polynomial evaluation helpers (for AcbField coefficients)
# ============================================================

"""
    eval_poly_acb_at_complex(f, vals) -> ComplexF64

Evaluate an Oscar polynomial `f` with AcbField coefficients at a
ComplexF64 point. Hexagon equations have AcbField coefficients because
F_values were baked in as constants.
"""
function eval_poly_acb_at_complex(f, vals::Vector{ComplexF64})
    iszero(f) && return ComplexF64(0.0)
    result = 0.0 + 0.0im
    for (c, m) in zip(coefficients(f), monomials(f))
        # c :: AcbFieldElem -> ComplexF64 via midpoints
        cre = Float64(real(c))
        cim = Float64(imag(c))
        coef = ComplexF64(cre, cim)
        degs = degrees(m)
        term = coef
        for i in eachindex(degs)
            if degs[i] > 0
                term *= vals[i]^degs[i]
            end
        end
        result += term
    end
    return result
end

# ============================================================
#  Pentagon residuals
# ============================================================

"""
    pentagon_residuals(F_values, Nijk) -> Vector{Float64}

Evaluate every pentagon equation at the given F-values and return the
absolute residuals. For a correct solution all entries are ~ 1e-14.

Arguments:
- `F_values::Vector{ComplexF64}`: flat F-symbol values in the same order
                                  as produced by pentagon_equations (= Phase 4
                                  pentagon solution)
- `Nijk::Array{Int,3}`:           fusion coefficients
"""
function pentagon_residuals(F_values::Vector{ComplexF64}, Nijk::Array{Int,3})
    r = size(Nijk, 1)
    _R, eqs, n = get_pentagon_system(Nijk, r)
    length(F_values) == n || error("F_values has length $(length(F_values)), " *
                                    "expected $n variables for pentagon")
    return [abs(eval_poly_complex(eq, F_values)) for eq in eqs]
end

# ============================================================
#  Hexagon residuals
# ============================================================

"""
    hexagon_residuals(F_values, R_values, Nijk) -> Vector{Float64}

Evaluate every hexagon equation at the given (F, R) values and return
the absolute residuals.

Arguments:
- `F_values::Vector{ComplexF64}`: pentagon solution
- `R_values::Vector{ComplexF64}`: hexagon solution, length = 2 * r_var_count
                                  (forward and reverse braidings concatenated,
                                   as returned by solve_hexagon_homotopy)
- `Nijk::Array{Int,3}`:           fusion coefficients
"""
function hexagon_residuals(F_values::Vector{ComplexF64},
                           R_values::Vector{ComplexF64},
                           Nijk::Array{Int,3})
    r = size(Nijk, 1)
    _R_ring, eqs, n_r = get_hexagon_system(Nijk, r, F_values)
    length(R_values) == n_r || error("R_values has length $(length(R_values)), " *
                                      "expected $n_r variables for hexagon")
    return [abs(eval_poly_acb_at_complex(eq, R_values)) for eq in eqs]
end

# ============================================================
#  R-block extraction
# ============================================================

"""
    block_positions_R(Nijk) -> Dict{Tuple{Int,Int,Int}, Vector{Int}}

Return a dictionary mapping (i, j, k) with `Nijk[i, j, k] > 0` to the
list of positions (within the first half = forward-R section) occupied
by the entries of R^{ij}_k, in row-major order after the `pop!`
assignment order used by `_fill_braiding!` in HexagonEquations.

This mirrors the internal convention of hexagon_equations and is the
inverse map you need when you have a flat R_values and want R^{ij}_k[a, b].
"""
function block_positions_R(Nijk::Array{Int,3})
    m = size(Nijk, 1)

    # Count r_var_count the same way HexagonEquations does.
    r_var_count = 0
    for i in 1:m, j in 1:m, k in 1:m
        Nijk[i, j, k] > 0 && (r_var_count += Nijk[i, j, k]^2)
    end

    positions = Dict{Tuple{Int,Int,Int}, Vector{Int}}()
    y_positions = collect(1:r_var_count)
    for i in 1:m, j in 1:m, k in 1:m
        N_ijk = Nijk[i, j, k]
        N_ijk == 0 && continue
        popped = Int[]
        for _ in 1:(N_ijk * N_ijk)
            push!(popped, pop!(y_positions))
        end
        positions[(i, j, k)] = popped
    end
    @assert isempty(y_positions)
    return positions
end

"""
    extract_R_block(R_values, Nijk, i, j, k) -> Matrix{ComplexF64}

Return the R^{ij}_k block from a flat R-values vector. `R_values` is
expected to be the forward-braiding values (the first half of what
`solve_hexagon_homotopy` returns) OR the full length-`2·r_var_count`
vector — in the latter case only the first half is used.
"""
function extract_R_block(R_values::Vector{ComplexF64}, Nijk::Array{Int,3},
                         i::Int, j::Int, k::Int)
    N_ijk = Nijk[i, j, k]
    N_ijk == 0 && return Matrix{ComplexF64}(undef, 0, 0)

    positions = block_positions_R(Nijk)
    haskey(positions, (i, j, k)) || error("($i, $j, $k) has no R-block in the Nijk")

    pos = positions[(i, j, k)]
    M = Matrix{ComplexF64}(undef, N_ijk, N_ijk)
    # In HexagonEquations `entries[p]` fills matrix row-major, so
    # entries[(a-1)*N_ijk + b] = matrix[a, b].
    for a in 1:N_ijk, b in 1:N_ijk
        M[a, b] = R_values[pos[(a - 1) * N_ijk + b]]
    end
    return M
end

# ============================================================
#  Ribbon (twist) consistency
# ============================================================

"""
    ribbon_residuals(R_values, T_expected, Nijk) -> Vector{Float64}

For every triple (i, j, k) with `Nijk[i, j, k] > 0`, check the
multiplicity-free ribbon relation

    (R^{ij}_k)² = θ_i · θ_j / θ_k

This follows from the ribbon axiom c_{Y,X} ∘ c_{X,Y} = θ_{X⊗Y} · (θ_X⁻¹ ⊗ θ_Y⁻¹),
restricted to the X_k summand of X_i ⊗ X_j. In the multiplicity-free
case each R-block is a scalar, so the relation becomes a single
quadratic equation per non-zero fusion entry.

This is **independent of the ± sign of R** (squared), so both of
Fibonacci's two braidings must satisfy it.

Arguments:
- `R_values::Vector{ComplexF64}`: hexagon solution (full or forward-half)
- `T_expected::Vector{ComplexF64}`: θ_i values (θ_k on RHS, θ_i θ_j on LHS)
- `Nijk::Array{Int,3}`:             fusion

Returns: `Vector{Float64}` one residual per non-zero (i, j, k) triple,
each entry |(R^{ij}_k)² · θ_k - θ_i · θ_j|.

Non-multiplicity-free case is NOT supported (Nijk[i,j,k] > 1 will
error).

Caveat: This is a NECESSARY condition only. It does not check that
θ_i matches the braiding (same formula holds for -θ_i too because of
the square).
"""
function ribbon_residuals(R_values::Vector{ComplexF64},
                          T_expected::Vector{ComplexF64},
                          Nijk::Array{Int,3})
    r = size(Nijk, 1)
    length(T_expected) == r || error("T_expected has wrong length")

    # Expected r_var_count for forward-only
    r_var_count = 0
    for i in 1:r, j in 1:r, k in 1:r
        Nijk[i, j, k] > 0 && (r_var_count += Nijk[i, j, k]^2)
    end
    # Accept either forward-only or full length
    if length(R_values) == 2 * r_var_count
        R_fwd = R_values[1:r_var_count]
    elseif length(R_values) == r_var_count
        R_fwd = R_values
    else
        error("R_values has length $(length(R_values)); expected " *
              "$r_var_count (forward only) or $(2 * r_var_count) (both).")
    end

    residuals = Float64[]
    for i in 1:r, j in 1:r, k in 1:r
        N_ijk = Nijk[i, j, k]
        N_ijk == 0 && continue
        N_ijk == 1 || error("Non-multiplicity-free (N[$i,$j,$k]=$N_ijk) " *
                             "not supported by ribbon_residuals")
        R_block = extract_R_block(R_fwd, Nijk, i, j, k)
        R = R_block[1, 1]
        lhs = R^2 * T_expected[k]
        rhs = T_expected[i] * T_expected[j]
        push!(residuals, abs(lhs - rhs))
    end
    return residuals
end

# ============================================================
#  Summary report
# ============================================================

"""
    VerifyReport

Summary struct for verification output.
"""
struct VerifyReport
    pentagon_max::Float64
    hexagon_max::Float64
    ribbon_max::Union{Float64, Nothing}
    n_pentagon_eqs::Int
    n_hexagon_eqs::Int
    rank::Int
end

function Base.show(io::IO, r::VerifyReport)
    ribstr = r.ribbon_max === nothing ? "n/a" : string(r.ribbon_max)
    print(io, "VerifyReport(rank=$(r.rank), ",
          "pentagon_max=$(r.pentagon_max) over $(r.n_pentagon_eqs) eqs, ",
          "hexagon_max=$(r.hexagon_max) over $(r.n_hexagon_eqs) eqs, ",
          "ribbon_max=$ribstr)")
end

"""
    verify_mtc(F_values, R_values, Nijk; T=nothing) -> VerifyReport

Run all three residual checks and return a summary.

If `T` is supplied, the ribbon check ((R^{ij}_k)² = θ_i θ_j / θ_k) is
included. Otherwise it is skipped (ribbon_max = nothing).
"""
function verify_mtc(F_values::Vector{ComplexF64},
                    R_values::Vector{ComplexF64},
                    Nijk::Array{Int,3};
                    T::Union{Vector{ComplexF64}, Nothing} = nothing)
    pent = pentagon_residuals(F_values, Nijk)
    hex = hexagon_residuals(F_values, R_values, Nijk)

    ribbon_max = nothing
    if T !== nothing
        rib = ribbon_residuals(R_values, T, Nijk)
        ribbon_max = maximum(rib)
    end

    return VerifyReport(
        maximum(pent),
        maximum(hex),
        ribbon_max,
        length(pent),
        length(hex),
        size(Nijk, 1),
    )
end

end # module Verify
