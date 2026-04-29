"""
    HexagonEquations

Exact hexagon equation generation over a cyclotomic field.
"""

using Oscar
using TensorCategories

function _coerce_exact(c, K)
    parent(c) === K && return c
    return K(c)
end

function assign_F_to_associator!(poly_C, F_values::Vector)
    K = base_ring(poly_C)
    m = poly_C.simples
    one_vec = poly_C.one
    y = copy(F_values)
    for i in 1:m, j in 1:m, k in 1:m, o in 1:m
        sum(one_vec[[i, j, k]]) > 0 && continue
        r, t = size(poly_C.ass[i, j, k, o])
        entries = [_coerce_exact(pop!(y), K) for _ in 1:(r * t)]
        poly_C.ass[i, j, k, o] = matrix(K, r, t, entries)
    end
    isempty(y) || error("F_values length does not match associator slots")
    return poly_C
end

function _number_of_variables_in_hexagon_equations(poly_C)
    m = poly_C.simples
    mult = poly_C.tensor_product
    n_r = 0
    for i in 1:m, j in 1:m, k in 1:m
        N_ijk = mult[i, j, k]
        n_r += N_ijk * N_ijk
    end
    return n_r
end

number_of_variables_in_hexagon_equations(poly_C) =
    _number_of_variables_in_hexagon_equations(poly_C)

function _braiding_block_positions(mult::Array{Int,3})
    m = size(mult, 1)
    total = sum(mult[i, j, k]^2 for i in 1:m, j in 1:m, k in 1:m)
    positions = Dict{Tuple{Int,Int,Int}, Vector{Int}}()
    y = collect(1:total)
    for i in 1:m, j in 1:m, k in 1:m
        N_ijk = mult[i, j, k]
        N_ijk == 0 && continue
        popped = Int[]
        for _ in 1:(N_ijk * N_ijk)
            push!(popped, pop!(y))
        end
        positions[(i, j, k)] = popped
    end
    return positions, total
end

function _fill_braiding!(poly_C, vars, mult::Array{Int,3}, R_ring)
    m = poly_C.simples
    braid_arr = Array{MatElem, 3}(undef, m, m, m)
    y = copy(vars)
    for i in 1:m, j in 1:m, k in 1:m
        N_ijk = mult[i, j, k]
        if N_ijk == 0
            braid_arr[i, j, k] = zero_matrix(R_ring, 0, 0)
        else
            entries = [pop!(y) for _ in 1:(N_ijk * N_ijk)]
            braid_arr[i, j, k] = matrix(R_ring, N_ijk, N_ijk, entries)
        end
    end
    TensorCategories.set_braiding!(poly_C, braid_arr)
    return poly_C
end

function _hexagon_equations_with_F_values(mult::Array{Int,3}, one_vec::Vector{Int},
                                          F_values::Vector; context::CyclotomicContext)
    K = field(context)
    r = size(mult, 1)
    dummy = TensorCategories.six_j_category(K, mult)
    dummy.one = one_vec
    r_var_count = _number_of_variables_in_hexagon_equations(dummy)

    R_ring, xs = polynomial_ring(K, 2 * r_var_count, :r)
    r_vars = xs[1:r_var_count]
    s_vars = xs[(r_var_count + 1):(2 * r_var_count)]

    poly_C_fwd = TensorCategories.six_j_category(R_ring, mult)
    poly_C_fwd.one = one_vec
    assign_F_to_associator!(poly_C_fwd, F_values)
    _fill_braiding!(poly_C_fwd, r_vars, mult, R_ring)

    eqs = elem_type(R_ring)[]
    Ss = simples(poly_C_fwd)
    for X in Ss, Y in Ss, Z in Ss
        lhs = (braiding(X, Y) ⊗ id(Z)) ∘ associator(Y, X, Z) ∘ (id(Y) ⊗ braiding(X, Z))
        rhs = associator(Y, Z, X) ∘ braiding(X, Y ⊗ Z) ∘ associator(X, Y, Z)
        append!(eqs, collect(matrix(lhs - rhs))[:])
    end

    r_pos, _ = _braiding_block_positions(mult)
    s_pos = r_pos
    for i in 1:r, j in 1:r, k in 1:r
        N_ijk = mult[i, j, k]
        N_ijk == 0 && continue
        haskey(s_pos, (j, i, k)) || continue
        R_block = [r_vars[r_pos[(i, j, k)][(a - 1) * N_ijk + b]]
                   for a in 1:N_ijk, b in 1:N_ijk]
        S_block = [s_vars[s_pos[(j, i, k)][(a - 1) * N_ijk + b]]
                   for a in 1:N_ijk, b in 1:N_ijk]
        for a in 1:N_ijk, c in 1:N_ijk
            target = a == c ? one(R_ring) : zero(R_ring)
            push!(eqs, sum(R_block[a, b] * S_block[b, c] for b in 1:N_ijk) - target)
        end
    end

    # Unit braiding normalization: c_{1,X} and c_{X,1}, together with the
    # explicit inverse variables, are identity on every nonzero channel.
    for i in 1:r, j in 1:r, k in 1:r
        (i == 1 || j == 1) || continue
        N_ijk = mult[i, j, k]
        N_ijk == 0 && continue
        for a in 1:N_ijk, b in 1:N_ijk
            target = a == b ? one(R_ring) : zero(R_ring)
            pos = r_pos[(i, j, k)][(a - 1) * N_ijk + b]
            push!(eqs, r_vars[pos] - target)
            push!(eqs, s_vars[pos] - target)
        end
    end

    return poly_C_fwd, nothing, filter(e -> !iszero(e), unique(eqs))
end

function _hexagon_equations_with_symbolic_FR(mult::Array{Int,3}, one_vec::Vector{Int};
                                             context = nothing)
    r = size(mult, 1)
    K = context === nothing ? QQ : field(context)
    _, _, f_var_count = get_pentagon_system(mult, r)

    dummy = TensorCategories.six_j_category(K, mult)
    dummy.one = one_vec
    r_var_count = _number_of_variables_in_hexagon_equations(dummy)

    H, xs = polynomial_ring(K, f_var_count + 2 * r_var_count, :h)
    f_vars = xs[1:f_var_count]
    r_vars = xs[(f_var_count + 1):(f_var_count + r_var_count)]
    s_vars = xs[(f_var_count + r_var_count + 1):(f_var_count + 2 * r_var_count)]

    poly_C = TensorCategories.six_j_category(H, mult)
    poly_C.one = one_vec
    assign_F_to_associator!(poly_C, collect(f_vars))
    _fill_braiding!(poly_C, r_vars, mult, H)

    eqs = elem_type(H)[]
    Ss = simples(poly_C)
    for X in Ss, Y in Ss, Z in Ss
        lhs = (braiding(X, Y) ⊗ id(Z)) ∘ associator(Y, X, Z) ∘ (id(Y) ⊗ braiding(X, Z))
        rhs = associator(Y, Z, X) ∘ braiding(X, Y ⊗ Z) ∘ associator(X, Y, Z)
        append!(eqs, collect(matrix(lhs - rhs))[:])
    end

    r_pos, _ = _braiding_block_positions(mult)
    s_pos = r_pos
    for i in 1:r, j in 1:r, k in 1:r
        N_ijk = mult[i, j, k]
        N_ijk == 0 && continue
        haskey(s_pos, (j, i, k)) || continue
        R_block = [r_vars[r_pos[(i, j, k)][(a - 1) * N_ijk + b]]
                   for a in 1:N_ijk, b in 1:N_ijk]
        S_block = [s_vars[s_pos[(j, i, k)][(a - 1) * N_ijk + b]]
                   for a in 1:N_ijk, b in 1:N_ijk]
        for a in 1:N_ijk, c in 1:N_ijk
            target = a == c ? one(H) : zero(H)
            push!(eqs, sum(R_block[a, b] * S_block[b, c] for b in 1:N_ijk) - target)
        end
    end

    for i in 1:r, j in 1:r, k in 1:r
        (i == 1 || j == 1) || continue
        N_ijk = mult[i, j, k]
        N_ijk == 0 && continue
        for a in 1:N_ijk, b in 1:N_ijk
            target = a == b ? one(H) : zero(H)
            pos = r_pos[(i, j, k)][(a - 1) * N_ijk + b]
            push!(eqs, r_vars[pos] - target)
            push!(eqs, s_vars[pos] - target)
        end
    end

    eqs = filter(e -> !iszero(e), unique(eqs))
    return H, eqs, f_var_count, r_var_count
end

function get_hexagon_fr_system(Nijk::Array{Int,3}, r::Int; context = nothing)
    size(Nijk) == (r, r, r) || error("Nijk must have shape ($r, $r, $r), got $(size(Nijk))")
    one_vec = zeros(Int, r)
    one_vec[1] = 1
    H, eqs, nF, nR = _hexagon_equations_with_symbolic_FR(Nijk, one_vec; context = context)
    isempty(eqs) && error("Hexagon produced no equations")
    return H, eqs, nF + 2 * nR
end

function get_hexagon_system(Nijk::Array{Int,3}, r::Int, F_values::Vector;
                            context::CyclotomicContext)
    one_vec = zeros(Int, r)
    one_vec[1] = 1
    _, _, eqs = _hexagon_equations_with_F_values(Nijk, one_vec, F_values; context = context)
    isempty(eqs) && error("Hexagon produced no equations")
    R_ring = parent(eqs[1])
    return R_ring, eqs, nvars(R_ring)
end
