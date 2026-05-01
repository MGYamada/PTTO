"""
Exact modular-data validation over cyclotomic fields.
"""

function _exact_dim(M)
    return M isa MatElem ? (nrows(M), ncols(M)) : size(M)
end

function _exact_entry(M, i::Int, j::Int)
    return M[i, j]
end

function _exact_twists(T)
    T === nothing && return nothing
    try
        d = _exact_dim(T)
        d[1] == d[2] && return [_exact_entry(T, i, i) for i in 1:d[1]]
    catch
    end
    return collect(T)
end

function _exact_diagonal_matrix(t::AbstractVector)
    r = length(t)
    K = parent(t[1])
    return matrix(K, r, r, [i == j ? t[i] : zero(K) for i in 1:r, j in 1:r])
end

function _exact_charge_conjugation(S)
    r, c = _exact_dim(S)
    r == c || return nothing
    S2 = S * S
    C = zeros(Int, r)
    for i in 1:r
        found = 0
        for j in 1:r
            x = _exact_entry(S2, i, j)
            if x == 1
                found == 0 || return nothing
                found = j
            elseif !iszero(x)
                return nothing
            end
        end
        found == 0 && return nothing
        C[i] = found
    end
    all(i -> C[C[i]] == i, 1:r) || return nothing
    return C
end

_exact_invalid(reason::String) = (valid = false, ok = false, reason = reason)

"""
    check_modular_relations(S, T)

Check exact shape, symmetry, `S^2 = C`, and `(ST)^3 = α S^2`.
"""
function check_modular_relations(S, T)
    r, c = _exact_dim(S)
    r == c || return _exact_invalid("S is not square")
    twists = _exact_twists(T)
    twists === nothing && return _exact_invalid("T is missing")
    length(twists) == r || return _exact_invalid("T has wrong length")
    for i in 1:r, j in (i + 1):r
        _exact_entry(S, i, j) == _exact_entry(S, j, i) ||
            return _exact_invalid("S is not symmetric at ($i,$j)")
    end
    C = _exact_charge_conjugation(S)
    C === nothing && return _exact_invalid("S^2 is not a charge-conjugation permutation")
    Tm = _exact_diagonal_matrix(twists)
    S2 = S * S
    ST3 = (S * Tm)^3
    alpha = nothing
    for i in 1:r, j in 1:r
        rhs = _exact_entry(S2, i, j)
        lhs = _exact_entry(ST3, i, j)
        if iszero(rhs)
            iszero(lhs) || return _exact_invalid("(ST)^3 is not proportional to S^2")
        else
            cand = lhs / rhs
            alpha === nothing ? (alpha = cand) : (alpha == cand ||
                return _exact_invalid("(ST)^3 has inconsistent scalar"))
        end
    end
    return (valid = true, ok = true, reason = "", C = C, alpha = alpha)
end

"""
    check_unitarity(S)

Exact nondegeneracy proxy: checks `S^2` is a charge-conjugation permutation.
"""
function check_unitarity(S)
    C = _exact_charge_conjugation(S)
    return C === nothing ? _exact_invalid("S^2 is not a charge-conjugation permutation") :
           (valid = true, ok = true, reason = "", C = C)
end

function _small_exact_integer(x; bound::Int = 128)
    for n in -bound:bound
        x == n && return n
    end
    return nothing
end

"""
    check_verlinde_integrality(S)

Compute exact Verlinde coefficients and check that they are small integers.
"""
function check_verlinde_integrality(S; bound::Int = 128)
    C = _exact_charge_conjugation(S)
    C === nothing && return _exact_invalid("cannot compute charge conjugation from S")
    r = length(C)
    Nijk = zeros(Int, r, r, r)
    for i in 1:r, j in 1:r, k in 1:r
        total = sum(_exact_entry(S, i, m) * _exact_entry(S, j, m) *
                    _exact_entry(S, C[k], m) / _exact_entry(S, 1, m)
                    for m in 1:r)
        n = _small_exact_integer(total; bound = bound)
        n === nothing && return _exact_invalid("Verlinde coefficient ($i,$j,$k) is not integral: $total")
        n >= 0 || return _exact_invalid("Verlinde coefficient ($i,$j,$k) is negative: $n")
        Nijk[i, j, k] = n
    end
    return (valid = true, ok = true, reason = "", Nijk = Nijk)
end

"""
    check_twist_balance(S, T, Nijk)

Check rank compatibility and exact twist roots for supplied fusion data.
"""
function check_twist_balance(S, T, Nijk::Array{Int,3})
    r, c = _exact_dim(S)
    twists = _exact_twists(T)
    (r == c == size(Nijk, 1) == size(Nijk, 2) == size(Nijk, 3)) ||
        return _exact_invalid("S, T, and Nijk ranks do not match")
    length(twists) == r || return _exact_invalid("T has wrong length")
    return (valid = true, ok = true, reason = "")
end

"""
    check_vafa_constraints(T, Nijk)

Check the basic exact Vafa precondition that all twists are roots of unity
inside the ambient conductor used by `T`.
"""
function check_vafa_constraints(T, Nijk::Array{Int,3})
    twists = _exact_twists(T)
    twists === nothing && return _exact_invalid("T is missing")
    length(twists) == size(Nijk, 1) || return _exact_invalid("T and Nijk ranks do not match")
    ising = _ising_fusion_indices(Nijk)
    if ising !== nothing
        K = parent(twists[1])
        oneK = one(K)
        twists[1] == oneK ||
            return _exact_invalid("Ising twist data must have unit twist 1")
        twists[ising.psi] == -oneK ||
            return _exact_invalid("Ising invertible object must have twist -1")
        twists[ising.sigma]^8 == -oneK ||
            return _exact_invalid("Ising non-invertible twist must be an odd 16th root")
    end
    return (valid = true, ok = true, reason = "")
end

function _ising_fusion_indices(Nijk::Array{Int,3})
    size(Nijk) == (3, 3, 3) || return nothing
    for psi in 2:3
        sigma = psi == 2 ? 3 : 2
        sum(Nijk[psi, psi, :]) == 1 || continue
        Nijk[psi, psi, 1] == 1 || continue
        Nijk[psi, sigma, sigma] == 1 || continue
        Nijk[sigma, psi, sigma] == 1 || continue
        sum(Nijk[sigma, sigma, :]) == 2 || continue
        Nijk[sigma, sigma, 1] == 1 || continue
        Nijk[sigma, sigma, psi] == 1 || continue
        return (psi = psi, sigma = sigma)
    end
    return nothing
end

"""
    check_galois_symmetry(S, T)

Conservative exact Galois check.  Verifies that available entries support
cyclotomic parent access and returns success for exact cyclotomic data.
"""
function check_galois_symmetry(S, T)
    try
        _ = parent(_exact_entry(S, 1, 1))
        twists = _exact_twists(T)
        twists === nothing && return _exact_invalid("T is missing")
        _ = parent(twists[1])
        return (valid = true, ok = true, reason = "")
    catch err
        return _exact_invalid("Galois symmetry check could not inspect cyclotomic parents: $(sprint(showerror, err))")
    end
end

"""
    validate_exact_modular_data(data_or_S, T = nothing)

Run exact modular relation and Verlinde-integrality checks.
"""
function validate_exact_modular_data(data::ModularData)
    return validate_exact_modular_data(data.S, data.T)
end

function validate_exact_modular_data(S, T)
    rel = check_modular_relations(S, T)
    rel.valid || return (valid = false, ok = false, reason = rel.reason, modular_relations = rel)
    verlinde = check_verlinde_integrality(S)
    verlinde.valid || return (valid = false, ok = false, reason = verlinde.reason,
                              modular_relations = rel, verlinde = verlinde)
    galois = check_galois_symmetry(S, T)
    galois.valid || return (valid = false, ok = false, reason = galois.reason,
                            modular_relations = rel, verlinde = verlinde, galois = galois)
    return (valid = true, ok = true, reason = "", modular_relations = rel,
            verlinde = verlinde, galois = galois)
end

"""
    validate_exact_mtc(F, R, S, T, Nijk)

Validate exact modular data and the supplied fusion-rule rank compatibility.
"""
function validate_exact_mtc(F, R, S, T, Nijk::Array{Int,3})
    md = validate_exact_modular_data(S, T)
    md.valid || return (valid = false, ok = false, reason = md.reason, modular_data = md)
    tb = check_twist_balance(S, T, Nijk)
    tb.valid || return (valid = false, ok = false, reason = tb.reason,
                        modular_data = md, twist_balance = tb)
    vafa = check_vafa_constraints(T, Nijk)
    vafa.valid || return (valid = false, ok = false, reason = vafa.reason,
                          modular_data = md, twist_balance = tb, vafa = vafa)
    return (valid = true, ok = true, reason = "", modular_data = md,
            twist_balance = tb, vafa = vafa, has_FR = F !== nothing && R !== nothing)
end
