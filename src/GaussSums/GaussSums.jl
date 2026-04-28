"""
Exact Gauss sums for cyclotomic modular data.

The convention is `d_i = S[1,i] / S[1,1]` and
`τ_n^+ = sum_i θ_i^n d_i^2`, where `θ_i` is the diagonal entry of `T`.
All arithmetic stays in the cyclotomic parent carried by the input.
"""

struct _GaussSumInput
    context::CyclotomicContext
    S::Any
    T::Any
end

_gauss_input(data::ModularData) = _GaussSumInput(data.context, data.S, data.T)

function _gauss_input(result::ClassifiedMTC)
    return _GaussSumInput(CyclotomicContext(result.N), result.S_cyclotomic, result.T_cyclotomic)
end

function _gauss_input(x)
    error("expected ModularData or ClassifiedMTC with exact cyclotomic S/T data, got $(typeof(x))")
end

_matrix_nrows(M::MatElem) = nrows(M)
_matrix_nrows(M::AbstractMatrix) = size(M, 1)

function _twist_at(T::MatElem, i::Int)
    return T[i, i]
end

function _twist_at(T::AbstractMatrix, i::Int)
    return T[i, i]
end

function _twist_at(T::AbstractVector, i::Int)
    return T[i]
end

function _rank_of(input::_GaussSumInput)
    return _matrix_nrows(input.S)
end

function _check_gauss_input(input::_GaussSumInput)
    r = _rank_of(input)
    r >= 1 || error("modular data must have positive rank")
    input.S[1, 1] != zero(input.context.field) || error("S[1,1] must be nonzero")
    return r
end

function quantum_dimensions(category_or_modular_data)
    input = _gauss_input(category_or_modular_data)
    r = _check_gauss_input(input)
    s00 = input.S[1, 1]
    return [input.S[1, i] / s00 for i in 1:r]
end

"""
    total_quantum_dimension_squared(category_or_modular_data)

Return `D^2 = sum_i d_i^2` exactly in the input cyclotomic field.
"""
function total_quantum_dimension_squared(category_or_modular_data)
    input = _gauss_input(category_or_modular_data)
    _check_gauss_input(input)
    K = input.context.field
    return sum(d^2 for d in quantum_dimensions(category_or_modular_data); init = zero(K))
end

function _gauss_sum(category_or_modular_data; n::Integer = 1, sign::Int = 1)
    input = _gauss_input(category_or_modular_data)
    r = _check_gauss_input(input)
    ds = quantum_dimensions(category_or_modular_data)
    K = input.context.field
    total = zero(K)
    for i in 1:r
        θ = _twist_at(input.T, i)
        total += θ^(sign * Int(n)) * ds[i]^2
    end
    return total
end

"""
    gauss_sum_plus(category_or_modular_data; n = 1)

Return the exact higher positive Gauss sum `τ_n^+ = sum_i θ_i^n d_i^2`.
"""
gauss_sum_plus(category_or_modular_data; n::Integer = 1) =
    _gauss_sum(category_or_modular_data; n = n, sign = 1)

"""
    gauss_sum_minus(category_or_modular_data; n = 1)

Return the exact higher negative Gauss sum `τ_n^- = sum_i θ_i^(-n) d_i^2`.
"""
gauss_sum_minus(category_or_modular_data; n::Integer = 1) =
    _gauss_sum(category_or_modular_data; n = n, sign = -1)

"""
    normalized_gauss_sum(category_or_modular_data; n = 1, normalization = :D)

Normalize `τ_n^+`.

- `:raw`: return `τ_n^+`
- `:D`: return `τ_n^+ / D`, with `D = 1 / S[1,1]`
- `:D2`: return `τ_n^+ / D^2`
"""
function normalized_gauss_sum(category_or_modular_data;
                              n::Integer = 1,
                              normalization::Symbol = :D)
    input = _gauss_input(category_or_modular_data)
    _check_gauss_input(input)
    τ = gauss_sum_plus(category_or_modular_data; n = n)
    normalization == :raw && return τ
    D = inv(input.S[1, 1])
    normalization == :D && return τ / D
    normalization == :D2 && return τ / D^2
    error("unknown Gauss-sum normalization: $(normalization); expected :raw, :D, or :D2")
end
