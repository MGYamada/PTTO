"""
Higher central charges from exact cyclotomic Gauss sums.
"""

struct HigherCentralChargeResult
    ok::Bool
    n::Int
    normalization::Symbol
    value::Any
    gauss_sum::Any
    denominator::Any
    D_squared::Any
    conductor::Union{Int, Nothing}
    status::Symbol
    message::String
end

function Base.show(io::IO, r::HigherCentralChargeResult)
    if r.ok
        print(io, "HigherCentralChargeResult(n=$(r.n), normalization=$(r.normalization), value=$(r.value))")
    else
        print(io, "HigherCentralChargeResult(n=$(r.n), status=$(r.status), message=$(repr(r.message)))")
    end
end

function _higher_failure(category_or_modular_data, n::Integer, normalization::Symbol,
                         status::Symbol, message::String)
    input = try
        _gauss_input(category_or_modular_data)
    catch
        nothing
    end
    conductor_value = input === nothing ? nothing : conductor(input.context)
    return HigherCentralChargeResult(false, Int(n), normalization, nothing, nothing, nothing,
                                     nothing, conductor_value, status, message)
end

"""
    higher_central_charge(category_or_modular_data; n = 1, normalization = :galois)

Return a `HigherCentralChargeResult`.

For `normalization = :galois`, this computes `τ_n^+ / σ_n(D)` and requires
`n` to be a unit modulo the cyclotomic conductor.  Use `normalization = :D`
for the raw `τ_n^+ / D` ratio at arbitrary integer `n`, or `:D2` for
`τ_n^+ / D^2`.
"""
function higher_central_charge(category_or_modular_data;
                               n::Integer = 1,
                               normalization::Symbol = :galois,
                               method::Symbol = :modular_data,
                               p::Union{Integer, Nothing} = nothing,
                               kwargs...)
    if method == :finite_field
        p === nothing && error("finite-field higher central charge requires keyword p")
        solution = solve_fr_mod_p(category_or_modular_data, Int(p); kwargs...)
        fp_normalization = normalization == :galois ? :D : normalization
        return higher_central_charge(solution, Int(n); normalization = fp_normalization)
    elseif method != :modular_data
        error("unknown higher central charge method: $method; expected :modular_data or :finite_field")
    end
    isempty(kwargs) || error("unsupported keyword arguments: $(collect(keys(kwargs)))")
    input = try
        _gauss_input(category_or_modular_data)
    catch err
        return _higher_failure(category_or_modular_data, n, normalization,
                               :missing_data, sprint(showerror, err))
    end

    try
        _check_gauss_input(input)
        τ = gauss_sum_plus(category_or_modular_data; n = n)
        D = inv(input.S[1, 1])
        D2 = total_quantum_dimension_squared(category_or_modular_data)
        denom = if normalization == :galois
            gcd(Int(n), input.context.N) == 1 ||
                return HigherCentralChargeResult(false, Int(n), normalization, nothing, τ,
                                                 nothing, D2, conductor(input.context),
                                                 :not_galois,
                                                 "n=$(Int(n)) is not a unit modulo conductor $(input.context.N)")
            galois_action(input.context, D, Int(n))
        elseif normalization == :D
            D
        elseif normalization == :D2
            D2
        elseif normalization == :raw
            one(input.context.field)
        else
            return HigherCentralChargeResult(false, Int(n), normalization, nothing, τ,
                                             nothing, D2, conductor(input.context),
                                             :unknown_normalization,
                                             "expected :galois, :D, :D2, or :raw")
        end
        value = τ / denom
        return HigherCentralChargeResult(true, Int(n), normalization, value, τ, denom,
                                         D2, conductor(input.context), :ok, "ok")
    catch err
        return HigherCentralChargeResult(false, Int(n), normalization, nothing, nothing,
                                         nothing, nothing, conductor(input.context),
                                         :error, sprint(showerror, err))
    end
end

"""
    higher_central_charges(category_or_modular_data, ns; normalization = :galois)

Compute higher central charges for each integer in `ns`.
"""
function higher_central_charges(category_or_modular_data, ns;
                                normalization::Symbol = :galois,
                                method::Symbol = :modular_data,
                                p::Union{Integer, Nothing} = nothing,
                                kwargs...)
    return [higher_central_charge(category_or_modular_data; n = n,
                                  normalization = normalization,
                                  method = method,
                                  p = p,
                                  kwargs...) for n in ns]
end

higher_central_charge(category_or_modular_data, n::Integer; kwargs...) =
    higher_central_charge(category_or_modular_data; n = n, kwargs...)

"""
    central_charge(category_or_modular_data)

Return the ordinary exact topological central charge as the `n = 1`
Galois-normalized higher central charge value.
"""
function central_charge(category_or_modular_data)
    result = higher_central_charge(category_or_modular_data; n = 1,
                                   normalization = :galois)
    result.ok || error(result.message)
    return result.value
end
