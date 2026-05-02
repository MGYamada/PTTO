"""
Finite-field higher central charges for solved `FRData{FpElem}`.

The finite-field F/R solver itself lives in `FR/FiniteFieldFRData.jl`.
This file deliberately avoids defining public solver entry points or public
result types; it only adapts solved finite-field F/R data to the HCC helpers.
"""

function _fp_hcc_result(; ok::Bool, n::Integer, value::Integer, gauss_sum::Integer,
                        denominator::Integer, D_squared::Integer, p::Integer,
                        conductor_value, normalization::Symbol, status::Symbol,
                        message::String)
    return (ok = ok,
            n = Int(n),
            value = Int(value),
            gauss_sum = Int(gauss_sum),
            denominator = Int(denominator),
            D_squared = Int(D_squared),
            p = Int(p),
            conductor = conductor_value,
            normalization = normalization,
            status = status,
            message = message)
end

function _fp_hcc_pow_unit(a::Integer, n::Integer, p::Int)
    a = mod(a, p)
    n == 0 && return 1
    n > 0 && return powermod(a, n, p)
    a == 0 && error("cannot raise zero to a negative power in F_$p")
    return powermod(invmod(a, p), -n, p)
end

function _fp_hcc_quantum_dimensions(data::FRData{FpElem}, p::Int)
    info = _known_fr_modp_conductor_and_labels(fusion_rule(data))
    dims = _known_quantum_dimensions_mod_p(info, p)
    dims === nothing &&
        error("finite-field higher central charge requires known built-in modular metadata")
    return dims
end

function _fp_hcc_twists(data::FRData{FpElem}, p::Int)
    info = _known_fr_modp_conductor_and_labels(fusion_rule(data))
    twists = _known_target_twists_mod_p(info, p)
    twists === nothing &&
        error("finite-field higher central charge requires known built-in twist metadata")
    return twists
end

function _fp_hcc_conductor(data::FRData{FpElem})
    meta = fr_metadata(data)
    haskey(meta, :conductor) && return Int(meta[:conductor])
    info = _known_fr_modp_conductor_and_labels(fusion_rule(data))
    return info.conductor
end

function _fp_hcc_gauss_sum(data::FRData{FpElem}, n::Integer, p::Int)
    ds = _fp_hcc_quantum_dimensions(data, p)
    twists = _fp_hcc_twists(data, p)
    length(ds) == length(twists) ||
        error("dimension/twist metadata length mismatch")
    total = 0
    for i in eachindex(ds)
        term = _fp_hcc_pow_unit(twists[i], n, p)
        term = _mul_mod_int(term, ds[i], p)
        term = _mul_mod_int(term, ds[i], p)
        total = mod(total + term, p)
    end
    return total
end

function _fp_hcc_total_dimension_squared(data::FRData{FpElem}, p::Int)
    return mod(sum(mod(d * d, p) for d in _fp_hcc_quantum_dimensions(data, p)), p)
end

"""
    higher_central_charge(frdata::FRData{FpElem}, n; normalization = :D2)

Compute the finite-field higher central charge associated with solved
finite-field `FRData`.  This compatibility method is intentionally thin:
finite-field F/R solving is handled by `solve_fr_mod_p` in the FR layer, and
exact modular-data reductions should use `higher_central_charge_modp`.
"""
function higher_central_charge(data::FRData{FpElem}, n::Integer;
                               normalization::Symbol = :D2,
                               method::Symbol = :frdata)
    method in (:frdata, :finite_field) ||
        error("unknown FRData higher central charge method: $method; expected :frdata or :finite_field")
    p = _frdata_prime(data)
    τ = _fp_hcc_gauss_sum(data, n, p)
    D2 = _fp_hcc_total_dimension_squared(data, p)
    denom = if normalization == :D2
        D2
    elseif normalization == :raw
        1
    else
        return _fp_hcc_result(ok = false, n = n, value = 0, gauss_sum = τ,
                              denominator = 0, D_squared = D2, p = p,
                              conductor_value = _fp_hcc_conductor(data),
                              normalization = normalization,
                              status = :unknown_normalization,
                              message = "expected :D2 or :raw")
    end
    denom != 0 || error("bad prime $p for HCC reduction: denominator is zero in F_$p")
    value = mod(τ * invmod(denom, p), p)
    return _fp_hcc_result(ok = true, n = n, value = value, gauss_sum = τ,
                          denominator = denom, D_squared = D2, p = p,
                          conductor_value = _fp_hcc_conductor(data),
                          normalization = normalization,
                          status = :ok, message = "ok")
end

higher_central_charge(data::FRData{FpElem}; n::Integer = 1,
                      normalization::Symbol = :D2,
                      method::Symbol = :frdata) =
    higher_central_charge(data, n; normalization = normalization,
                          method = method)

higher_central_charges(data::FRData{FpElem}, ns;
                       normalization::Symbol = :D2,
                       method::Symbol = :frdata) =
    [higher_central_charge(data, n; normalization = normalization,
                           method = method) for n in ns]
