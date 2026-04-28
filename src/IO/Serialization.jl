"""
JSON serialization and lightweight reporting for classification outputs.

The loader returns the stored JSON payload as a `Dict`.  Exact Oscar objects
are exported through stable textual forms so results can be archived and
inspected without pretending that arbitrary cyclotomic field parents can be
losslessly reconstructed from JSON alone.
"""

using JSON

const ACMG_IO_SCHEMA_VERSION = 1

function _io_format_from_path(path::AbstractString)
    ext = lowercase(splitext(path)[2])
    ext == ".json" && return :json
    isempty(ext) && return :json
    error("unsupported file extension '$ext'; supported classification format is .json")
end

function _io_require_json(format::Symbol)
    format == :json || error("unsupported format: $format; supported format is :json")
    return nothing
end

_io_scalar(x::Nothing) = nothing
_io_scalar(x::Bool) = x
_io_scalar(x::Integer) = Int(x)
_io_scalar(x::AbstractFloat) = x
_io_scalar(x::AbstractString) = String(x)
_io_scalar(x::Symbol) = String(x)
_io_scalar(x) = string(x)

function _io_vector_payload(v)
    return [_io_payload(x) for x in collect(v)]
end

function _io_array_payload(A::AbstractArray)
    nd = ndims(A)
    if nd == 1
        return [_io_payload(A[i]) for i in axes(A, 1)]
    elseif nd == 2
        return [[_io_payload(A[i, j]) for j in axes(A, 2)] for i in axes(A, 1)]
    elseif nd == 3
        return [[[_io_payload(A[i, j, k]) for k in axes(A, 3)]
                 for j in axes(A, 2)] for i in axes(A, 1)]
    end
    return string(A)
end

function _io_matrix_payload(M)
    if M isa MatElem
        return [[_io_payload(M[i, j]) for j in 1:ncols(M)] for i in 1:nrows(M)]
    end
    return _io_array_payload(M)
end

function _io_payload(x)
    x === nothing && return nothing
    x isa Bool && return x
    x isa Integer && return Int(x)
    x isa AbstractFloat && return x
    x isa AbstractString && return String(x)
    x isa Symbol && return String(x)
    x isa AbstractArray && return _io_array_payload(x)
    x isa MatElem && return _io_matrix_payload(x)
    return string(x)
end

function _io_stratum_payload(s::Stratum)
    pairs = sort(collect(s.multiplicities); by = first)
    return Dict{String, Any}(
        "multiplicities" => Dict(string(k) => v for (k, v) in pairs),
        "total_dim" => s.total_dim,
        "display" => sprint(show, s),
    )
end

function _io_report_payload(r::FRRoundtripReport)
    return Dict{String, Any}(
        "ok" => r.ok,
        "S_error" => _io_payload(r.S_error),
        "T_error" => _io_payload(r.T_error),
        "S_max" => _io_payload(r.S_max),
        "T_max" => _io_payload(r.T_max),
        "best_perm" => r.best_perm,
        "S_roundtrip" => _io_matrix_payload(r.S_roundtrip),
        "T_roundtrip" => _io_payload(r.T_roundtrip),
        "candidate_index" => r.candidate_index,
        "galois_exponent" => r.galois_exponent,
    )
end

_io_report_payload(::Nothing) = nothing

function _io_modular_data_payload(data::ModularData)
    return Dict{String, Any}(
        "type" => "ModularData",
        "N" => conductor(data.context),
        "labels" => [String(x) for x in data.labels],
        "S" => _io_matrix_payload(data.S),
        "T" => _io_matrix_payload(data.T),
        "cond_S" => data.cond_S,
        "cond_T" => data.cond_T,
        "cond_F" => data.cond_F,
    )
end

function _io_classified_payload(m::ClassifiedMTC)
    return Dict{String, Any}(
        "type" => "ClassifiedMTC",
        "N" => m.N,
        "N_input" => m.N_input,
        "rank" => m.rank,
        "stratum" => _io_stratum_payload(m.stratum),
        "Nijk" => _io_array_payload(m.Nijk),
        "fusion_rule_key" => fusion_rule_key(m.Nijk),
        "scale_factor" => m.scale_factor,
        "used_primes" => m.used_primes,
        "fresh_primes" => m.fresh_primes,
        "verify_fresh" => m.verify_fresh,
        "verify_exact_lift" => m.verify_exact_lift,
        "S_cyclotomic" => _io_matrix_payload(m.S_cyclotomic),
        "T_cyclotomic" => _io_payload(m.T_cyclotomic),
        "F_values" => _io_payload(m.F_values),
        "R_values" => _io_payload(m.R_values),
        "galois_sector" => m.galois_sector,
        "fr_status" => string(m.fr_status),
        "verify_report" => _io_report_payload(m.verify_report),
    )
end

function _io_classification_payload(result)
    items = result isa AbstractVector ? result : [result]
    return Dict{String, Any}(
        "schema" => "ACMG.classification",
        "schema_version" => ACMG_IO_SCHEMA_VERSION,
        "count" => length(items),
        "items" => [_io_classified_payload(m) for m in items],
    )
end

function _io_json_string(payload)
    return JSON.json(payload)
end

"""
    save_classification(path, result; format = nothing)

Write one `ClassifiedMTC` or a vector of them to JSON.  Exact cyclotomic
values are stored as textual coordinates suitable for archiving and review.
"""
function save_classification(path::AbstractString, result; format = nothing)
    fmt = format === nothing ? _io_format_from_path(path) : Symbol(format)
    _io_require_json(fmt)
    payload = _io_classification_payload(result)
    open(path, "w") do io
        JSON.print(io, payload, 2)
        write(io, '\n')
    end
    return path
end

"""
    load_classification(path; format = nothing)

Read a classification JSON file written by `save_classification`.
Returns the JSON payload as a `Dict`.
"""
function load_classification(path::AbstractString; format = nothing)
    fmt = format === nothing ? _io_format_from_path(path) : Symbol(format)
    _io_require_json(fmt)
    return JSON.parsefile(path)
end

function _io_export(payload; format::Symbol)
    format == :dict && return payload
    _io_require_json(format)
    return _io_json_string(payload)
end

"""
    export_modular_data(result; format = :json)

Export exact S/T data from `ClassifiedMTC` or `ModularData`.
Use `format = :dict` for the raw Julia payload.
"""
function export_modular_data(result; format = :json)
    fmt = Symbol(format)
    payload = if result isa ModularData
        _io_modular_data_payload(result)
    elseif result isa ClassifiedMTC
        Dict{String, Any}(
            "type" => "ModularData",
            "N" => result.N,
            "rank" => result.rank,
            "S" => _io_matrix_payload(result.S_cyclotomic),
            "T" => _io_payload(result.T_cyclotomic),
            "galois_sector" => result.galois_sector,
        )
    else
        error("export_modular_data expects ModularData or ClassifiedMTC")
    end
    return _io_export(payload; format = fmt)
end

"""
    export_fusion_rule(result; format = :json)

Export the fusion tensor `Nijk` from `ClassifiedMTC`, `FusionRule`, or a
rank-3 integer array.
"""
function export_fusion_rule(result; format = :json)
    fmt = Symbol(format)
    Nijk = result isa ClassifiedMTC ? result.Nijk :
           result isa FusionRule ? result.N :
           result isa AbstractArray ? result :
           error("export_fusion_rule expects ClassifiedMTC, FusionRule, or an array")
    payload = Dict{String, Any}(
        "type" => "FusionRule",
        "rank" => size(Nijk, 1),
        "Nijk" => _io_array_payload(Nijk),
        "fusion_rule_key" => canonical_rule(Nijk),
    )
    return _io_export(payload; format = fmt)
end

"""
    export_FR(result; format = :json)

Export the exact F/R layer from a `ClassifiedMTC`.
"""
function export_FR(result::ClassifiedMTC; format = :json)
    fmt = Symbol(format)
    payload = Dict{String, Any}(
        "type" => "FRData",
        "N" => result.N,
        "rank" => result.rank,
        "F_values" => _io_payload(result.F_values),
        "R_values" => _io_payload(result.R_values),
        "verify_report" => _io_report_payload(result.verify_report),
        "galois_sector" => result.galois_sector,
        "fr_status" => string(result.fr_status),
    )
    return _io_export(payload; format = fmt)
end

function _io_report_lines(result::ClassifiedMTC)
    fr = result.verify_report
    lines = String[
        "## ClassifiedMTC",
        "",
        "- N: $(result.N)",
        "- input N: $(result.N_input)",
        "- rank: $(result.rank)",
        "- galois sector: $(result.galois_sector)",
        "- stratum: $(sprint(show, result.stratum))",
        "- used primes: $(join(result.used_primes, ", "))",
        "- fresh primes: $(isempty(result.fresh_primes) ? "none" : join(result.fresh_primes, ", "))",
        "- fresh verification: $(result.verify_fresh)",
        "- exact lift verification: $(result.verify_exact_lift)",
        "- F/R status: $(result.fr_status)",
        "- F/R attached: $(result.F_values !== nothing && result.R_values !== nothing)",
    ]
    if fr !== nothing
        push!(lines, "- F/R roundtrip ok: $(fr.ok)")
        push!(lines, "- S error: $(fr.S_error)")
        push!(lines, "- T error: $(fr.T_error)")
        push!(lines, "- best permutation: $(fr.best_perm)")
        push!(lines, "- candidate index: $(fr.candidate_index)")
        push!(lines, "- galois exponent: $(fr.galois_exponent)")
    end
    push!(lines, "")
    return lines
end

"""
    write_report(path, result)

Write a concise Markdown report for one classification result or a vector of
results.
"""
function write_report(path::AbstractString, result)
    items = result isa AbstractVector ? result : [result]
    lines = String["# ACMG Classification Report", "", "Results: $(length(items))", ""]
    for (i, item) in enumerate(items)
        length(items) > 1 && push!(lines, "### Result $i", "")
        append!(lines, _io_report_lines(item))
    end
    open(path, "w") do io
        write(io, join(lines, '\n'))
    end
    return path
end
