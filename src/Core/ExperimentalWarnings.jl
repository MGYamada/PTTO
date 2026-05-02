const _EXPERIMENTAL_WARNED = Set{String}()

_experimental_warnings_enabled() =
    get(ENV, "ACMG_WARN_EXPERIMENTAL", "false") == "true"

"""
    warn_experimental(name::AbstractString)

Emit an opt-in warning for an experimental ACMG API.

Warnings are enabled only when `ENV["ACMG_WARN_EXPERIMENTAL"] == "true"`.
Each API name is warned at most once per Julia session so experimental paths
can call this helper without flooding normal computations.
"""
function warn_experimental(name::AbstractString)
    _experimental_warnings_enabled() || return false
    name in _EXPERIMENTAL_WARNED && return false
    push!(_EXPERIMENTAL_WARNED, String(name))
    @warn "$name is an experimental ACMG API; inputs, outputs, and semantics may change. Pin ACMG.jl when depending on this interface." maxlog=1
    return true
end
