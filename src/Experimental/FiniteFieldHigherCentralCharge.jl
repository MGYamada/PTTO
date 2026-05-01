"""
Finite-field prototype for higher central charges from F/R solutions.

This is deliberately small: for the built-in examples it constructs a
finite-field realization by reducing exact cyclotomic modular data, and it
can also keep a reduced reference F/R solution for examples where Phase 4 is
practical.  The higher central charge is then evaluated as the
gauge-invariant Gauss sum over F_p.
"""

struct FRSolutionModP
    category::Symbol
    p::Int
    conductor::Int
    labels::Vector{Symbol}
    S::Matrix{Int}
    T::Vector{Int}
    F::Vector{Int}
    R::Vector{Int}
    fusion::Union{Array{Int, 3}, Nothing}
    source::Symbol
end

struct HigherCentralChargeModPResult
    ok::Bool
    n::Int
    value::Int
    gauss_sum::Int
    denominator::Int
    D_squared::Int
    p::Int
    conductor::Int
    normalization::Symbol
    status::Symbol
    message::String
end

function Base.show(io::IO, r::HigherCentralChargeModPResult)
    if r.ok
        print(io, "HigherCentralChargeModPResult(n=$(r.n), p=$(r.p), value=$(r.value))")
    else
        print(io, "HigherCentralChargeModPResult(n=$(r.n), p=$(r.p), status=$(r.status), message=$(repr(r.message)))")
    end
end

function _category_symbol(name::Symbol)
    name == :Semion && return :semion
    name == :Fibonacci && return :fibonacci
    name == :Ising && return :ising
    (name == :ToricCode || name == :toric || name == :toric_code) && return :toric_code
    return name
end

_category_symbol(name::AbstractString) =
    _category_symbol(Symbol(lowercase(replace(name, "-" => "_", " " => "_"))))

function _category_symbol(data::ModularData)
    data.labels == [:one, :s] && return :semion
    data.labels == [:one, :tau] && return :fibonacci
    data.labels == [:one, :sigma, :psi] && return :ising
    data.labels == [:one, :e, :m, :epsilon] && return :toric_code
    return :custom
end

function _builtin_data_for_finite_field(category)
    category isa ModularData && return category
    sym = _category_symbol(category)
    return modular_data(sym)
end

function _semion_fusion_rule()
    N = zeros(Int, 2, 2, 2)
    N[1, 1, 1] = 1
    N[1, 2, 2] = 1
    N[2, 1, 2] = 1
    N[2, 2, 1] = 1
    return N
end

function _fibonacci_fusion_rule()
    N = _semion_fusion_rule()
    N[2, 2, 2] = 1
    return N
end

function _ising_fusion_rule()
    N = zeros(Int, 3, 3, 3)
    N[1, 1, 1] = 1
    N[1, 2, 2] = 1
    N[1, 3, 3] = 1
    N[2, 1, 2] = 1
    N[3, 1, 3] = 1
    N[2, 2, 1] = 1
    N[2, 2, 3] = 1
    N[2, 3, 2] = 1
    N[3, 2, 2] = 1
    N[3, 3, 1] = 1
    return N
end

function _toric_code_fusion_rule()
    N = zeros(Int, 4, 4, 4)
    for a in 0:1, b in 0:1, c in 0:1, d in 0:1
        i = 1 + a + 2b
        j = 1 + c + 2d
        k = 1 + xor(a, c) + 2 * xor(b, d)
        N[i, j, k] = 1
    end
    return N
end

function _builtin_fusion_rule(sym::Symbol)
    sym == :semion && return _semion_fusion_rule()
    sym == :fibonacci && return _fibonacci_fusion_rule()
    sym == :ising && return _ising_fusion_rule()
    sym == :toric_code && return _toric_code_fusion_rule()
    return nothing
end

function _reduce_vector_mod_p(ctx::CyclotomicContext, values::AbstractVector, p::Int)
    return [reduce_mod_p(ctx, x, p) for x in values]
end

"""
    solve_fr_mod_p(category_spec, p; compute_fr = category == :fibonacci)

Experimental API.

Construct a finite-field prototype F/R solution for a built-in category.

Inputs are a built-in category name or `ModularData` object and a prime `p`.
The output is an `FRSolutionModP` containing reduced S/T data and, when
requested and available, reduced reference F/R coordinates.

The current experimental path supports `:semion`, `:fibonacci`, `:ising`,
and `:toric_code`.  It reduces exact cyclotomic S/T data to F_p and, when
`compute_fr` is true, stores a reduced reference Phase-4 F/R solution.
For higher central charges only the gauge-invariant reduced S/T data is used.

Mathematical caveats: finite-field residues are computational evidence and
are not by themselves a characteristic-zero reconstruction or uniqueness
proof.  API inputs and outputs may change before v1.0.
"""
function solve_fr_mod_p(category_spec::Union{Symbol, AbstractString, ModularData}, p::Integer;
                        compute_fr::Union{Bool, Nothing} = nothing,
                        primes::Vector{Int} = Int[],
                        kwargs...)
    warn_experimental("solve_fr_mod_p finite-field modular-data prototype")
    isempty(kwargs) || error("unsupported keyword arguments: $(collect(keys(kwargs)))")
    data = _builtin_data_for_finite_field(category_spec)
    sym = _category_symbol(category_spec isa ModularData ? data : category_spec)
    p_int = Int(p)
    reduced = reduce_mod_p(data, p_int)
    T_vec = [reduced.T[i, i] for i in 1:length(data.labels)]
    fusion = _builtin_fusion_rule(sym)

    do_fr = compute_fr === nothing ? (sym == :fibonacci) : compute_fr
    Fp = Int[]
    Rp = Int[]
    if do_fr
        fusion === nothing && error("no built-in fusion rule is available for $sym")
        isempty(primes) && (primes = sym == :fibonacci && p_int == 41 ? [41, 61] : [p_int])
        twists = [data.T[i, i] for i in 1:length(data.labels)]
        fr = compute_FR_from_ST(fusion;
                                conductor = data.context.N,
                                primes = primes,
                                S = data.S,
                                T = twists)
        Fp = _reduce_vector_mod_p(data.context, fr.F, p_int)
        Rp = _reduce_vector_mod_p(data.context, fr.R, p_int)
    end

    return FRSolutionModP(sym, p_int, data.context.N, copy(data.labels),
                          reduced.S, T_vec, Fp, Rp, fusion,
                          do_fr ? :reduced_phase4_reference : :reduced_modular_data)
end

function solve_FR_mod_p(category_spec, p::Integer; kwargs...)
    Base.depwarn("solve_FR_mod_p is deprecated; use solve_fr_mod_p instead.",
                 :solve_FR_mod_p)
    return solve_fr_mod_p(category_spec, p; kwargs...)
end

function _fp_pow_unit(a::Integer, n::Integer, p::Int)
    a = mod(a, p)
    n == 0 && return 1
    n > 0 && return powermod(a, n, p)
    a == 0 && error("cannot raise zero to a negative power in F_$p")
    return powermod(invmod(a, p), -n, p)
end

function _fp_quantum_dimensions(solution::FRSolutionModP)
    s00_inv = invmod(solution.S[1, 1], solution.p)
    return [mod(solution.S[1, i] * s00_inv, solution.p)
            for i in 1:length(solution.labels)]
end

function _fp_gauss_sum(solution::FRSolutionModP, n::Integer)
    p = solution.p
    ds = _fp_quantum_dimensions(solution)
    total = 0
    for i in eachindex(ds)
        term = _fp_pow_unit(solution.T[i], n, p)
        term = mod(term * ds[i], p)
        term = mod(term * ds[i], p)
        total = mod(total + term, p)
    end
    return total
end

function _fp_total_quantum_dimension_squared(solution::FRSolutionModP)
    p = solution.p
    return mod(sum(mod(d * d, p) for d in _fp_quantum_dimensions(solution)), p)
end

"""
    higher_central_charge(solution::FRSolutionModP, n; normalization = :D)

Compute `sum_i d_i^2 theta_i^n / D` in F_p from a finite-field F/R
solution prototype.  The default normalization matches the usual higher
central charge formula.
"""
function higher_central_charge(solution::FRSolutionModP, n::Integer;
                               normalization::Symbol = :D)
    p = solution.p
    τ = _fp_gauss_sum(solution, n)
    D = invmod(solution.S[1, 1], p)
    D2 = _fp_total_quantum_dimension_squared(solution)
    denom = if normalization == :D
        D
    elseif normalization == :D2
        D2
    elseif normalization == :raw
        1
    else
        return HigherCentralChargeModPResult(false, Int(n), 0, τ, 0, D2,
                                             p, solution.conductor, normalization,
                                             :unknown_normalization,
                                             "expected :D, :D2, or :raw")
    end
    value = mod(τ * invmod(denom, p), p)
    return HigherCentralChargeModPResult(true, Int(n), value, τ, denom, D2,
                                         p, solution.conductor, normalization,
                                         :ok, "ok")
end

higher_central_charge(solution::FRSolutionModP; n::Integer = 1,
                      normalization::Symbol = :D) =
    higher_central_charge(solution, n; normalization = normalization)

higher_central_charges(solution::FRSolutionModP, ns;
                       normalization::Symbol = :D) =
    [higher_central_charge(solution, n; normalization = normalization) for n in ns]

function _frdata_modp_hcc_solution(data::FRData{FpElem})
    p = _frdata_prime(data)
    info = _known_fr_modp_conductor_and_labels(fusion_rule(data))
    modular = _known_fr_modp_modular_data(info)
    modular === nothing &&
        error("finite-field higher central charge from FRData requires known modular metadata")
    reduced = reduce_mod_p(modular, p)
    T_vec = [reduced.T[i, i] for i in 1:length(modular.labels)]
    return FRSolutionModP(info.name, p, info.conductor, copy(modular.labels),
                          reduced.S, T_vec,
                          [x.value for x in F_values(data)],
                          [x.value for x in R_values(data)],
                          fusion_rule(data).N,
                          :frdata_modp)
end

"""
    higher_central_charge(frdata::FRData{FpElem}, n; normalization=:D)

Compute the finite-field higher central charge associated with solved
finite-field `FRData`.  Gauge fixing happens before F/R solving; this
function consumes the solved representative directly.
"""
function higher_central_charge(data::FRData{FpElem}, n::Integer;
                               normalization::Symbol = :D,
                               method::Symbol = :frdata)
    method in (:frdata, :finite_field) ||
        error("unknown FRData higher central charge method: $method; expected :frdata or :finite_field")
    return higher_central_charge(_frdata_modp_hcc_solution(data), n;
                                 normalization = normalization)
end

higher_central_charge(data::FRData{FpElem}; n::Integer = 1,
                      normalization::Symbol = :D,
                      method::Symbol = :frdata) =
    higher_central_charge(data, n; normalization = normalization,
                          method = method)

higher_central_charges(data::FRData{FpElem}, ns;
                       normalization::Symbol = :D,
                       method::Symbol = :frdata) =
    [higher_central_charge(data, n; normalization = normalization,
                           method = method) for n in ns]

"""
    lift_higher_central_charge(value_mod_p, cyclotomic_field_info)

Experimental API.

Placeholder for cyclotomic reconstruction from modular residues.

Inputs are a finite-field residue and metadata describing the intended
cyclotomic field.  The current output is the supplied residue.

For now this returns the supplied residue unless `cyclotomic_field_info` is a
`ModularData` object, in which case callers should compare by reducing the
exact cyclotomic value with `reduce_mod_p`.

Mathematical caveats: this helper does not prove a lift or perform complete
cyclotomic reconstruction.  API inputs and outputs may change before v1.0.
"""
function lift_higher_central_charge(value_mod_p, cyclotomic_field_info)
    warn_experimental("lift_higher_central_charge")
    return value_mod_p
end
