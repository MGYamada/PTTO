"""
Backend-neutral F/R equation infrastructure for multiplicity-free fusion rules.

This layer is deliberately lightweight: it records variables and sparse
polynomial equations without committing to Oscar, Symbolics, Nemo, or any
particular solver backend.  The first supported scope is multiplicity-free
fusion rules with tensor unit at object `1`.
"""

const FRIndex = Tuple{Vararg{Int}}

struct EquationVariable
    name::Symbol
    kind::Symbol
    indices::FRIndex
    domain::Symbol
end

struct EquationTerm
    coeff::Any
    powers::Vector{Pair{EquationVariable, Int}}
end

struct EquationExpr
    terms::Vector{EquationTerm}
end

struct PolynomialEquation
    lhs::EquationExpr
    rhs::EquationExpr
    metadata::Dict{Symbol, Any}
end

struct EquationSystem
    variables::Vector{EquationVariable}
    equations::Vector{Any}
    metadata::Dict{Symbol, Any}
    assumptions::Vector{String}
    base_ring::Symbol
end

struct GaugeVariable
    variable::EquationVariable
    a::Int
    b::Int
    c::Int
end

struct FREquationSystem
    rules::FusionRule
    variables::Vector{EquationVariable}
    equations::Vector{Any}
    metadata::Dict{Symbol, Any}
    assumptions::Vector{String}
    base_ring::Symbol
end

struct FiniteFieldEquationSystem
    system::FREquationSystem
    p::Int
    variables::Vector{EquationVariable}
    equations::Vector{Any}
    metadata::Dict{Symbol, Any}
end

_fusion_rule(fr::FusionRule) = fr
_fusion_rule(Nijk::Array{Int, 3}) = FusionRule(Nijk)
_fusion_tensor(fr::FusionRule) = fr.N
_fusion_tensor(Nijk::Array{Int, 3}) = Nijk

_fr_name(prefix::AbstractString, xs::Int...) = Symbol(prefix * "_" * join(xs, "_"))

_zero_expr() = EquationExpr(EquationTerm[])
_one_expr() = EquationExpr([EquationTerm(1, Pair{EquationVariable, Int}[])])
_const_expr(c) = iszero(c) ? _zero_expr() : EquationExpr([EquationTerm(c, Pair{EquationVariable, Int}[])])
_var_expr(v::EquationVariable) = EquationExpr([EquationTerm(1, [v => 1])])

function _normalize_powers(powers::Vector{Pair{EquationVariable, Int}})
    acc = Dict{EquationVariable, Int}()
    for pair in powers
        pair.second == 0 && continue
        acc[pair.first] = get(acc, pair.first, 0) + pair.second
    end
    out = [v => e for (v, e) in acc if e != 0]
    sort!(out; by = p -> string(p.first.name))
    return out
end

function _combine_terms(terms::Vector{EquationTerm})
    acc = Dict{Tuple{Vararg{Pair{EquationVariable, Int}}}, Any}()
    key_order = Vector{Tuple{Vararg{Pair{EquationVariable, Int}}}}()
    for t in terms
        iszero(t.coeff) && continue
        powers = _normalize_powers(t.powers)
        key = Tuple(powers)
        if !haskey(acc, key)
            push!(key_order, key)
            acc[key] = t.coeff
        else
            acc[key] += t.coeff
        end
    end
    out = EquationTerm[]
    for key in key_order
        coeff = acc[key]
        iszero(coeff) || push!(out, EquationTerm(coeff, collect(key)))
    end
    return out
end

Base.:+(a::EquationExpr, b::EquationExpr) = EquationExpr(_combine_terms(vcat(a.terms, b.terms)))
Base.:-(a::EquationExpr) = EquationExpr([EquationTerm(-t.coeff, copy(t.powers)) for t in a.terms])
Base.:-(a::EquationExpr, b::EquationExpr) = a + (-b)

function Base.:*(a::EquationExpr, b::EquationExpr)
    terms = EquationTerm[]
    for ta in a.terms, tb in b.terms
        push!(terms, EquationTerm(ta.coeff * tb.coeff, vcat(ta.powers, tb.powers)))
    end
    return EquationExpr(_combine_terms(terms))
end

Base.iszero(e::EquationExpr) = isempty(e.terms)

_equation(lhs::EquationExpr, rhs::EquationExpr; metadata = Dict{Symbol, Any}()) =
    PolynomialEquation(lhs, rhs, metadata)

function Base.show(io::IO, v::EquationVariable)
    print(io, v.name)
end

function Base.show(io::IO, t::EquationTerm)
    if isempty(t.powers)
        print(io, t.coeff)
        return
    end
    if t.coeff != 1
        print(io, t.coeff, "*")
    end
    print(io, join([string(p.first.name) * (p.second == 1 ? "" : "^$(p.second)")
                    for p in t.powers], "*"))
end

function Base.show(io::IO, e::EquationExpr)
    isempty(e.terms) && (print(io, "0"); return)
    print(io, join(string.(e.terms), " + "))
end

function Base.show(io::IO, eq::PolynomialEquation)
    print(io, eq.lhs, " = ", eq.rhs)
end

function Base.show(io::IO, ::MIME"text/plain", sys::EquationSystem)
    print(io, "EquationSystem($(length(sys.variables)) variables, ",
          "$(length(sys.equations)) equations, base_ring = :$(sys.base_ring))")
end

function Base.show(io::IO, ::MIME"text/plain", sys::FREquationSystem)
    print(io, "FREquationSystem(rank = $(sys.rules.rank), variables = $(length(sys.variables)), ",
          "equations = $(length(sys.equations)))")
end

function Base.show(io::IO, ::MIME"text/plain", sys::FiniteFieldEquationSystem)
    print(io, "FiniteFieldEquationSystem(p = $(sys.p), variables = ",
          "$(length(sys.variables)), equations = $(length(sys.equations)))")
end

"""
    simple_objects(rules)

Return the deterministic object labels `1:rank` for a fusion rule.  The unit
object is assumed to be label `1`, matching the rest of ACMG.
"""
simple_objects(rules) = collect(1:_fusion_rule(rules).rank)

"""
    fusion_product(rules, a, b)

Return pairs `(c, multiplicity)` occurring in `a ⊗ b`.
"""
function fusion_product(rules, a::Int, b::Int)
    fr = _fusion_rule(rules)
    _check_object(fr, a); _check_object(fr, b)
    return [(c, fr.N[a, b, c]) for c in 1:fr.rank if fr.N[a, b, c] != 0]
end

fusion_channels(rules, a::Int, b::Int) = [c for (c, _) in fusion_product(rules, a, b)]

function is_admissible(rules, a::Int, b::Int, c::Int)
    fr = _fusion_rule(rules)
    _check_object(fr, a); _check_object(fr, b); _check_object(fr, c)
    return fr.N[a, b, c] != 0
end

function admissible_triples(rules)
    fr = _fusion_rule(rules)
    return [(a, b, c) for a in 1:fr.rank, b in 1:fr.rank, c in 1:fr.rank
            if fr.N[a, b, c] != 0]
end

function admissible_quadruples(rules)
    fr = _fusion_rule(rules)
    return [(a, b, c, d) for a in 1:fr.rank, b in 1:fr.rank,
            c in 1:fr.rank, d in 1:fr.rank
            if any(e -> fr.N[a, b, e] != 0 && fr.N[e, c, d] != 0, 1:fr.rank)]
end

is_multiplicity_free(rules) = all(x -> x == 0 || x == 1, _fusion_tensor(rules))

function require_multiplicity_free(rules)
    is_multiplicity_free(rules) ||
        error("v0.8 F/R equation infrastructure currently supports multiplicity-free fusion rules only")
    return true
end

function _check_object(fr::FusionRule, a::Int)
    1 <= a <= fr.rank || error("object label $a is outside 1:$(fr.rank)")
    return true
end

_vexpr(x::GaugeVariable) = _var_expr(x.variable)

function _forward_r_var_count_for_fusion(rules)
    fr = _fusion_rule(rules)
    return sum(fr.N[a, b, c]^2 for a in 1:fr.rank, b in 1:fr.rank, c in 1:fr.rank)
end

"""
    pentagon_equations(rules)

Return the TensorCategories-backed pentagon equations used by
`get_pentagon_system`.  Variable ordering is the TensorCategories/internal
`PentagonEquations.jl` ordering.
"""
function pentagon_equations(rules::Union{FusionRule, Array{Int, 3}};
                            normalize_unit::Bool = true)
    normalize_unit || @warn "normalize_unit=false is ignored by TensorCategories-backed pentagon_equations"
    fr = _fusion_rule(rules)
    require_multiplicity_free(fr)
    _, eqs, _ = get_pentagon_system(fr.N, fr.rank)
    return eqs
end

function _unique_equations(eqs::AbstractVector)
    seen = Set{String}()
    out = Any[]
    for eq in eqs
        key = eq isa PolynomialEquation ? string(eq.lhs - eq.rhs) : string(eq)
        key in seen && continue
        push!(seen, key)
        push!(out, eq)
    end
    return out
end

function hexagon_equations(rules::Union{FusionRule, Array{Int, 3}}; context = nothing)
    fr = _fusion_rule(rules)
    require_multiplicity_free(fr)
    _, eqs, _ = get_hexagon_fr_system(fr.N, fr.rank; context = context)
    return eqs
end

"""
    fr_equation_system(rules; include_pentagon=true, include_hexagon=true)

Build a backend-neutral F/R equation system for a multiplicity-free fusion
rule.
"""
function fr_equation_system(rules; include_pentagon::Bool = true,
                            include_hexagon::Bool = true,
                            normalize_unit::Bool = true,
                            base_ring::Symbol = :ZZ)
    fr = _fusion_rule(rules)
    require_multiplicity_free(fr)
    eqs = Any[]
    Rpent, _, nF = get_pentagon_system(fr.N, fr.rank)
    nR = _forward_r_var_count_for_fusion(fr)
    include_pentagon && append!(eqs, pentagon_equations(fr;
                                                        normalize_unit = normalize_unit))
    if include_hexagon
        append!(eqs, hexagon_equations(fr))
    end
    variables = EquationVariable[]
    assumptions = ["multiplicity-free fusion rule", "unit object has label 1",
                   "pentagon and hexagon equations are TensorCategories-backed"]
    metadata = Dict{Symbol, Any}(:rank => fr.rank,
                                 :f_variables => nF,
                                 :r_variables => nR,
                                 :hexagon_inverse_r_variables => include_hexagon ? nR : 0,
                                 :include_pentagon => include_pentagon,
                                 :include_hexagon => include_hexagon,
                                 :hexagon_variables_include_F_and_R => include_hexagon,
                                 :normalizations => normalize_unit)
    return FREquationSystem(fr, variables, _unique_equations(eqs),
                            metadata, assumptions, base_ring)
end

function validate_fr_system(system::FREquationSystem)
    require_multiplicity_free(system.rules)
    vars = Set(system.variables)
    for eq in system.equations
        eq isa PolynomialEquation || continue
        for side in (eq.lhs, eq.rhs), term in side.terms, p in term.powers
        p.first in vars || error("equation references variable $(p.first.name) not present in system")
        p.second > 0 || error("negative or zero exponents are not supported in polynomial equations")
        end
    end
    return true
end

validate(system::FREquationSystem) = validate_fr_system(system)

function gauge_variables(rules)
    fr = _fusion_rule(rules)
    require_multiplicity_free(fr)
    out = GaugeVariable[]
    for (a, b, c) in admissible_triples(fr)
        v = EquationVariable(_fr_name("u", a, b, c), :gauge, (a, b, c), :nonzero)
        push!(out, GaugeVariable(v, a, b, c))
    end
    return out
end

function _gauge_map(gvars::Vector{GaugeVariable})
    return Dict((v.a, v.b, v.c) => v for v in gvars)
end

function gauge_transform(system::FREquationSystem)
    gvars = gauge_variables(system.rules)
    return (gauge_variables = gvars, F = NamedTuple[], R = NamedTuple[],
            note = "symbolic F/R variable coordinate actions were removed with fsymbol_variables/rsymbol_variables")
end

function gauge_fix(system::FREquationSystem; strategy::Symbol = :safe)
    strategy == :safe || error("only strategy=:safe is implemented for v0.8 gauge_fix")
    fixed_eqs = copy(system.equations)
    fixed = Vector{NamedTuple}()
    meta = copy(system.metadata)
    meta[:gauge_fix_strategy] = :safe
    meta[:fixed_symbols] = fixed
    meta[:gauge_fix_note] = "TensorCategories-backed equations already carry their variable convention; obsolete F/R coordinate fixes are not added"
    meta[:residual_gauge_variables] = [v.variable.name for v in gauge_variables(system.rules)
                                       if v.a != 1 && v.b != 1]
    return FREquationSystem(system.rules, system.variables,
                            _unique_equations(fixed_eqs), meta, system.assumptions,
                            system.base_ring)
end

function _mod_coeff(c, p::Int)
    c isa Integer && return mod(c, p)
    if c isa Rational
        den = mod(denominator(c), p)
        den == 0 && error("bad prime $p: denominator $(denominator(c)) vanishes modulo p")
        return mod(numerator(c), p) * invmod(den, p) % p
    end
    return c
end

function _reduce_expr_mod_p(e::EquationExpr, p::Int)
    return EquationExpr(_combine_terms([EquationTerm(_mod_coeff(t.coeff, p), copy(t.powers))
                                        for t in e.terms]))
end

function reduce_mod_p(system::FREquationSystem, p::Integer)
    p > 1 && isprime(Int(p)) || error("p must be prime, got $p")
    pp = Int(p)
    eqs = Any[eq isa PolynomialEquation ?
              PolynomialEquation(_reduce_expr_mod_p(eq.lhs, pp),
                                 _reduce_expr_mod_p(eq.rhs, pp),
                                 copy(eq.metadata)) : eq
              for eq in system.equations]
    meta = copy(system.metadata)
    meta[:base_field] = Symbol("F_$pp")
    meta[:tensorcategories_equations_copied] = any(eq -> !(eq isa PolynomialEquation), system.equations)
    if haskey(meta, :conductor)
        meta[:p_mod_conductor] = mod(pp, meta[:conductor])
    end
    return FiniteFieldEquationSystem(system, pp, copy(system.variables), eqs, meta)
end

"""
    solve_finite_field(system::FiniteFieldEquationSystem; kwargs...)

Experimental API.

Attempt to solve a reduced F/R equation system over a finite field.

The input is a `FiniteFieldEquationSystem`, normally produced by
`reduce_mod_p(::FREquationSystem, p)`.  A future implementation may return
finite-field solution records; v0.8.6 only validates the interface and reports
that the general solver is not implemented.

Mathematical caveats: finite-field solutions require lifting and exact
verification before they should be treated as cyclotomic F/R data.  API inputs
and outputs may change before v1.0.
"""
function solve_finite_field(system::FiniteFieldEquationSystem; kwargs...)
    warn_experimental("solve_finite_field")
    error("solve_finite_field is not implemented for general FR systems; v0.8 only provides reduction infrastructure")
end

"""
    cyclotomic_reconstruct(fp_solution; conductor::Integer)

Experimental API.

Reconstruct cyclotomic data from finite-field solution data.

The input is a finite-field solution-like object and a positive conductor.
The intended output is a cyclotomic lift, but v0.8.6 only validates the
conductor and reports that the reconstruction backend is incomplete.

Mathematical caveats: modular residues are not a proof of a cyclotomic lift
without CRT bounds and exact equation verification.  API inputs and outputs
may change before v1.0.
"""
function cyclotomic_reconstruct(fp_solution; conductor::Integer)
    warn_experimental("cyclotomic_reconstruct")
    conductor > 0 || error("conductor must be positive")
    error("cyclotomic_reconstruct is experimental and not implemented beyond interface validation")
end

function frobenius_metadata(p::Integer, conductor::Integer)
    p > 1 && isprime(Int(p)) || error("p must be prime, got $p")
    conductor > 0 || error("conductor must be positive")
    return Dict(:p => Int(p), :conductor => Int(conductor),
                :p_mod_conductor => mod(Int(p), Int(conductor)),
                :split_prime_hint => mod(Int(p), Int(conductor)) == 1)
end

function check_modular_data(candidate, known_data)
    return candidate == known_data
end

function semion_fusion_rules()
    N = zeros(Int, 2, 2, 2)
    N[1, 1, 1] = 1
    N[1, 2, 2] = 1
    N[2, 1, 2] = 1
    N[2, 2, 1] = 1
    return FusionRule(N)
end

function fibonacci_fusion_rules()
    N = zeros(Int, 2, 2, 2)
    N[1, 1, 1] = 1
    N[1, 2, 2] = 1
    N[2, 1, 2] = 1
    N[2, 2, 1] = 1
    N[2, 2, 2] = 1
    return FusionRule(N)
end

function toric_code_fusion_rules()
    N = zeros(Int, 4, 4, 4)
    for a in 1:4, b in 1:4
        c = xor(a - 1, b - 1) + 1
        N[a, b, c] = 1
    end
    return FusionRule(N)
end

function ising_fusion_rules()
    N = zeros(Int, 3, 3, 3)
    N[1, 1, 1] = 1
    N[1, 2, 2] = 1; N[2, 1, 2] = 1
    N[1, 3, 3] = 1; N[3, 1, 3] = 1
    N[2, 2, 1] = 1
    N[2, 2, 3] = 1
    N[2, 3, 2] = 1; N[3, 2, 2] = 1
    N[3, 3, 1] = 1
    return FusionRule(N)
end
