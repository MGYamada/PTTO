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
    equations::Vector{PolynomialEquation}
    metadata::Dict{Symbol, Any}
    assumptions::Vector{String}
    base_ring::Symbol
end

struct FSymbolVariable
    variable::EquationVariable
    a::Int
    b::Int
    c::Int
    d::Int
    e::Int
    f::Int
end

struct RSymbolVariable
    variable::EquationVariable
    a::Int
    b::Int
    c::Int
end

struct GaugeVariable
    variable::EquationVariable
    a::Int
    b::Int
    c::Int
end

struct FREquationSystem
    rules::FusionRule
    fvars::Vector{FSymbolVariable}
    rvars::Vector{RSymbolVariable}
    variables::Vector{EquationVariable}
    equations::Vector{PolynomialEquation}
    metadata::Dict{Symbol, Any}
    assumptions::Vector{String}
    base_ring::Symbol
end

struct FiniteFieldEquationSystem
    system::FREquationSystem
    p::Int
    variables::Vector{EquationVariable}
    equations::Vector{PolynomialEquation}
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
    print(io, "FREquationSystem(rank = $(sys.rules.rank), F = $(length(sys.fvars)), ",
          "R = $(length(sys.rvars)), equations = $(length(sys.equations)))")
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

"""
    fsymbol_variables(rules)

Generate valid multiplicity-free F-symbol variables
`F^{a,b,c}_d[e,f]`.
"""
function fsymbol_variables(rules)
    fr = _fusion_rule(rules)
    require_multiplicity_free(fr)
    out = FSymbolVariable[]
    for a in 1:fr.rank, b in 1:fr.rank, c in 1:fr.rank, d in 1:fr.rank
        left = [e for e in 1:fr.rank if fr.N[a, b, e] != 0 && fr.N[e, c, d] != 0]
        right = [f for f in 1:fr.rank if fr.N[b, c, f] != 0 && fr.N[a, f, d] != 0]
        for e in left, f in right
            v = EquationVariable(_fr_name("F", a, b, c, d, e, f), :F,
                                 (a, b, c, d, e, f), :nonzero)
            push!(out, FSymbolVariable(v, a, b, c, d, e, f))
        end
    end
    return out
end

"""
    rsymbol_variables(rules)

Generate valid multiplicity-free R-symbol variables `R^{a,b}_c`.
"""
function rsymbol_variables(rules)
    fr = _fusion_rule(rules)
    require_multiplicity_free(fr)
    out = RSymbolVariable[]
    for a in 1:fr.rank, b in 1:fr.rank, c in 1:fr.rank
        fr.N[a, b, c] == 0 && continue
        v = EquationVariable(_fr_name("R", a, b, c), :R, (a, b, c), :nonzero)
        push!(out, RSymbolVariable(v, a, b, c))
    end
    return out
end

function _fmap(fvars::Vector{FSymbolVariable})
    return Dict((v.a, v.b, v.c, v.d, v.e, v.f) => v for v in fvars)
end

function _rmap(rvars::Vector{RSymbolVariable})
    return Dict((v.a, v.b, v.c) => v for v in rvars)
end

_fv(fdict, a, b, c, d, e, f) = get(fdict, (a, b, c, d, e, f), nothing)
_rv(rdict, a, b, c) = get(rdict, (a, b, c), nothing)
_vexpr(x::FSymbolVariable) = _var_expr(x.variable)
_vexpr(x::RSymbolVariable) = _var_expr(x.variable)
_vexpr(x::GaugeVariable) = _var_expr(x.variable)

function unit_fsymbol_normalizations(rules, fvars = fsymbol_variables(rules))
    return [v for v in fvars if v.a == 1 || v.b == 1 || v.c == 1]
end

function unit_rsymbol_normalizations(rules, rvars = rsymbol_variables(rules))
    return [v for v in rvars if v.a == 1 || v.b == 1]
end

"""
    pentagon_equations(rules, fvars)

Generate multiplicity-free pentagon equations in a sparse backend-neutral
representation.  The convention records both paths from `(((a*b)*c)*d)` to
`a*(b*(c*d))`; impossible channels are skipped.
"""
function pentagon_equations(rules, fvars::Vector{FSymbolVariable};
                            normalize_unit::Bool = true)
    fr = _fusion_rule(rules)
    require_multiplicity_free(fr)
    fdict = _fmap(fvars)
    eqs = PolynomialEquation[]
    for a in 1:fr.rank, b in 1:fr.rank, c in 1:fr.rank, d in 1:fr.rank
        for x in fusion_channels(fr, a, b),
            y in fusion_channels(fr, x, c),
            z in fusion_channels(fr, y, d),
            u in fusion_channels(fr, b, c),
            v in fusion_channels(fr, a, u)

            fr.N[v, d, z] == 0 && continue
            lhs1 = _fv(fdict, a, b, c, v, x, u)
            lhs1 === nothing && continue
            for w in fusion_channels(fr, u, d)
                fr.N[a, w, z] == 0 && continue
                lhs2 = _fv(fdict, a, u, d, z, v, w)
                lhs2 === nothing && continue
                lhs = _vexpr(lhs1) * _vexpr(lhs2)
                rhs = _zero_expr()
                for n in fusion_channels(fr, c, d)
                    fr.N[x, n, z] == 0 && continue
                    fr.N[b, n, w] == 0 && continue
                    r1 = _fv(fdict, x, c, d, z, y, n)
                    r2 = _fv(fdict, a, b, n, z, x, w)
                    r3 = _fv(fdict, b, c, d, w, u, n)
                    (r1 === nothing || r2 === nothing || r3 === nothing) && continue
                    rhs = rhs + (_vexpr(r1) * _vexpr(r2) * _vexpr(r3))
                end
                iszero(lhs - rhs) || push!(eqs, _equation(lhs, rhs;
                    metadata = Dict(:kind => :pentagon, :objects => (a, b, c, d),
                                    :target => z, :channels => (x, y, u, v, w))))
            end
        end
    end
    if normalize_unit
        for v in unit_fsymbol_normalizations(fr, fvars)
            push!(eqs, _equation(_vexpr(v), _one_expr();
                metadata = Dict(:kind => :unit_F_normalization, :indices => v.variable.indices)))
        end
    end
    return _unique_equations(eqs)
end

pentagon_equations(rules::Union{FusionRule, Array{Int, 3}};
                   normalize_unit::Bool = true) =
    pentagon_equations(rules, fsymbol_variables(rules);
                       normalize_unit = normalize_unit)

function _unique_equations(eqs::Vector{PolynomialEquation})
    seen = Set{String}()
    out = PolynomialEquation[]
    for eq in eqs
        key = string(eq.lhs - eq.rhs)
        key in seen && continue
        push!(seen, key)
        push!(out, eq)
    end
    return out
end

function left_hexagon_equations(rules, fvars = fsymbol_variables(rules),
                                rvars = rsymbol_variables(rules);
                                normalize_unit::Bool = true)
    fr = _fusion_rule(rules)
    require_multiplicity_free(fr)
    fdict = _fmap(fvars); rdict = _rmap(rvars)
    eqs = PolynomialEquation[]
    for fv in fvars
        for g in fusion_channels(fr, fv.a, fv.c)
            fr.N[fv.b, g, fv.d] == 0 && continue
            lf = _rv(rdict, fv.a, fv.f, fv.d)
            r1 = _rv(rdict, fv.a, fv.b, fv.e)
            r2 = _rv(rdict, fv.a, fv.c, g)
            f2 = _fv(fdict, fv.b, fv.a, fv.c, fv.d, fv.e, g)
            (lf === nothing || r1 === nothing || r2 === nothing || f2 === nothing) && continue
            push!(eqs, _equation(_vexpr(lf) * _vexpr(fv),
                                 _vexpr(r1) * _vexpr(r2) * _vexpr(f2);
                                 metadata = Dict(:kind => :left_hexagon,
                                                 :indices => fv.variable.indices,
                                                 :channel => g)))
        end
    end
    if normalize_unit
        for v in unit_rsymbol_normalizations(fr, rvars)
            push!(eqs, _equation(_vexpr(v), _one_expr();
                metadata = Dict(:kind => :unit_R_normalization, :indices => v.variable.indices)))
        end
    end
    return _unique_equations(eqs)
end

function right_hexagon_equations(rules, fvars = fsymbol_variables(rules),
                                 rvars = rsymbol_variables(rules);
                                 normalize_unit::Bool = true)
    fr = _fusion_rule(rules)
    require_multiplicity_free(fr)
    fdict = _fmap(fvars); rdict = _rmap(rvars)
    eqs = PolynomialEquation[]
    for fv in fvars
        for g in fusion_channels(fr, fv.a, fv.c)
            fr.N[g, fv.b, fv.d] == 0 && continue
            lf = _rv(rdict, fv.e, fv.c, fv.d)
            r1 = _rv(rdict, fv.b, fv.c, fv.f)
            r2 = _rv(rdict, fv.a, fv.c, g)
            f2 = _fv(fdict, fv.a, fv.c, fv.b, fv.d, g, fv.f)
            (lf === nothing || r1 === nothing || r2 === nothing || f2 === nothing) && continue
            push!(eqs, _equation(_vexpr(lf) * _vexpr(fv),
                                 _vexpr(r1) * _vexpr(r2) * _vexpr(f2);
                                 metadata = Dict(:kind => :right_hexagon,
                                                 :indices => fv.variable.indices,
                                                 :channel => g)))
        end
    end
    if normalize_unit
        for v in unit_rsymbol_normalizations(fr, rvars)
            push!(eqs, _equation(_vexpr(v), _one_expr();
                metadata = Dict(:kind => :unit_R_normalization, :indices => v.variable.indices)))
        end
    end
    return _unique_equations(eqs)
end

hexagon_equations(rules::Union{FusionRule, Array{Int, 3}},
                  fvars::Vector{FSymbolVariable},
                  rvars::Vector{RSymbolVariable}; kwargs...) =
    _unique_equations(vcat(left_hexagon_equations(rules, fvars, rvars; kwargs...),
                          right_hexagon_equations(rules, fvars, rvars; kwargs...)))

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
    fvars = fsymbol_variables(fr)
    rvars = rsymbol_variables(fr)
    eqs = PolynomialEquation[]
    include_pentagon && append!(eqs, pentagon_equations(fr, fvars;
                                                        normalize_unit = normalize_unit))
    include_hexagon && append!(eqs, hexagon_equations(fr, fvars, rvars;
                                                      normalize_unit = normalize_unit))
    variables = vcat([v.variable for v in fvars], [v.variable for v in rvars])
    assumptions = ["multiplicity-free fusion rule", "unit object has label 1",
                   "hexagon equations use scalar multiplicity-free channels"]
    metadata = Dict{Symbol, Any}(:rank => fr.rank,
                                 :f_variables => length(fvars),
                                 :r_variables => length(rvars),
                                 :include_pentagon => include_pentagon,
                                 :include_hexagon => include_hexagon,
                                 :normalizations => normalize_unit)
    return FREquationSystem(fr, fvars, rvars, variables, _unique_equations(eqs),
                            metadata, assumptions, base_ring)
end

function validate_fr_system(system::FREquationSystem)
    require_multiplicity_free(system.rules)
    vars = Set(system.variables)
    for eq in system.equations, side in (eq.lhs, eq.rhs), term in side.terms, p in term.powers
        p.first in vars || error("equation references variable $(p.first.name) not present in system")
        p.second > 0 || error("negative or zero exponents are not supported in polynomial equations")
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

function gauge_transform_fsymbol(var::FSymbolVariable,
                                 gvars::Vector{GaugeVariable} = GaugeVariable[])
    gdict = _gauge_map(gvars)
    factor = _one_expr()
    for (ch, exp) in (((var.a, var.b, var.e), 1),
                      ((var.e, var.c, var.d), 1),
                      ((var.b, var.c, var.f), -1),
                      ((var.a, var.f, var.d), -1))
        exp < 0 && continue
        gv = get(gdict, ch, nothing)
        gv === nothing && continue
        factor = factor * _vexpr(gv)
    end
    # Negative exponents are recorded as metadata rather than forced into a
    # polynomial expression; solver backends can introduce inverse variables.
    return (symbol = var, transformed = _vexpr(var) * factor,
            denominator = [ch for (ch, exp) in (((var.a, var.b, var.e), 1),
                                               ((var.e, var.c, var.d), 1),
                                               ((var.b, var.c, var.f), -1),
                                               ((var.a, var.f, var.d), -1)) if exp < 0])
end

function gauge_transform_rsymbol(var::RSymbolVariable,
                                 gvars::Vector{GaugeVariable} = GaugeVariable[])
    gdict = _gauge_map(gvars)
    numerator = get(gdict, (var.b, var.a, var.c), nothing)
    denominator = get(gdict, (var.a, var.b, var.c), nothing)
    transformed = _vexpr(var)
    numerator === nothing || (transformed = transformed * _vexpr(numerator))
    return (symbol = var, transformed = transformed,
            denominator = denominator === nothing ? Tuple{Int,Int,Int}[] : [(var.a, var.b, var.c)])
end

function gauge_transform(system::FREquationSystem)
    gvars = gauge_variables(system.rules)
    f_actions = [gauge_transform_fsymbol(v, gvars) for v in system.fvars]
    r_actions = [gauge_transform_rsymbol(v, gvars) for v in system.rvars]
    return (gauge_variables = gvars, F = f_actions, R = r_actions)
end

function gauge_fix(system::FREquationSystem; strategy::Symbol = :safe)
    strategy == :safe || error("only strategy=:safe is implemented for v0.8 gauge_fix")
    fixed_eqs = copy(system.equations)
    fixed = Vector{NamedTuple}()
    for v in unit_fsymbol_normalizations(system.rules, system.fvars)
        push!(fixed_eqs, _equation(_vexpr(v), _one_expr();
            metadata = Dict(:kind => :safe_gauge_fix, :reason => :unit_F,
                            :indices => v.variable.indices)))
        push!(fixed, (kind = :F, reason = :unit, indices = v.variable.indices))
    end
    for v in unit_rsymbol_normalizations(system.rules, system.rvars)
        push!(fixed_eqs, _equation(_vexpr(v), _one_expr();
            metadata = Dict(:kind => :safe_gauge_fix, :reason => :unit_R,
                            :indices => v.variable.indices)))
        push!(fixed, (kind = :R, reason = :unit, indices = v.variable.indices))
    end
    meta = copy(system.metadata)
    meta[:gauge_fix_strategy] = :safe
    meta[:fixed_symbols] = fixed
    meta[:residual_gauge_variables] = [v.variable.name for v in gauge_variables(system.rules)
                                       if v.a != 1 && v.b != 1]
    return FREquationSystem(system.rules, system.fvars, system.rvars, system.variables,
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
    eqs = [PolynomialEquation(_reduce_expr_mod_p(eq.lhs, pp),
                              _reduce_expr_mod_p(eq.rhs, pp),
                              copy(eq.metadata)) for eq in system.equations]
    meta = copy(system.metadata)
    meta[:base_field] = Symbol("F_$pp")
    if haskey(meta, :conductor)
        meta[:p_mod_conductor] = mod(pp, meta[:conductor])
    end
    return FiniteFieldEquationSystem(system, pp, copy(system.variables), eqs, meta)
end

function solve_finite_field(system::FiniteFieldEquationSystem; kwargs...)
    error("solve_finite_field is not implemented for general FR systems; v0.8 only provides reduction infrastructure")
end

function cyclotomic_reconstruct(fp_solution; conductor::Integer)
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
    N[2, 3, 3] = 1; N[3, 2, 3] = 1
    N[3, 3, 1] = 1; N[3, 3, 2] = 1
    return FusionRule(N)
end
