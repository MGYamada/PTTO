"""
Exact F/R layer orchestration helpers for the pipeline.

These helpers solve, score, and attach F/R data to already lifted modular
data while keeping the main pipeline focused on conductor orchestration.
"""

_is_reconstruction_unstable_message(msg::AbstractString) = begin
    low = lowercase(msg)
    return occursin("reconstruct", low) ||
           occursin("inconsistent", low) ||
           occursin("fresh-prime", low) ||
           occursin("crt", low)
end

function _normalize_toric_gauge_mode(gauge_fixing::Symbol,
                                     toric_gauge_fixing::Bool)
    gauge_fixing in (:auto, :toric, :general, :none) ||
        error("gauge_fixing must be one of :auto, :toric, :general, or :none; got $(repr(gauge_fixing))")
    !toric_gauge_fixing && gauge_fixing == :auto && return :none
    !toric_gauge_fixing && gauge_fixing == :toric &&
        error("conflicting gauge-fixing options: toric_gauge_fixing=false but gauge_fixing=:toric")
    return gauge_fixing
end

function _toric_pre_reconstruction_data(Nijk::Array{Int,3};
                                        gauge_fixing::Symbol = :auto,
                                        toric_gauge_fixing::Bool = true,
                                        field = nothing,
                                        conventions = :tensorcategories)
    mode = _normalize_toric_gauge_mode(gauge_fixing, toric_gauge_fixing)
    general = general_gauge_data(Nijk; field = field, conventions = conventions)
    mode == :none &&
        return (apply = false, mode = mode, data = nothing,
                reason = :disabled, general = general)
    if !is_toric(general)
        mode == :toric &&
            throw(ToricGaugeFixingError("gauge_fixing=:toric was requested, but the GeneralGauge data has non-GL(1) factors"))
        return (apply = false, mode = mode, data = nothing,
                reason = :not_multiplicity_free, general = general)
    end
    data = toric_gauge_data(Nijk; field = field, conventions = conventions)
    return (apply = mode != :general, mode = mode, data = data,
            reason = mode == :general ? :general_metadata_only : :multiplicity_free,
            general = general)
end

function _fr_reconstruction_no_solution_error()
    error("exact F/R reconstruction did not find a cyclotomic solution for this input")
end

function _select_fr_for_st(candidates, Nijk, S_cyc, T_cyc, N)
    isempty(candidates) && _fr_reconstruction_no_solution_error()
    all_scores = NamedTuple[]
    best_idx = 0
    best_score = nothing
    best_candidate = nothing
    ctx = CyclotomicContext(N)
    galois_units = [a for a in 1:N if gcd(a, N) == 1]
    for (ci, cand) in enumerate(candidates)
        for a in galois_units
            F = a == 1 ? cand.F : galois_action(ctx, cand.F, a)
            R = a == 1 ? cand.R : galois_action(ctx, cand.R, a)
            score = _score_fr_st_match(F, R, Nijk, S_cyc, T_cyc, N;
                                       candidate_index = ci)
            report = _with_fr_roundtrip_metadata(score.report;
                                                 candidate_index = ci,
                                                 galois_exponent = a)
            score = merge(score, (galois_exponent = a,
                                  report = report,
                                  order_key = (score.order_key..., a)))
            push!(all_scores, score)
            if best_score === nothing || score.order_key < best_score.order_key
                best_score = score
                best_idx = ci
                best_candidate = (F = F, R = R, report = report)
            end
        end
    end
    return (selected = best_candidate,
            score = best_score,
            selected_index = best_idx,
            selected_ok = best_score.ok,
            all_scores = all_scores)
end

function _normalize_T_by_unit(T_cyc::AbstractVector)
    isempty(T_cyc) && return collect(T_cyc)
    t0 = T_cyc[1]
    iszero(t0) && return collect(T_cyc)
    return [t / t0 for t in T_cyc]
end


function _matrix_dim(M)
    return M isa MatElem ? (nrows(M), ncols(M)) : size(M)
end

_matrix_entry(M, i::Int, j::Int) = M[i, j]

function _forward_r_var_count(Nijk::Array{Int,3})
    r = size(Nijk, 1)
    return sum(Nijk[i, j, k]^2 for i in 1:r, j in 1:r, k in 1:r)
end

function _forward_R_values(R_values::Vector, Nijk::Array{Int,3})
    n = _forward_r_var_count(Nijk)
    length(R_values) == n && return R_values
    length(R_values) == 2n && return R_values[1:n]
    error("R_values has length $(length(R_values)); expected $n or $(2n)")
end

function _extract_R_block_exact(R_values::Vector, Nijk::Array{Int,3},
                                i::Int, j::Int, k::Int)
    n = Nijk[i, j, k]
    n == 0 && return nothing
    positions, _ = _braiding_block_positions(Nijk)
    pos = positions[(i, j, k)]
    K = parent(R_values[1])
    M = zero_matrix(K, n, n)
    for a in 1:n, b in 1:n
        M[a, b] = R_values[pos[(a - 1) * n + b]]
    end
    return M
end

function _trace_exact(M)
    n = nrows(M)
    t = zero(base_ring(M))
    for i in 1:n
        t += M[i, i]
    end
    return t
end

function _fusion_automorphisms_fixing_unit(Nijk::Array{Int,3})
    r = size(Nijk, 1)
    r == 1 && return [Int[1]]
    autos = Vector{Vector{Int}}()
    for rest in _all_permutations(collect(2:r))
        perm = vcat(1, rest)
        ok = true
        @inbounds for i in 1:r, j in 1:r, k in 1:r
            if Nijk[perm[i], perm[j], perm[k]] != Nijk[i, j, k]
                ok = false
                break
            end
        end
        ok && push!(autos, perm)
    end
    isempty(autos) && push!(autos, collect(1:r))
    return autos
end

function _modular_data_from_FR(R_values::Vector,
                               Nijk::Array{Int,3},
                               S_target)
    r = size(Nijk, 1)
    dims = _matrix_dim(S_target)
    dims == (r, r) || error("S_target has shape $dims; expected ($r, $r)")

    K = parent(_matrix_entry(S_target, 1, 1))
    R_fwd = _forward_R_values(R_values, Nijk)
    d = [_matrix_entry(S_target, a, 1) / _matrix_entry(S_target, 1, 1) for a in 1:r]
    D = inv(_matrix_entry(S_target, 1, 1))

    S_from_R = zero_matrix(K, r, r)
    for a in 1:r, b in 1:r
        acc = zero(K)
        for c in 1:r
            Nijk[a, b, c] == 0 && continue
            Rab = _extract_R_block_exact(R_fwd, Nijk, a, b, c)
            Rba = _extract_R_block_exact(R_fwd, Nijk, b, a, c)
            acc += d[c] * _trace_exact(Rba * Rab)
        end
        S_from_R[a, b] = acc / D
    end

    T_from_R = Vector{typeof(K(1))}(undef, r)
    for a in 1:r
        acc = zero(K)
        for c in 1:r
            Nijk[a, a, c] == 0 && continue
            Raa = _extract_R_block_exact(R_fwd, Nijk, a, a, c)
            acc += d[c] * _trace_exact(Raa)
        end
        T_from_R[a] = acc / d[a]
    end
    if !iszero(T_from_R[1])
        t0 = T_from_R[1]
        T_from_R = [t / t0 for t in T_from_R]
    end
    return (S = S_from_R, T = T_from_R)
end

function _modular_data_roundtrip(F_values::Vector,
                                 R_values::Vector,
                                 Nijk::Array{Int,3},
                                 S_target,
                                 T_target,
                                 N::Int)
    r = size(Nijk, 1)
    length(T_target) == r || error("T_target has length $(length(T_target)); expected $r")
    md_from_fr = _modular_data_from_FR(R_values, Nijk, S_target)
    S_from_R = md_from_fr.S
    T_from_R = md_from_fr.T
    K = parent(_matrix_entry(S_target, 1, 1))

    automorphisms = _fusion_automorphisms_fixing_unit(Nijk)
    best = nothing
    best_key = nothing
    for perm in automorphisms
        S_diffs = [S_from_R[perm[i], perm[j]] - _matrix_entry(S_target, i, j)
                   for i in 1:r, j in 1:r]
        T_diffs = [T_from_R[perm[i]] - T_target[i] for i in 1:r]
        S_ok = all(iszero, S_diffs)
        T_ok = all(iszero, T_diffs)
        S_mismatches = count(!iszero, S_diffs)
        T_mismatches = count(!iszero, T_diffs)
        S_err = S_ok ? zero(K) : first(x for x in S_diffs if !iszero(x))
        T_err = T_ok ? zero(K) : first(x for x in T_diffs if !iszero(x))
        score = FRRoundtripReport(ok = S_ok && T_ok,
                                  S_error = S_err,
                                  T_error = T_err,
                                  best_perm = perm,
                                  S_roundtrip = S_from_R,
                                  T_roundtrip = T_from_R)
        key = (S_mismatches,
               T_mismatches,
               string(S_err),
               string(T_err),
               perm)
        if best === nothing || key < best_key
            best = score
            best_key = key
        end
    end
    return best
end

"""
    _fr_roundtrip_attachable(report)

Pipeline acceptance predicate for exact `(F, R)` data. A candidate is
attachable only when the reconstructed modular data matches both target
`S` and target `T` exactly.
"""
_fr_roundtrip_attachable(report::FRRoundtripReport) = report.ok

function _score_fr_st_match(F_values::Vector,
                            R_values::Vector,
                            Nijk::Array{Int,3},
                            S_target,
                            T_target,
                            N::Int;
                            candidate_index::Int)
    md = _modular_data_roundtrip(F_values, R_values, Nijk, S_target, T_target, N)
    return (ok = md.ok,
            S_error = md.S_error,
            T_error = md.T_error,
            S_max = md.S_max,
            T_max = md.T_max,
            best_perm = md.best_perm,
            S_roundtrip = md.S_roundtrip,
            T_roundtrip = md.T_roundtrip,
            report = _with_fr_roundtrip_metadata(md;
                                                 candidate_index = candidate_index),
            candidate_index = candidate_index,
            order_key = (md.ok ? 0 : (iszero(md.S_error) ? 1 : 2),
                         string(md.S_max),
                         string(md.T_max),
                         candidate_index))
end

function _twists_from_T_arg(T)
    T === nothing && return nothing
    try
        r = size(T, 1)
        size(T, 2) == r && return [T[i, i] for i in 1:r]
    catch
    end
    return collect(T)
end

function _select_fr_candidate(candidates::Vector{NamedTuple},
                              Nijk::Array{Int,3},
                              ctx::CyclotomicContext,
                              S_arg, T_arg)
    S_target = S_arg
    T_target = _twists_from_T_arg(T_arg)
    (S_target === nothing || T_target === nothing) && return candidates[1]

    scored = [_score_fr_st_match(c.F, c.R, Nijk, S_target, T_target, ctx.N;
                                 candidate_index = i)
              for (i, c) in enumerate(candidates)]
    best = first(sort(scored; by = s -> s.order_key))
    return merge(candidates[best.candidate_index], (report = best.report,))
end

# ============================================================
#  compute_FR_from_ST: exact F/R reconstruction over Q(ζ_N)
# ============================================================

function compute_FR_from_ST(Nijk::Array{Int,3};
                            context = nothing,
                            conductor = nothing,
                            N = nothing,
                            S = nothing,
                            T = nothing,
                            return_all::Bool = false,
                            primes::Vector{Int} = [101, 103, 107, 109],
                            max_solutions::Int = 32,
                            max_modular_solutions::Int = max(1024, 16 * max_solutions),
                            reconstruction_bound::Int = 4,
                            denominator_bound::Int = 4,
                            max_crt_tuples::Int = 4096,
                            max_ambiguous_crt_coords::Int = 4096,
                            exact_fallback::Bool = true,
                            pentagon_gauge_fixing::Bool = false,
                            toric_gauge_data = nothing,
                            canonicalize_gauge::Bool = true,
                            verbose::Bool = false,
                            kwargs...)
    isempty(kwargs) || error("unsupported keyword arguments: $(collect(keys(kwargs)))")
    ctx = _default_context_from_kwargs(context = context, conductor = conductor, N = N)
    r = size(Nijk, 1)
    if r == 1
        K = field(ctx)
        F = elem_type(K)[]
        R = elem_type(K)[one(K), one(K)]
        candidates = NamedTuple[(F = F, R = R, report = nothing)]
        selected = _select_fr_candidate(candidates, Nijk, ctx, S, T)
        fixed = canonicalize_gauge ? canonical_gauge(selected.F, selected.R, Nijk) :
                (F = selected.F, R = selected.R, gauge = nothing)
        result = (F = fixed.F,
                  R = fixed.R,
                  report = selected.report)
        return return_all ? merge(result, (candidates = candidates,)) : result
    end
    _, pentagon_eqs, nF = get_pentagon_system(Nijk, r)
    gauge_fixed = if pentagon_gauge_fixing
        if toric_gauge_data !== nothing
            Int[i for i in toric_gauge_data.slice.fixed_indices if i <= nF]
        else
            _select_pentagon_gauge_fixed_indices(Nijk, r, nF)
        end
    else
        Int[]
    end
    reduced_pentagon = _substitute_fixed_one_polys(pentagon_eqs, nF, gauge_fixed;
                                                   var_prefix = :x)
    verbose && !isempty(gauge_fixed) &&
        println("  pentagon gauge fixing: fixed $(length(gauge_fixed)) / $nF variables; " *
                "reduced nF=$(reduced_pentagon.n)")

    F_reduced_solutions = solve_pentagon_modular_crt(reduced_pentagon.eqs,
                                                     reduced_pentagon.n;
                                             Nijk = Nijk,
                                             context = ctx,
                                             primes = primes,
                                             max_solutions = max_solutions,
                                             max_modular_solutions = max_modular_solutions,
                                             reconstruction_bound = reconstruction_bound,
                                             denominator_bound = denominator_bound,
                                             max_crt_tuples = max_crt_tuples,
                                             max_ambiguous_crt_coords = max_ambiguous_crt_coords,
                                             exact_fallback = exact_fallback,
                                             show_progress = verbose)
    F_solutions = [_expand_fixed_one_solution(F, nF, reduced_pentagon.free_indices, field(ctx))
                   for F in F_reduced_solutions]
    candidates = NamedTuple[]
    for F in F_solutions
        _, hex_eqs, nR = get_hexagon_system(Nijk, r, F; context = ctx)
        R_solutions = solve_hexagon_modular_crt(hex_eqs, nR;
                                                Nijk = Nijk,
                                                context = ctx,
                                                primes = primes,
                                                max_solutions = max_solutions,
                                                max_modular_solutions = max_modular_solutions,
                                                reconstruction_bound = reconstruction_bound,
                                                denominator_bound = denominator_bound,
                                                max_crt_tuples = max_crt_tuples,
                                                max_ambiguous_crt_coords = max_ambiguous_crt_coords,
                                                exact_fallback = exact_fallback,
                                                show_progress = verbose)
        for R in R_solutions
            push!(candidates, (F = F, R = R, report = nothing))
        end
    end
    isempty(candidates) && _fr_reconstruction_no_solution_error()
    selected = _select_fr_candidate(candidates, Nijk, ctx, S, T)
    fixed = canonicalize_gauge ? canonical_gauge(selected.F, selected.R, Nijk) :
            (F = selected.F, R = selected.R, gauge = nothing)
    result = (F = fixed.F,
              R = fixed.R,
              report = selected.report)
    return return_all ? merge(result, (candidates = candidates,)) : result
end
