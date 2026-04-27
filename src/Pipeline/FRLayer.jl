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

function _phase4_removed_error()
    error("exact Phase 4 did not find a cyclotomic F/R solution for this input")
end

function _select_fr_for_st(candidates, Nijk, S_cyc, T_cyc, N)
    isempty(candidates) && _phase4_removed_error()
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
            score = merge(score, (galois_exponent = a,
                                  order_key = (score.ok ? 0 : 1,
                                               string(score.T_max),
                                               string(score.S_max),
                                               ci,
                                               a)))
            push!(all_scores, score)
            if best_score === nothing || score.order_key < best_score.order_key
                best_score = score
                best_idx = ci
                best_candidate = (F = F, R = R, report = cand.report)
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
    for perm in automorphisms
        S_diffs = [S_from_R[perm[i], perm[j]] - _matrix_entry(S_target, i, j)
                   for i in 1:r, j in 1:r]
        T_diffs = [T_from_R[perm[i]] - T_target[i] for i in 1:r]
        S_ok = all(iszero, S_diffs)
        T_ok = all(iszero, T_diffs)
        S_err = S_ok ? zero(K) : first(x for x in S_diffs if !iszero(x))
        T_err = T_ok ? zero(K) : first(x for x in T_diffs if !iszero(x))
        score = (ok = S_ok && T_ok,
                 S_max = S_err,
                 T_max = T_err,
                 best_perm = perm,
                 S_roundtrip = S_from_R,
                 T_roundtrip = T_from_R)
        if best === nothing || (score.ok && !best.ok)
            best = score
        end
    end
    return best
end

function _score_fr_st_match(F_values::Vector,
                            R_values::Vector,
                            Nijk::Array{Int,3},
                            S_target,
                            T_target,
                            N::Int;
                            candidate_index::Int)
    md = _modular_data_roundtrip(F_values, R_values, Nijk, S_target, T_target, N)
    return merge(md, (candidate_index = candidate_index,
                      order_key = (md.ok ? 0 : 1,
                                   string(md.T_max),
                                   string(md.S_max),
                                   candidate_index)))
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

function _select_phase4_candidate(candidates::Vector{NamedTuple},
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
    return merge(candidates[best.candidate_index], (report = best,))
end

# ============================================================
#  compute_FR_from_ST: exact Phase 4 over Q(ζ_N)
# ============================================================

function compute_FR_from_ST(Nijk::Array{Int,3};
                            context = nothing,
                            conductor = nothing,
                            N = nothing,
                            S = nothing,
                            T = nothing,
                            return_all::Bool = false,
                            primes::Vector{Int} = [101, 103, 107, 109],
                            verbose::Bool = false,
                            kwargs...)
    isempty(kwargs) || error("unsupported keyword arguments: $(collect(keys(kwargs)))")
    ctx = _default_context_from_kwargs(context = context, conductor = conductor, N = N)
    r = size(Nijk, 1)
    _, pentagon_eqs, nF = get_pentagon_system(Nijk, r)
    F_solutions = solve_pentagon_modular_crt(pentagon_eqs, nF;
                                             Nijk = Nijk,
                                             context = ctx,
                                             primes = primes,
                                             show_progress = verbose)
    candidates = NamedTuple[]
    for F in F_solutions
        _, hex_eqs, nR = get_hexagon_system(Nijk, r, F; context = ctx)
        R_solutions = solve_hexagon_modular_crt(hex_eqs, nR;
                                                Nijk = Nijk,
                                                context = ctx,
                                                primes = primes,
                                                show_progress = verbose)
        for R in R_solutions
            push!(candidates, (F = F, R = R, report = nothing))
        end
    end
    isempty(candidates) && _phase4_removed_error()
    selected = _select_phase4_candidate(candidates, Nijk, ctx, S, T)
    return (F = selected.F,
            R = selected.R,
            report = selected.report,
            candidates = candidates)
end
