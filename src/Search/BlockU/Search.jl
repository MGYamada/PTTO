# ============================================================
#  Step 6: Top-level single-prime driver
# ============================================================

function _candidate_from_fixed_S(S_atomic::Matrix{Int}, T_atomic::Vector{Int},
                                 p::Int, U_params;
                                 verlinde_threshold::Int)
    r = size(S_atomic, 1)
    result = verlinde_find_unit(S_atomic, p; threshold = verlinde_threshold)
    result === nothing && return MTCCandidate[]
    (u, N_tensor) = result
    Suu_inv = invmod(S_atomic[u, u], p)
    d = [mod(S_atomic[u, i] * Suu_inv, p) for i in 1:r]
    D2 = sum(mod(d[i] * d[i], p) for i in 1:r) % p
    return [MTCCandidate(p, U_params, copy(S_atomic), copy(T_atomic),
                         u, N_tensor, d, D2)]
end

"""
    find_mtcs_at_prime(catalog::Vector{AtomicIrrep}, stratum::Stratum,
                       p::Int; verlinde_threshold::Int = 3,
                       max_block_dim::Int = 3,
                       search_mode::Symbol = :groebner,
                       max_units_for_groebner::Int = typemax(Int),
                       groebner_allow_fallback::Bool = false,
                       precheck_unit_axiom::Bool = true)
        -> Vector{MTCCandidate}

Top-level Phase 2 driver at a single prime.

Steps:
1. Build block-diagonal atomic (S, T) from stratum + catalog.
2. Reduce to F_p (requires N | p-1 where N = catalog[].N).
3. Decompose T into eigenspaces; compute parameter_dim.
4. For each degenerate eigenspace of dimension n_θ ≥ 2, solve for
   block-U candidates using the selected backend, then take the Cartesian
   product across all degenerate eigenspaces.
5. For each block-U candidate, apply transformation to S, then check
   Verlinde integrality.
6. Return all MTC candidates found.

`max_block_dim` is a safety guard against very large block systems: if any
degenerate eigenspace has n_θ > max_block_dim, an error is raised.

`search_mode` selects the block-U candidate backend:
- `:groebner` (default): Gröbner basis plus finite-field point extraction
- `:exhaustive`: Cayley+reflection sweep

In `:groebner` mode, the driver performs best-effort Gröbner system
construction for fixed-unit variants and attempts direct U-block point
extraction. If that produces no candidates, it next tries Gröbner point
extraction for the orthogonality block itself.

`max_units_for_groebner` can cap how many unit indices are used for
fixed-unit Gröbner extraction (default: all).

`groebner_allow_fallback` controls whether `:groebner` mode is allowed to
fall back to the explicit `:exhaustive` Cayley/reflection backend when
solver extraction returns no candidates. It is disabled by default.

`precheck_unit_axiom` enables a fast unit-axiom prefilter before full
`verlinde_find_unit` tensor construction.
"""
function find_mtcs_at_prime(catalog::Vector{AtomicIrrep}, stratum::Stratum,
                            p::Int; verlinde_threshold::Int = 3,
                            max_block_dim::Int = 3,
                            search_mode::Symbol = :groebner,
                            max_units_for_groebner::Int = typemax(Int),
                            groebner_allow_fallback::Bool = false,
                            precheck_unit_axiom::Bool = true)
    validate_search_mode(search_mode)
    max_units_for_groebner >= 1 || error("max_units_for_groebner must be ≥ 1")
    # Common N
    N = catalog[first(keys(stratum.multiplicities))].N

    # Build block-diagonal atomic data
    S_K, T_K, K, blocks = build_block_diagonal(catalog, stratum)

    # Reduce to F_p
    zeta_N_Fp = find_zeta_in_Fp(N, p)
    S_atomic = [cyclotomic_to_Fp(S_K[i, j], zeta_N_Fp, p) for i in 1:size(S_K, 1), j in 1:size(S_K, 2)]
    T_atomic = [cyclotomic_to_Fp(t, zeta_N_Fp, p) for t in T_K]

    # T-eigenspace decomposition
    eigenspaces = t_eigenspace_decomposition(T_atomic, p)
    p_dim = parameter_dim(eigenspaces)

    degenerate = sort([(theta, indices) for (theta, indices) in eigenspaces if length(indices) >= 2];
                      by = x -> x[1])

    # Case 1: no degeneracy → atomic S is already an MTC candidate (modulo signs)
    if isempty(degenerate)
        return _candidate_from_fixed_S(S_atomic, T_atomic, p, :atomic;
                                       verlinde_threshold = verlinde_threshold)
    end

    if _orthogonal_block_product_fixes_S(S_atomic, degenerate, p)
        return _candidate_from_fixed_S(S_atomic, T_atomic, p, :block_invariant;
                                       verlinde_threshold = verlinde_threshold)
    end

    # Case 2+: one or more degenerate eigenspaces.
    block_candidate_sets = Vector{Vector{NamedTuple}}()
    use_verlinde_warmup = length(degenerate) == 1
    for (theta_deg, indices_deg) in degenerate
        n_block = length(indices_deg)
        n_block <= max_block_dim || error(
            "Degenerate eigenspace dim n=$n_block exceeds max_block_dim=$max_block_dim. " *
            "Increase max_block_dim only if the selected block-U backend can handle it.")

        U_blocks = _block_U_candidates_for_eigenspace(S_atomic, indices_deg, p;
                                                      search_mode = search_mode,
                                                      max_units_for_groebner = max_units_for_groebner,
                                                      groebner_allow_fallback = groebner_allow_fallback,
                                                      use_verlinde_warmup = use_verlinde_warmup)
        isempty(U_blocks) && return MTCCandidate[]
        push!(block_candidate_sets,
              [(theta = theta_deg, indices = indices_deg, U = U_block) for U_block in U_blocks])
    end

    candidates = MTCCandidate[]
    r = size(S_atomic, 1)
    for block_tuple in Iterators.product(block_candidate_sets...)
        block_choices = collect(block_tuple)
        S_prime = _apply_block_U_product(S_atomic, block_choices, p)
        if precheck_unit_axiom
            any_unit = any(u -> passes_unit_axiom(S_prime, p, u), 1:r)
            any_unit || continue
        end
        result = verlinde_find_unit(S_prime, p; threshold = verlinde_threshold)
        if result !== nothing
            (u_idx, N_tensor) = result
            Suu_inv = invmod(S_prime[u_idx, u_idx], p)
            d = [mod(S_prime[u_idx, i] * Suu_inv, p) for i in 1:r]
            D2 = sum(mod(d[m] * d[m], p) for m in 1:r) % p
            U_params = length(block_choices) == 1 ? block_choices[1].U : block_choices
            push!(candidates, MTCCandidate(
                p, U_params,
                copy(S_prime), copy(T_atomic),
                u_idx, N_tensor, d, D2))
        end
    end
    return candidates
end
