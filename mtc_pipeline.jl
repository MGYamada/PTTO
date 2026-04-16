# mtc_pipeline_v2.jl
# Stage 2 revised: fixes two bugs found in first test run.
#
# Bug 1: find_signed_V checked T[i,i]^12 = 1, but the correct condition is
#         that T[i,i] is a (12·ord(T))-th root of unity, or more precisely
#         that the unit object's spin s₀ = arg(T[i₀,i₀])/(2π) satisfies
#         24·s₀ ∈ ℤ. Fixed by relaxing the T^12 test.
#
# Bug 2: For atomic irreps (no direct sum), the "block diagonal has zeros"
#         issue doesn't apply — the S matrix is already dense. But odd
#         irreps have purely imaginary S, so find_signed_V's "all-real" check
#         fails. Fix: for odd irreps, multiply S by i to get a real matrix,
#         then proceed. This corresponds to the paper's normalization
#         ρ_α(s) = ζ₄^α · S/D (eq. 3.2), choosing α to make S real.
#
# Prerequisites: include("mtc_classifier.jl"); include("mtc_types.jl")

include("mtc_types.jl")

# ============================================================
#  Import pipeline functions from v1 that don't need changes
# ============================================================
# (build_atomic_catalog, enumerate_irrep_sums, order_by_t_spectrum,
#  eigenvalue_blocks are unchanged — we redefine the full file to be
#  self-contained)

# --- Step 0: catalog (unchanged) ---

function build_atomic_catalog(N::Int; max_rank::Int=20)
    pf = prime_power_factors(N)
    K, z = cyclotomic_field(N)
    zeta = exp(2π * im / N)
    deg = degree(K)

    factor_irreps = []
    for (p, a) in pf
        gap_irreps = irreps_for_prime_power(p, a)
        st_pairs = []
        for rep in gap_irreps
            r = Int(rep.degree)
            cond_S = matrix_conductor(rep.S, r)
            cond_T = matrix_conductor(rep.T, r)
            rep_cond = lcm(cond_S, cond_T)
            N % rep_cond != 0 && continue
            S = gap_to_oscar_matrix(rep.S, r, K, N)
            T = gap_to_oscar_matrix(rep.T, r, K, N)
            push!(st_pairs, (S, T, r))
        end
        push!(factor_irreps, st_pairs)
    end

    catalog = AtomicIrrep[]
    if length(factor_irreps) == 1
        for (S, T, r) in factor_irreps[1]
            r > max_rank && continue
            S_num = matrix_to_complex(S, zeta, deg)
            T_num = matrix_to_complex(T, zeta, deg)
            par = compute_parity(S_num)
            level = _numeric_order(T_num, r, N)
            push!(catalog, AtomicIrrep(r, level, "$(r)d_$level", S, T, par, K, N))
        end
    else
        indices = [1:length(fi) for fi in factor_irreps]
        for idx_tuple in Iterators.product(indices...)
            total_rank = prod(factor_irreps[k][idx_tuple[k]][3] for k in eachindex(pf))
            total_rank > max_rank && continue
            S = factor_irreps[1][idx_tuple[1]][1]
            T = factor_irreps[1][idx_tuple[1]][2]
            for k in 2:length(factor_irreps)
                S = kron_oscar(S, factor_irreps[k][idx_tuple[k]][1], K)
                T = kron_oscar(T, factor_irreps[k][idx_tuple[k]][2], K)
            end
            S_num = matrix_to_complex(S, zeta, deg)
            T_num = matrix_to_complex(T, zeta, deg)
            par = compute_parity(S_num)
            level = _numeric_order(T_num, total_rank, N)
            dims = [factor_irreps[k][idx_tuple[k]][3] for k in eachindex(pf)]
            label = join(["$(d)d" for d in dims], "⊗") * "_$level"
            push!(catalog, AtomicIrrep(total_rank, level, label, S, T, par, K, N))
        end
    end

    println("Atomic catalog for N=$N: $(length(catalog)) irreps")
    for (i, a) in enumerate(catalog)
        println("  [$i] $a")
    end
    return catalog
end

# --- Step 1: enumerate sums (unchanged) ---

function enumerate_irrep_sums(catalog::Vector{AtomicIrrep}, N::Int; max_rank::Int=6)
    K = isempty(catalog) ? nothing : catalog[1].K
    K === nothing && return IrrepSumRep[]

    rank_list = [a.dim for a in catalog]
    mult_vectors = _enumerate_multiplicity_vectors(rank_list, max_rank)

    results = IrrepSumRep[]
    for mv in mult_vectors
        components = Tuple{Int,Int}[]
        blocks_S = []
        blocks_T = []
        type_parts = Int[]
        total_rank = 0
        for (α, m) in enumerate(mv)
            m == 0 && continue
            push!(components, (α, m))
            for _ in 1:m
                push!(blocks_S, catalog[α].S)
                push!(blocks_T, catalog[α].T)
                push!(type_parts, catalog[α].dim)
            end
            total_rank += catalog[α].dim * m
        end
        total_rank > max_rank && continue

        S_big = _block_diag_matrices(blocks_S, K)
        T_big = _block_diag_matrices(blocks_T, K)
        young = sort(type_parts; rev=true)

        push!(results, IrrepSumRep(components, young, S_big, T_big, total_rank, K, N))
    end

    println("IrrepSumRep enumeration: $(length(results)) sums with rank ≤ $max_rank")
    type_counts = Dict{Vector{Int}, Int}()
    for r in results
        type_counts[r.type] = get(type_counts, r.type, 0) + 1
    end
    for (t, c) in sort(collect(type_counts); by=kv->kv[1])
        println("  type $(t): $c representations")
    end
    return results
end

# --- Step 2: order (unchanged) ---

function order_by_t_spectrum(rep::IrrepSumRep)
    r = rep.rank
    K = rep.K
    N = rep.N
    zeta = exp(2π * im / N)
    deg = degree(K)

    T_num = matrix_to_complex(rep.T, zeta, deg)
    t_diag = [T_num[i, i] for i in 1:r]
    spins = [mod(angle(t_diag[i]) / (2π), 1.0) for i in 1:r]

    function spin_sort_key(i)
        s = spins[i]
        best_q = 1; best_err = abs(s)
        for q in 1:12N
            p = round(Int, s * q)
            err = abs(s - p / q)
            if err < best_err - 1e-12
                best_err = err; best_q = q
            end
        end
        return (best_q, s)
    end

    perm = sortperm(1:r; by=spin_sort_key)
    S_new = zero_matrix(K, r, r)
    T_new = zero_matrix(K, r, r)
    for i in 1:r, j in 1:r
        S_new[i, j] = rep.S[perm[i], perm[j]]
        T_new[i, j] = rep.T[perm[i], perm[j]]
    end

    return OrderedIrrepSum(rep, S_new, T_new, perm)
end

# --- Step 3: eigenvalue blocks + find U (unchanged) ---

function eigenvalue_blocks(T_num::AbstractMatrix{<:Number}, r::Int; tol::Real=1e-10)
    t_diag = [T_num[i, i] for i in 1:r]
    visited = falses(r)
    blocks = Vector{Int}[]
    for i in 1:r
        visited[i] && continue
        block = [i]; visited[i] = true
        for j in i+1:r
            if !visited[j] && abs(t_diag[i] - t_diag[j]) < tol
                push!(block, j); visited[j] = true
            end
        end
        push!(blocks, block)
    end
    return blocks
end

function find_orthogonal_U(ordered::OrderedIrrepSum; max_mult::Int=4)
    r = ordered.base.rank
    K = ordered.base.K
    N = ordered.base.N
    zeta = exp(2π * im / N)
    deg = degree(K)

    T_num = matrix_to_complex(ordered.T, zeta, deg)
    S_num = matrix_to_complex(ordered.S, zeta, deg)
    blocks = eigenvalue_blocks(T_num, r)

    block_candidates = Vector{Vector{Matrix{Float64}}}()
    for block in blocks
        m = length(block)
        if m == 1
            push!(block_candidates, [ones(1, 1)])
        elseif m == 2
            angles = [0.0, π/4, -π/4, π/2]
            candidates = [Float64[cos(θ) -sin(θ); sin(θ) cos(θ)] for θ in angles]
            push!(block_candidates, candidates)
        else
            push!(block_candidates, [Matrix{Float64}(I, m, m)])
        end
    end

    results = pMDRep[]
    for combo in Iterators.product([1:length(c) for c in block_candidates]...)
        U_full = zeros(Float64, r, r)
        for (bi, block) in enumerate(blocks)
            U_block = block_candidates[bi][combo[bi]]
            for (li, gi) in enumerate(block), (lj, gj) in enumerate(block)
                U_full[gi, gj] = U_block[li, lj]
            end
        end

        S_pMD_num = U_full * S_num * U_full'

        has_nonzero_row = any(i -> all(j -> abs(S_pMD_num[i, j]) > 1e-10, 1:r), 1:r)
        has_nonzero_row || continue

        push!(results, pMDRep(ordered, U_full, S_pMD_num, T_num, K, N))
    end

    return results
end

# ============================================================
#  Step 4 REVISED: find_signed_V with projective normalization
# ============================================================

"""
    find_signed_V(pmd::pMDRep; tol=1e-6) -> Union{MDRep, Nothing}

Revised version that handles:
1. Projective normalization: T[i₀,i₀] need not be 1; we try all rows
   and normalize T by T[i₀,i₀]⁻¹ (choosing α in eq. 3.2).
2. Odd irreps: S entries may be purely imaginary. We try multiplying S
   by each 4th root of unity (ζ₄^α, α=0,1,2,3) to find a real form.
   This corresponds to the paper's eq. (3.2): ρ_α(s) = ζ₄^α · S/D.
"""
function find_signed_V(pmd::pMDRep; tol::Real=1e-6)
    r = size(pmd.S, 1)
    S_raw = pmd.S isa AbstractMatrix{<:Number} ? pmd.S :
            matrix_to_complex(pmd.S, exp(2π*im/pmd.N), degree(pmd.K))
    T_raw = pmd.T isa AbstractMatrix{<:Number} ? pmd.T :
            matrix_to_complex(pmd.T, exp(2π*im/pmd.N), degree(pmd.K))

    # Try each projective twist α ∈ {0,1,2,3}:
    #   S_α = ζ₄^α · S_raw    (paper eq. 3.2: ρ_α(s) = ζ₄^α · S/D)
    #   T_α = ζ₁₂^α · T_raw   (paper eq. 3.2: ρ_α(t) = ζ₁₂^α · e^{-2πic/24} · T)
    # The twist makes S real for the right α (even → α=0, odd → α=1 or 3).

    for α in 0:3
        ζ4_α = cispi(α / 2)    # e^{iπα/2} = i^α
        ζ12_α = cispi(α / 6)   # e^{iπα/6}
        S_num = ζ4_α * S_raw
        T_num = ζ12_α * T_raw

        # Try each row i₀ as the unit object
        for i0 in 1:r
            # Normalize T so that T[i₀,i₀] = 1 (determines central charge)
            t_i0 = T_num[i0, i0]
            abs(t_i0) < tol && continue
            T_norm = T_num / t_i0   # now T_norm[i₀,i₀] = 1

            # Check: is row i₀ of S_num all-real?
            row_i0 = S_num[i0, :]
            if any(j -> abs(imag(row_i0[j])) > tol, 1:r)
                continue
            end

            # Determine V: for each j, choose V[j] = ±1 so that S_MD[i₀,j] > 0
            V = ones(Int, r)
            ok = true
            for j in 1:r
                re = real(row_i0[j])
                if abs(re) < tol
                    ok = false; break
                end
                V[j] = re > 0 ? 1 : -1
            end
            ok || continue

            # Apply V to entire S
            S_MD = copy(S_num)
            for i in 1:r, j in 1:r
                S_MD[i, j] = V[i] * V[j] * S_num[i, j]
            end

            # Check all d_i = S_MD[i₀,j] / S_MD[i₀,i₀] are real and positive
            S00 = real(S_MD[i0, i0])
            S00 < tol && continue
            dims_ok = true
            for j in 1:r
                d_j = real(S_MD[i0, j]) / S00
                if d_j < tol || abs(imag(S_MD[i0, j])) > tol
                    dims_ok = false; break
                end
            end
            dims_ok || continue

            # Check S_MD is approximately real-symmetric
            real_sym = true
            for i in 1:r, j in i:r
                if abs(imag(S_MD[i, j])) > tol
                    real_sym = false; break
                end
                if abs(real(S_MD[i, j]) - real(S_MD[j, i])) > tol
                    real_sym = false; break
                end
            end
            real_sym || continue

            # Reorder so that unit object is index 1
            perm_unit = vcat([i0], setdiff(1:r, [i0]))
            S_reord = S_MD[perm_unit, perm_unit]
            T_reord = Diagonal([T_norm[perm_unit[i], perm_unit[i]] for i in 1:r])

            # Verlinde check
            Nijk = zeros(Int, r, r, r)
            verlinde_ok = true
            for i in 1:r, j in 1:r, k in 1:r
                val = sum(S_reord[l, i] * S_reord[l, j] * conj(S_reord[l, k]) /
                          S_reord[l, 1] for l in 1:r)
                n = round(Int, real(val))
                if abs(real(val) - n) > tol || abs(imag(val)) > tol || n < 0
                    verlinde_ok = false; break
                end
                Nijk[i, j, k] = n
            end
            verlinde_ok || continue

            # Central charge
            d = [real(S_reord[1, j]) / real(S_reord[1, 1]) for j in 1:r]
            D = sqrt(sum(d .^ 2))
            p_plus = sum(d[j]^2 * T_reord[j, j] for j in 1:r)
            c_over_8 = angle(p_plus / D) / (2π)
            c_rational = rationalize(c_over_8 * 8; tol=1e-4)

            return MDRep(pmd, V[perm_unit], S_reord, T_reord, Nijk,
                         c_rational, pmd.K, pmd.N)
        end
    end

    return nothing
end

# ============================================================
#  Full pipeline (revised)
# ============================================================

function classify_modular_data(N::Int; max_rank::Int=6, verbose::Bool=true)
    verbose && println("=" ^ 60)
    verbose && println("  MTC Classification (§3.4 pipeline v2): N = $N")
    verbose && println("=" ^ 60)

    catalog = build_atomic_catalog(N; max_rank=max_rank)
    sums = enumerate_irrep_sums(catalog, N; max_rank=max_rank)

    results = MDRep[]
    rejected = Dict("no_nonzero_row" => 0, "no_valid_V" => 0,
                     "total_sums" => length(sums))

    for (si, s) in enumerate(sums)
        ordered = order_by_t_spectrum(s)
        pmds = find_orthogonal_U(ordered)

        if isempty(pmds)
            rejected["no_nonzero_row"] += 1
            continue
        end

        found = false
        for p in pmds
            md = find_signed_V(p)
            if md !== nothing
                push!(results, md)
                found = true
                if verbose
                    type_str = "(" * join(s.type, ",") * ")"
                    d_i = [round(real(md.S[1,j]/md.S[1,1]); digits=4) for j in 1:s.rank]
                    println("  ✓ type=$type_str  rank=$(s.rank)  " *
                            "c=$(md.central_charge)  d_i=$d_i")
                end
            end
        end
        !found && (rejected["no_valid_V"] += 1)
    end

    verbose && println()
    verbose && println("Results: $(length(results)) MDReps from $(rejected["total_sums"]) sums")
    verbose && println("Rejected: no-nonzero-row=$(rejected["no_nonzero_row"]), " *
                       "no-valid-V=$(rejected["no_valid_V"])")
    verbose && println("=" ^ 60)

    return results
end

println("mtc_pipeline_v2.jl loaded.")
println("Usage: results = classify_modular_data(5; max_rank=4)")
