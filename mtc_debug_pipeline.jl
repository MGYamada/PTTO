# mtc_debug_pipeline.jl
# Debug helper: trace what happens at each step for every IrrepSum
# Run after: include("mtc_classifier.jl"); include("mtc_pipeline.jl")

"""
    debug_classify(N; max_rank=4)

Like classify_modular_data but with verbose tracing at each step.
"""
function debug_classify(N::Int; max_rank::Int=4)
    catalog = build_atomic_catalog(N; max_rank=max_rank)
    sums = enumerate_irrep_sums(catalog, N; max_rank=max_rank)

    K = catalog[1].K
    zeta = exp(2π * im / N)
    deg = degree(K)

    for (si, s) in enumerate(sums)
        type_str = "(" * join(s.type, ",") * ")"
        comp_str = join(["$(catalog[α].label)×$m" for (α, m) in s.components], " ⊕ ")
        println("\n", "─" ^ 60)
        println("IrrepSum [$si]: type=$type_str rank=$(s.rank)  [$comp_str]")

        # Step 2: order
        ordered = order_by_t_spectrum(s)
        T_num = matrix_to_complex(ordered.T, zeta, deg)
        S_num = matrix_to_complex(ordered.S, zeta, deg)
        t_diag = [T_num[i, i] for i in 1:s.rank]
        spins = [mod(angle(t_diag[i]) / (2π), 1.0) for i in 1:s.rank]
        println("  T spins (sorted): ", [round(s; digits=4) for s in spins])

        # Check parity of direct sum
        par = compute_parity(S_num)
        println("  Parity of ⊕: $(par == 1 ? "even" : par == -1 ? "odd" : "mixed(0)")")

        # Step 3: eigenvalue blocks
        blocks = eigenvalue_blocks(T_num, s.rank)
        block_sizes = [length(b) for b in blocks]
        println("  Eigenvalue blocks: $block_sizes ($(length(blocks)) blocks)")

        # Check: which rows of S have all nonzero entries BEFORE U?
        nonzero_rows_before = Int[]
        for i in 1:s.rank
            if all(j -> abs(S_num[i, j]) > 1e-10, 1:s.rank)
                push!(nonzero_rows_before, i)
            end
        end
        println("  S rows all-nonzero BEFORE U: $nonzero_rows_before")

        # Show S matrix
        println("  S (numeric, ordered):")
        for i in 1:s.rank
            entries = [abs(S_num[i, j]) < 1e-10 ? "  ·  " :
                       string(round(S_num[i, j]; digits=3)) for j in 1:s.rank]
            println("    ", join(lpad.(entries, 12)))
        end

        # Try each U candidate
        pmds = find_orthogonal_U(ordered)
        println("  U candidates passing nonzero-row filter: $(length(pmds))")

        if isempty(pmds)
            println("  → REJECTED: no U gives a nonzero row (block-diagonal zero issue)")
            # Show why: for each angle combo, show the resulting S[i,:]
            if all(bs == 1 for bs in block_sizes)
                println("    All blocks are 1×1 → U=I. The block-diagonal S itself has no")
                println("    all-nonzero row. This is expected for nontrivial direct sums")
                println("    — the paper's U is not just ±1 in the mult-1 blocks; the")
                println("    off-diagonal zeros require mixing from mult≥2 blocks.")
            end
            continue
        end

        for (pi, p) in enumerate(pmds)
            md = find_signed_V(p)
            if md !== nothing
                d_i = [round(real(md.S[1, j] / md.S[1, 1]); digits=4) for j in 1:s.rank]
                println("  → MDRep found! d_i = $d_i, c = $(md.central_charge)")
            else
                # Why did V fail?
                S_p = p.S
                println("  pMD [$pi]: V search failed. Checking details...")
                for i in 1:s.rank
                    t_i = p.T[i, i]
                    is_12th = abs(t_i^12 - 1) < 1e-6
                    row_real = all(j -> abs(imag(S_p[i, j])) < 1e-6, 1:s.rank)
                    row_pos = row_real && all(j -> real(S_p[i, j]) > 1e-6, 1:s.rank)
                    println("    row $i: T^12≈1? $is_12th  all-real? $row_real  all-pos? $row_pos  " *
                            "vals=$(round.(real.(S_p[i, :]); digits=3))")
                end
            end
        end
    end
end

println("debug helper loaded. Run: debug_classify(5; max_rank=4)")
