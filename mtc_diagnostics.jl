# mtc_diagnostics.jl
# Diagnostic helpers that run alongside the existing mtc_classifier.jl to
# verify our understanding against paper's predictions, before we touch
# the main pipeline.
#
# This file assumes mtc_classifier.jl is already loaded (so that
# inspect_atomic_irreps, matrix_to_complex, etc. are in scope).

include("mtc_types.jl")

"""
    diagnose_atomic_irreps(N; verbose=true)

For each atomic irrep of SL₂(ℤ/NZ) produced by the current pipeline,
compute and report:
  - parity (ρ(s²) = ±I ?)
  - whether S[1,i]/S[1,1] is real (Theorem 2.1(3) compliance of *atomic*)
  - numeric level (order of T)
  - non-degenerate block and its spectrum

This is a read-only probe — it does NOT modify any existing code path.
Its job is to make the implicit structure visible, so we can see whether
our claim "2_5^1 is odd; its raw S is purely imaginary" matches reality
before we build the normalization machinery.

Returns the list of AtomicIrrep wrappers (one per entry in the catalog),
for further programmatic use.
"""
function diagnose_atomic_irreps(N::Int; verbose::Bool=true)
    # Build the catalog via the existing inspector, but suppress its
    # voluminous printing by redirecting stdout temporarily.
    catalog = redirect_stdout(devnull) do
        inspect_atomic_irreps(N; max_rank=20)
    end

    K, _ = cyclotomic_field(N)
    zeta = exp(2π * im / N)
    deg  = degree(K)

    atomics = AtomicIrrep[]

    verbose && println("=" ^ 70)
    verbose && println("Atomic irreps of SL₂(ℤ/$(N)ℤ) — parity & Thm 2.1(3) check")
    verbose && println("=" ^ 70)

    for (α, entry) in enumerate(catalog)
        S_num = matrix_to_complex(entry.S, zeta, deg)
        T_num = matrix_to_complex(entry.T, zeta, deg)

        par    = compute_parity(S_num)
        level  = _numeric_order(T_num, entry.rank, N)
        dim_ok = has_dim_real_numeric(S_num)
        perron = fp_ratio_nonnegative(S_num)
        indices, block = non_degenerate_block(S_num, T_num)

        par_str = par == +1 ? "even" :
                  par == -1 ? "odd " : "mixd"

        atomic = AtomicIrrep(
            entry.rank, level, entry.name,
            entry.S, entry.T, par, K, N,
        )
        push!(atomics, atomic)

        if verbose
            println()
            println("[α=$α] $(entry.name)  rank=$(entry.rank)  level=$level  parity=$par_str")
            println("      S[1,1] = $(round(S_num[1,1]; digits=4))  " *
                    "(imag=$(round(imag(S_num[1,1]); sigdigits=3)))")
            println("      T diagonal = ", [round(T_num[i,i]; digits=4) for i in 1:entry.rank])
            println("      d_i = S[1,i]/S[1,1] all-real? $dim_ok    FP-nonnegative? $perron")
            println("      non-degenerate indices: $indices   (block size $(size(block)))")
            if !isempty(indices)
                block_str = "      non-deg S-block =\n"
                for i in 1:size(block, 1)
                    block_str *= "        " *
                        join([string(round(block[i, j]; digits=4)) for j in 1:size(block, 2)], "  ") *
                        "\n"
                end
                print(block_str)
            end
        end
    end
    verbose && println("=" ^ 70)
    return atomics
end

"""
    check_existing_sign_fix(N; catalog_index=1)

Take a specific atomic irrep and compare:
  (a) the naive `fix_signs_exact` output, vs
  (b) what Theorem 3.4 / Remark 3.8 would predict.

For an *atomic* irrep there is no basis freedom U (the T is not always
non-degenerate but the irrep is already in its canonical form); only V
survives. So this test isolates whether `fix_signs_exact` is doing the
Remark 3.8 step correctly.

Intent: spot-check that, e.g., for 2_5^1 (the Fibonacci atomic), the
sign-fix should *not* succeed without also handling the odd-parity
(s_5^1)^{-1} prefactor — which the current code ignores.
"""
function check_existing_sign_fix(N::Int; catalog_index::Int=1, verbose::Bool=true)
    catalog = redirect_stdout(devnull) do
        inspect_atomic_irreps(N; max_rank=20)
    end
    catalog_index > length(catalog) && error(
        "Only $(length(catalog)) atomic irreps for N=$N; asked for [$catalog_index].")

    entry = catalog[catalog_index]
    K, _ = cyclotomic_field(N)
    zeta = exp(2π * im / N)
    deg  = degree(K)

    S_num_before = matrix_to_complex(entry.S, zeta, deg)
    S_fixed, signs = fix_signs_exact(entry.S, K, N, entry.rank)
    S_num_after  = matrix_to_complex(S_fixed, zeta, deg)

    par = compute_parity(S_num_before)
    dim_ok_before = has_dim_real_numeric(S_num_before)
    dim_ok_after  = has_dim_real_numeric(S_num_after)

    if verbose
        println("=" ^ 70)
        println("sign-fix diagnostic on $(entry.name)  (atomic [$catalog_index], N=$N)")
        println("=" ^ 70)
        println("parity = $(par == +1 ? "even" : par == -1 ? "odd" : "mixed")")
        println("signs chosen: $signs")
        println("d_i all-real before fix? $dim_ok_before")
        println("d_i all-real after  fix? $dim_ok_after")
        println()
        println("S before:")
        _print_complex_matrix(S_num_before)
        println("S after fix_signs_exact:")
        _print_complex_matrix(S_num_after)
        println()
        if par == -1 && !dim_ok_after
            println("DIAGNOSIS: this irrep is odd and fix_signs_exact was NOT")
            println("sufficient to make d_i real. Remark 3.8 V alone is not enough;")
            println("the odd-rep scalar prefactor (s_n^m)⁻¹ must be absorbed too.")
            println("This is the 2_5^1 / Fibonacci issue identified in the summary.")
        elseif par == -1 && dim_ok_after
            println("NOTE: odd irrep passed — fix_signs_exact happens to work here.")
        end
        println("=" ^ 70)
    end

    return (entry=entry, S_before=S_num_before, S_after=S_num_after,
            signs=signs, parity=par,
            dim_ok_before=dim_ok_before, dim_ok_after=dim_ok_after)
end

function _print_complex_matrix(M::AbstractMatrix{<:Number})
    for i in 1:size(M, 1)
        print("  ")
        for j in 1:size(M, 2)
            x = M[i, j]
            re = round(real(x); digits=4)
            im_ = round(imag(x); digits=4)
            if abs(im_) < 1e-10
                print(lpad(string(re), 12))
            else
                print(lpad(string(re) * (im_ ≥ 0 ? "+" : "") * string(im_) * "i", 20))
            end
            print("  ")
        end
        println()
    end
end

println("mtc_diagnostics.jl loaded.")
println("Try: diagnose_atomic_irreps(5)")
println("     check_existing_sign_fix(5; catalog_index=1)")
