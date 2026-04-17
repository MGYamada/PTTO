"""
Phase 3 integration: multi-prime CRT reconstruction of SU(2)_4.

Runs the Phase 2 pipeline at multiple primes, groups the resulting
MTC candidates by fusion-tensor invariants, and reconstructs the
S-matrix entries as elements of Z[√3].

Expected outcome (matching Python's M6 result):
    2√3 · S' = M  with M ∈ Z[√3]^{5×5} having entries in {0, ±1, ±2, ±√3}

Run with:
    julia --project=. scripts/su24_crt.jl
"""

using Oscar
using ACMG

println("=== Phase 3 integration: multi-prime CRT for SU(2)_4 ===\n")

# Step 1: Build atomic catalog
println("Building atomic catalog for N = 24...")
catalog = ACMG.build_atomic_catalog(24; max_rank = 5, verbose = false)
println("  Catalog: $(length(catalog)) atomic irreps")

# Step 2: Find a known-good (ρ_3, ρ_2) pair (we know from Phase 2 script
# that e.g. d3_idx=81, d2_idx=49 works)
d3_idx = 81
d2_idx = 49
stratum = ACMG.Stratum(Dict(d3_idx => 1, d2_idx => 1), 5)
println("\nUsing stratum: catalog[$d3_idx] ⊕ catalog[$d2_idx]")
println("  ρ_3 = $(catalog[d3_idx])")
println("  ρ_2 = $(catalog[d2_idx])")

# Step 3: Run find_mtcs_at_prime at multiple primes
# Good primes with 24 | p-1 and 3 is a QR (so √3 ∈ F_p):
# p = 73, 97, 193, 241 all satisfy both
test_primes = [73, 97, 193, 241]
fresh_primes = [313, 337, 409]  # for cross-validation

println("\nRunning Phase 2 sweep at primes $(test_primes)...")
results_by_prime = Dict{Int, Vector{ACMG.MTCCandidate}}()
for p in test_primes
    cands = ACMG.find_mtcs_at_prime(catalog, stratum, p; verlinde_threshold = 3)
    results_by_prime[p] = cands
    println("  p=$p: $(length(cands)) candidates")
end

# Step 4: Group by fusion tensor
println("\nGrouping candidates by fusion tensor...")
groups = ACMG.group_mtcs_by_fusion(results_by_prime)
println("  Found $(length(groups)) distinct fusion-tensor groups")

for (gi, g) in enumerate(groups)
    n_primes = length(g)
    rep = first(values(g))
    d_signed = [ACMG.signed_Fp(x, rep.p) for x in rep.d]
    println("  Group $gi: present at $n_primes primes, unit=$(rep.unit_index), d=$d_signed")
end

# Step 5: For each complete group (one candidate per prime), reconstruct S in Z[√3]
println("\n=== Reconstructing S-matrix in Z[√3] for each group ===")

for (gi, group) in enumerate(groups)
    if length(group) != length(test_primes)
        println("\nGroup $gi: incomplete (only $(length(group))/$(length(test_primes)) primes), skipping")
        continue
    end

    println("\nGroup $gi — complete across all $(length(test_primes)) primes")
    rep = first(values(group))
    println("  Representative: $rep")

    # Reconstruct 2√3 · S (so that entries are integers in Z[√3])
    local recon
    try
        recon = ACMG.reconstruct_S_matrix(group; scale_d = 3, bound = 5)
    catch e
        println("  ✗ Reconstruction failed: $e")
        continue
    end

    println("\n  Reconstructed matrix 2√3 · S (entries in Z[√3]):")
    println(ACMG.describe_matrix(recon, 3))

    # Sanity check at used primes
    println("\n  Verifying at used primes (must pass):")
    for p in test_primes
        ok = ACMG.verify_reconstruction(recon, group[p], 3; scale = 2)
        println("    p=$p: $(ok ? "✓" : "✗")")
    end

    # Cross-validate at fresh primes
    println("\n  Cross-validation at FRESH primes (unused in reconstruction):")
    for p in fresh_primes
        cands = ACMG.find_mtcs_at_prime(catalog, stratum, p; verlinde_threshold = 3)
        # Find candidate with matching fusion tensor
        matching = nothing
        for c in cands
            if c.N == rep.N && c.unit_index == rep.unit_index
                matching = c
                break
            end
        end
        if matching === nothing
            println("    p=$p: no matching candidate found ⚠")
            continue
        end
        ok = ACMG.verify_reconstruction(recon, matching, 3; scale = 2)
        println("    p=$p: $(ok ? "✓" : "✗")")
    end
end

println("\n=== Phase 3 integration complete ===")
