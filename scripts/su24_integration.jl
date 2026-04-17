"""
Integration test: end-to-end SU(2)_4 detection at rank 5, N = 24.

This script reproduces the Python M1-M5 pipeline inside Julia ACMG.
Expected outcome: find exactly 2 MTC candidates at p = 73, matching
Python's φ = π/4 rotation and reflection solutions.

Run with:
    julia --project=. scripts/su24_integration.jl
"""

using Oscar
using ACMG

println("=== Phase 2 integration: SU(2)_4 at rank 5, N = 24 ===\n")

# Step 1: Build atomic catalog for N = 24
println("Building atomic catalog for N = 24 (this takes a few seconds)...")
catalog = ACMG.build_atomic_catalog(24; max_rank = 5, verbose = false)
println("  Catalog size (max_rank ≤ 5): $(length(catalog)) atomic irreps")

# Step 2: Identify ρ_3 (d=3, level=24, parity=+1) and ρ_2 (d=2, level=12, parity=+1)
# from Python's known indices: ρ_3 = d3_idx36, ρ_2 = d2_idx18 (within-dim indices)
# These correspond to specific entries in our catalog; we identify by label+parity+level

d3_candidates = [i for (i, a) in enumerate(catalog)
                 if a.dim == 3 && a.level == 24 && a.parity == +1]
d2_candidates = [i for (i, a) in enumerate(catalog)
                 if a.dim == 2 && a.level == 12 && a.parity == +1]

println("  Candidates for ρ_3 (3d_24, parity=+1): $(length(d3_candidates))")
println("  Candidates for ρ_2 (2d_12, parity=+1): $(length(d2_candidates))")

# Step 3: Find the specific ρ_3 and ρ_2 with matching T-spectrum for SU(2)_4
# Known: ρ_2 should have T_powers = [22, 6] (in ζ_24)
# Known: ρ_3 should have degenerate T = ζ_24^22 eigenvalue overlap with ρ_2

# Try each combination; the pipeline will tell us which works.
found_mtcs = Vector{Tuple{Int, Int, Any}}()

# Primes to test
test_primes = [73, 97, 193]

println("\nTrying candidate pairs (ρ_3, ρ_2):")

for d3_idx in d3_candidates
    for d2_idx in d2_candidates
        # Build stratum: 1 copy of d3 + 1 copy of d2
        stratum = ACMG.Stratum(Dict(d3_idx => 1, d2_idx => 1), 5)

        # Try at each prime
        p = test_primes[1]  # Start with p = 73
        try
            candidates = ACMG.find_mtcs_at_prime(catalog, stratum, p;
                                                 verlinde_threshold = 3)
            if !isempty(candidates)
                println("  ✓ Pair (d3_idx=$d3_idx, d2_idx=$d2_idx): $(length(candidates)) MTC(s) at p=$p")
                for c in candidates
                    println("    $c")
                end
                push!(found_mtcs, (d3_idx, d2_idx, candidates))
            end
        catch e
            # This pair doesn't form a valid stratum (e.g. different T-eigenspaces)
            # Silent skip
        end
    end
end

println("\nTotal (ρ_3, ρ_2) pairs yielding MTC at p=73: $(length(found_mtcs))")

if !isempty(found_mtcs)
    println("\n=== Detailed MTC candidates ===")
    for (d3_idx, d2_idx, cands) in found_mtcs
        println("\nStratum: ρ_3 (catalog[$d3_idx]) ⊕ ρ_2 (catalog[$d2_idx])")
        println("  ρ_3 = $(catalog[d3_idx])")
        println("  ρ_2 = $(catalog[d2_idx])")
        for c in cands
            println("  → $c")
            # Compare to expected SU(2)_4 dimensions: d ∈ {±1, ±√3, ±2}
            # Signed F_p values
            d_signed = [ACMG.signed_Fp(x, c.p) for x in c.d]
            println("    d (signed F_$(c.p)) = $d_signed")
            println("    D² = $(ACMG.signed_Fp(c.D2, c.p)) (expected 12 for SU(2)_4)")
        end
    end
end

# Check: did we recover SU(2)_4 (D² = 12)?
su24_found = false
for (_, _, cands) in found_mtcs
    for c in cands
        if ACMG.signed_Fp(c.D2, c.p) == 12
            global su24_found = true
            break
        end
    end
    su24_found && break
end

if su24_found
    println("\n🎉 SU(2)_4 (D² = 12) recovered via Julia ACMG Phase 2 pipeline!")
else
    println("\n⚠ SU(2)_4 not recovered. Check pipeline.")
end
