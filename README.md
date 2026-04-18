# ACMG — Arithmetic Condensed Matter Geometry

Modular Tensor Category (MTC) classification via cyclotomic representation
variety enumeration with multi-prime F_p validation.

## Overview

ACMG classifies MTCs by **fixing the conductor** N = max{cond(S), cond(T)}
and exploiting the following key structural insight:

> The moduli space of SL(2, ℤ/N)-modular data of rank r is a
> **stratified scheme** over Spec ℤ[1/2N]. Each stratum is indexed by the
> isomorphism type (m_λ) of the underlying SL(2, ℤ/N)-representation
> ρ = ⊕_λ ρ_λ^{m_λ}. Within each stratum, the moduli is parametrised by
> a block-U group O(n_θ) on each T-eigenspace V_θ (with dim n_θ = dim V_θ).
> MTC integrality (Verlinde non-negativity) then cuts out a 0-dimensional
> locus inside this continuous family.

This gives a pipeline:

```
Atomic SL(2, ℤ/N)-irreps catalog     (via Oscar + GAP/SL2Reps)
        ↓
Stratum enumeration (m_λ with Σ m_λ d_λ = r)
        ↓
T-spectrum filter (retain (m_λ) matching target T-eigenvalues)
        ↓
Block-U parametrisation on each T-eigenspace V_θ
        ↓
F_p multi-prime Verlinde sweep (fast, exact)
        ↓
CRT reconstruction to Z[ζ_N]
        ↓
Classified MTCs
```

The pipeline has been validated end-to-end on SU(2)_4 at rank 5, N = 24.

## Current status (v0.3 prototype)

### Phase 0: Atomic catalog — ✅ complete

- `FpArith.jl`: Tonelli-Shanks sqrt, Euler criterion, primitive roots,
  roots of unity, matrix ops mod p.
- `Types.jl`: `ModularDatumFp`, `FusionRule` with ring-agnostic axiom
  validation.
- `ModularData.jl`: (S, T) axiom checking, α (= e^{-2πi c/8}) extraction
  from (ST)³ = α · S², charge conjugation from S².
- `FusionExtract.jl`: F_p Verlinde coefficients, Z lift with
  non-negativity check.
- `SL2Reps.jl`: Oscar + GAP/SL2Reps integration, `build_atomic_catalog(N)`
  producing the irreducible SL(2, ℤ/N)-representations.

### Phase 1: Stratum enumeration — ✅ complete

- `StratumEnum.jl`: combinatorial enumeration of {m_λ} with
  Σ m_λ · dim(λ) = r.
- Note: `require_unit_summand=false` is the correct default. In non-
  abelian MTCs the unit object sits as a basis vector *inside* a larger
  irrep block (wherever T-eigenvalue 1 occurs), not as a separate 1d_1
  summand.

### Phase 2: Block-U parametrisation — ✅ complete

- `BlockU.jl`: general O(n) Cayley parametrisation, T-eigenspace
  decomposition, Verlinde integrality check with unit identification.
- `find_mtcs_at_prime(catalog, stratum, p)`: full per-prime driver.
- Single degenerate eigenspace of dimension n_θ ≤ 3 is currently
  supported (feasibility limit of naive Cayley sweep; see table below).
- Validated end-to-end on SU(2)_4 at rank 5, N = 24 — all matching
  strata recover exactly 2 Galois-conjugate MTC candidates per prime.

### Phase 3: CRT reconstruction — ✅ complete

- `CRT.jl`: multi-prime CRT, rational reconstruction, Galois-aware
  grouping, Z[√d] entry-wise reconstruction.
- `group_mtcs_galois_aware(results, anchor_prime)`: aligns candidates
  across primes into the same Galois sector via 2-prime trial
  reconstruction — a fusion-tensor match alone is insufficient because
  each prime returns Galois-conjugate pairs.
- `compute_sqrt3_cyclotomic_mod_p(p)`: uses ζ₂₄² + ζ₂₄²² for a
  Galois-consistent choice of √3 across primes (not "smaller
  representative", which is NOT Galois-invariant).
- Validated on SU(2)_4 N = 24: reconstructs 2√3·S ∈ Z[√3]^{5×5} from
  4 primes {73, 97, 193, 241}, cross-validated at 3 fresh primes
  {313, 337, 409}. Both Galois orbit members (Group 1, Group 2) recover
  7/7 ✓ each.

### Test statistics

**334 tests passing** (runtests.jl):
- 63 FpArith primitives
- 34 Fibonacci at p = 41 (N = 5, rank 2)
- 30 Ising at p = 17 (N = 16, rank 3)
- 9 SL2Reps helpers + atomic catalog construction (N ∈ {5, 8, 16})
- 15 StratumEnum (synthetic rank 1–5 cases)
- 72 BlockU pure F_p helpers (incl. Fibonacci and Ising Verlinde check)
- 40 CRT primitives (crt, rational_reconstruct, Z[√d] reconstruction,
  Galois-aware grouping)
- 71 General O(n) Cayley (inverse_mod_p, enumerate_so_n_Fp,
  enumerate_o_n_Fp, apply_block_U)


## Design notes

### Central conceptual shift from v5

The v5 pipeline fixed a **prime p** and solved modular-data equations
over F_p directly. The v6 / ACMG pipeline fixes the **conductor N** and
works primarily with SL(2, ℤ/N)-representation theory over Z[ζ_N], using
F_p only for fast computation and final validation. The reasons:

1. **Finiteness** of the search space. Fixing N restricts modular data
   entries to Z[ζ_N] ∩ {|z| ≤ D} — a finite set. The p-fix approach
   searches F_p instead, also finite but without the algebraic structure.

2. **Natural alignment with NRW(W)** (2203.14829, 2308.09670). Their
   classification stratifies by SL(2, ℤ/N) irrep content — we follow the
   same stratification.

3. **Conductor-first beats rank-first**. The NRW paper enumerates rank
   first then searches for fields, which forces unnecessary field
   extensions when the MTC actually lives in a smaller cyclotomic field.

### Modular data moduli: what the geometry actually is

Consider the moduli space

    M_{r, N}^{MTC} ⊂ M_{r, N}^{SL_2(ℤ/N)-rep}

of rank r modular data of conductor dividing N. The right-hand side is
the representation variety (up to equivalence), with the additional
conditions S symmetric and S² = charge-conjugation.

**Stratification.** Each isomorphism class ρ = ⊕_λ ρ_λ^{m_λ} defines a
locally closed stratum. Within a stratum, the choice of basis for V
(not changing the SL(2, ℤ/N)-module structure) is an element of
Aut_G(ρ) = ∏_λ GL(m_λ). For all m_λ = 1, this is just a torus.

**Continuous moduli are hidden in T-eigenspace overlaps.** The surprise
(discovered in this iteration) is that even when all m_λ = 1 — so
Aut_G(ρ) is just a torus acting trivially on S — there can still be a
continuous moduli if **different irreducibles share a T-eigenvalue θ**.
Specifically, the T-eigenspace

    V_θ = ⊕_λ V_θ^{(ρ_λ)}

can be 2-dimensional, and rotating its basis (an O(2) action) gives a
new S-matrix while preserving T and the SL(2, ℤ/N)-module structure.

**Block-U, algebraic definition.** The full group acting is

    𝒰(ρ) = {U ∈ GL(V) : [T, U] = 0,
                          U|_{V_θ} regular semisimple,
                          U Uᵀ ∈ End_G(ρ)}

(the last being projective orthogonality). Decomposing U = ⊕_θ U_θ on
T-eigenspaces, the continuous moduli has dimension

    parameter_dim = Σ_θ C(n_θ, 2)   where n_θ = dim V_θ

For SU(2)_4 at N = 24: T-spectrum (with multiplicities) has a single
degenerate eigenvalue (multiplicity 2), so parameter_dim = 1. This
matches the Python M4 sweep finding φ = π/4.

### Verlinde integrality: from continuous to discrete

On the continuous moduli the Verlinde fusion coefficients
N_{ij}^k(φ) are algebraic functions of the block-U parameters. The MTC
locus is the subset of parameter space where all N_{ij}^k are
non-negative integers — generically a **0-dimensional algebraic set**.

For SU(2)_4 this is φ = π/4 (a rotation / reflection pair). Not an
accident: the relevant 2×2 block of S on V_1 is rank-1 with all entries
equal (by Kac–Peterson), and its eigenvector is the Hadamard direction,
i.e., φ = π/4. **This is derivable a priori, without any numerical
sweep.**

### The étale picture (partial)

Strata where all n_θ = 1 are 0-dimensional: parameter_dim = 0. The
scheme 𝒰 reduces to ∏_θ μ_2 = μ_2^{#θ} on these strata, which is
étale over Z[1/2]. Rank ≤ 4 over small N tends to fall in this case,
and indeed: empirically most "simple" MTCs (Fibonacci, Ising family)
have no continuous moduli.

When some n_θ ≥ 2, continuous moduli appear and MTC integrality
becomes a non-trivial cutting condition. Rank 5 N = 24 is the smallest
case with this phenomenon.

Precisely: **the final set of MTCs** (after Verlinde integrality cuts
out a 0-dimensional locus) is still étale over Z[1/2N], but the
intermediate universe of SL(2, ℤ/N)-modular data is not. This is the
correct version of "MTC moduli is étale" — the statement applies to
rational points after integrality, not to the ambient variety.

### Degenerate perturbation theory as a mental model

The mathematical content of block-U is exactly **degenerate
perturbation theory** in algebraic-geometric disguise:

- T = unperturbed Hamiltonian with degenerate spectrum.
- V_θ (n_θ ≥ 2) = degenerate subspace.
- S|_{V_θ} = perturbation lifting the degeneracy.
- O(n_θ) freedom = choice of 0-th order eigenbasis.
- Hadamard / φ = π/4 = the correct eigenbasis for SU(2)_4.

Explicitly, the tangent direction δT that resolves the degeneracy of T
at V_1 is δT ∝ σ_x (off-diagonal Pauli). Its eigenvectors are exactly
the Hadamard basis, recovering φ = π/4 analytically.

### Multi-prime F_p validation

Practical computation uses **F_p reduction at multiple good primes**
(N | p-1), even though the underlying moduli lives in Z[ζ_N].
Rationale:

1. Block-U parametrisation can be computed exactly in F_p via the
   algebraic circle {(u, v) : u² + v² = 1} ⊂ F_p² — only ~p points.
2. Verlinde integrality check in F_p is an integer comparison.
3. CRT across primes reconstructs the algebraic entries in Z[ζ_N].
4. Sanity check: verify at fresh primes not used for reconstruction.

This has been validated end-to-end on SU(2)_4 rank 5 N = 24:

- 10 primes {73, 97, 193, 241, 313, 337, 409, 433, 457, 577}.
- Each prime gives 2 MTC solutions (φ = π/4, rotation & reflection).
- 4-prime CRT reconstructs S = (1/2√3) · M with M ∈ Z[√3]^{5×5}
  entries in {0, ±1, ±2, ±√3}.
- Verified at 6 unused primes: all ✓.
- Quantum dims d = (√3, √3, 1, 1, 2), D² = 12.

Fast (~1 second), exact (no floating-point), parallel (primes
independent).

### Galois-aware grouping across primes

A subtle but important point in multi-prime CRT reconstruction:
**two primes picking up the same abstract MTC may return different
Galois-conjugate representatives** of the (S, T)-matrix pair.

Concretely, `find_mtcs_at_prime` returns multiple candidates per prime
(e.g. for SU(2)_4, 2 candidates — a Galois-conjugate pair that differ
by √3 ↔ −√3). Across primes, the "first" candidate and the "second"
candidate in the ordered output may switch roles arbitrarily, because
(a) the F_p enumeration order of O(n)(F_p) points is prime-dependent,
and (b) the choice of primitive root in `find_zeta_in_Fp(N, p)` varies
across primes.

**Fusion tensors are Galois-invariant** — so matching candidates across
primes by N-tensor alone lumps together Galois-conjugate pairs into a
single group. But CRT requires **same Galois sector at every prime**,
otherwise entries won't lift to a consistent Z[√d] element.

**Solution: `group_mtcs_galois_aware(results, anchor_prime)`.** For each
candidate at the anchor prime, seed a group. For each other prime, pick
the candidate whose S-matrix, together with the anchor's, is consistent
under 2-prime Z[√d] reconstruction (i.e. actually lifts with small
integer coefficients). This picks the "correct" Galois-conjugate at each
prime, yielding one Galois-coherent group per orbit member.

For SU(2)_4 this correctly splits the collection into 2 groups (one per
Galois orbit element), each lifting to a consistent Z[√3] matrix that
verifies at fresh primes.

This construction is not in NRWW; it arises specifically because our
pipeline computes in F_p rather than an abstract number field and must
cope with the Galois ambiguity of F_p realisations.

### What becomes obsolete

- **v5 Part B**: the "2×2 rotation on T-overlapping atomic pairs"
  framework is partially correct in its setup (shared-eigenvalue O(2))
  but v5 got parametrisation details wrong or simply tested conductors
  (N = 8, 12, 15, 20) where no genuine multi-component MTC exists.
- **p-fix as primary search**: p-fix is retained as an accelerator, not
  as the primary search space.
- **Homotopy continuation / damped Newton**: replaced by the block-U
  formulation, which is a direct algebraic construction.

## Prime selection for F_p validation

For conductor N, select primes p with N | p − 1:

| N | Example good primes |
|-----|-----|
|  5 | 11, 31, 41, 61, 71 |
|  8 | 17, 41, 73, 89, 97 |
| 12 | 13, 37, 61, 73, 97 |
| 24 | 73, 97, 193, 241, 313, 337, 409 |

A prime p is **admissible** for a given MTC candidate if:

- N | p − 1 (so ζ_N ∈ F_p),
- √3, √2, and in general √(D²-factors) are in F_p (quadratic residue
  condition),
- all N_{ij}^k lift to small non-negative integers.

Chebotarev density ⇒ for each MTC, a positive-density set of primes is
admissible. Running at several primes in parallel gives robust
confirmation.

## Repository layout

```
src/
  ACMG.jl                 — module root, exports
  FpArith.jl              — F_p primitives (Tonelli-Shanks, matmul, etc.)
  Types.jl                — ModularDatumFp, FusionRule
  ModularData.jl          — (S, T) axiom validation
  FusionExtract.jl        — F_p Verlinde → Z lift
  Dimensions.jl           — (stub) quantum dimension enumeration
  Enumerator.jl           — (stub) top-level driver
  SL2Reps.jl              — Oscar + GAP/SL2Reps catalog builder (Phase 0)
  StratumEnum.jl          — {m_λ} combinatorial enumeration (Phase 1)
  BlockU.jl               — general O(n) Cayley + single-prime driver (Phase 2)
  CRT.jl                  — multi-prime CRT + Galois-aware grouping (Phase 3)

test/
  runtests.jl
  test_fparith.jl
  test_fibonacci.jl
  test_ising.jl
  test_sl2reps.jl
  test_stratum_enum.jl
  test_blocku.jl          — pure F_p block-U helpers
  test_crt.jl             — CRT primitives + Z[√d] reconstruction
  test_block_u_general.jl — general O(n) Cayley parametrisation

scripts/
  su24_integration.jl     — Phase 2 end-to-end for SU(2)_4
  su24_crt.jl             — Phase 3 end-to-end (multi-prime CRT)
  diagnose_crt.jl         — per-prime candidate inspection (debugging aid)
```

## References

- Bruillard, Ng, Rowell, Wang, *Rank-finiteness for modular categories*
  (2016). Admissibility criteria used in filtering.
- Ng, Rowell, Wang, Wen, *Reconstruction of modular data from
  SL(2, ℤ) representations* arXiv:2203.14829. The irrep-sum + block-U
  framework.
- Ng, Rowell, Wen, *Classification of modular data up to rank 12*
  arXiv:2308.09670. Ground-truth catalog (NsdGOL files) used for
  cross-validation.

## Next steps

**Near-term:**

1. **Multi-variant validation**: reproduce all 16 variants of
   NsdGOL[2] ∪ NsdGOL[3] (rank 5 N = 24) through the pipeline.
   Each (ρ_3, ρ_2) pair from Phase 2 gives a Galois orbit; verify
   all orbit members are covered.
2. **Phase 1.5**: `enumerate_strata_by_T(catalog, target_spins)` to
   pre-filter strata by T-spectrum multiset (reduces ~25k strata for
   rank 5 N = 24 to O(1) candidates).
3. **Rank 5 beyond N = 24**: try N ∈ {5, 7, 11, 15, 20} — most of these
   should have no multi-component MTC (consistent with empirical rank 5
   classification).

**Medium-term:**

4. **Rank 6 + n_θ = 3**: find the smallest real MTC with n_θ = 3 in
   NsdGOL6.g (no such example exists at rank 5). Current naive O(3)
   sweep handles this in ~20 s/prime at p = 73.
5. **Phase 4 (BNRW admissibility)**: implement Cauchy theorem
   (Galois norm of D²) and Frobenius-Schur indicators for candidate
   filtering beyond Verlinde integrality.
6. **Multiple degenerate eigenspaces**: extend `find_mtcs_at_prime` to
   handle ≥ 2 degenerate V_θ simultaneously (cartesian product of
   O(n_θ) sweeps).

**Long-term:**

7. **Algebraic solver for n_θ ≥ 4**: Cayley parametrisation + Gröbner
   basis over F_p for the Verlinde polynomial system. Target: Drinfeld
   centers D(Z_n), which have n_θ = n at rank n².
8. **Galois orbit organisation** of the final MTC catalog.
9. **Haagerup Z(H_3)** at rank 12, N = 39 — the original motivation for
   multi-component MTC search.

## Status summary (as of v0.3)

| Component | Status | Notes |
|-----|-----|-----|
| Phase 0 (catalog) | ✅ | N = 24, 103 irreps ≤ rank 5 |
| Phase 1 (strata) | ✅ | 25k strata for rank 5 N = 24 |
| Phase 2 (block-U) | ✅ | General O(n), n ≤ 3 naive, SU(2)_4 validated |
| Phase 3 (CRT) | ✅ | Galois-aware, 7/7 primes for SU(2)_4 |
| Phase 1.5 (T-filter) | 🔜 | Optimisation; not blocking |
| Phase 4 (admissibility) | 🔜 | BNRW Cauchy + FS indicators |
| n_θ ≥ 4 solver | 🔜 | Algebraic (Gröbner over F_p) required |

