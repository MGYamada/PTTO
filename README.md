# ACMG — Arithmetic Condensed Matter Geometry

Modular Tensor Category (MTC) classification by **conductor**, as a single
automated pipeline from an integer `N` to a list of fully verified MTCs
including their pentagon / hexagon `(F, R)` symbols.

```julia
using ACMG

classified = classify_mtcs_at_conductor(24;
                                         max_rank = 5,
                                         primes = [73, 97, 193, 241,
                                                   313, 337, 409],
                                         skip_FR = true)   # rank-5 HC infeasible
# → 2 ClassifiedMTCs, one per Galois sector of SU(2)_4.
```

The single call above runs all five phases of the pipeline:

```
 N ─┬─ Phase 0:  Atomic SL(2, ℤ/N)-irreps  (Oscar + GAP/SL2Reps)
    │
    ├─ Phase 1:  Stratum enumeration       (partitions m_λ of rank r)
    │
    ├─ Phase 2:  Block-U sweep over F_p    (O(n) Cayley + Verlinde check)
    │
    ├─ Phase 3:  CRT → ℤ[√d]               (Galois-aware, cross-validated)
    │
    └─ Phase 4:  Pentagon + Hexagon + Ribbon over ℂ
                                           (TensorCategories + HC)
  ↓
  List[ClassifiedMTC]
```

Each `ClassifiedMTC` carries the stratum decomposition, the exact
arithmetic S-matrix in ℤ[√d], its complex-lifted `(S, T)`, the
`(F, R)` symbols solving pentagon/hexagon over ℂ, and a `VerifyReport`
with pentagon, hexagon and ribbon residuals.

## Installation

```julia
julia> ]
pkg> activate .
pkg> instantiate
```

ACMG depends on:

- **Oscar.jl** — cyclotomic fields and GAP/SL2Reps access (Phase 0)
- **TensorCategories.jl** — pentagon equation generator (Phase 4)
- **HomotopyContinuation.jl** — pentagon/hexagon solver (Phase 4)
- **KrylovKit.jl** — damped Newton refinement (Phase 4)
- **Primes.jl**, **LinearAlgebra**, **SparseArrays**

Julia ≥ 1.9 required.

## Quick start

### 1. Fibonacci: `(S, T) → (F, R)` in 30 seconds

```julia
using ACMG

# Fibonacci fusion ring
Nijk = zeros(Int, 2, 2, 2)
Nijk[1,1,1] = Nijk[1,2,2] = Nijk[2,1,2] = Nijk[2,2,1] = Nijk[2,2,2] = 1

# Fibonacci T
T = ComplexF64[1.0, exp(4π * im / 5)]

result = compute_FR_from_ST(Nijk, T)
# result.F :: Vector{ComplexF64}  pentagon F-symbol values
# result.R :: Vector{ComplexF64}  hexagon R-symbol values
# result.report :: VerifyReport   pentagon/hexagon/ribbon residuals < 1e-14
```

### 2. SU(2)_4: full `N → MTCs` pipeline

```julia
classified = classify_mtcs_at_conductor(24;
                                         max_rank = 5,
                                         primes = [73, 97, 193, 241,
                                                   313, 337, 409],
                                         skip_FR = true)

for c in classified
    println(c)
    # Exact arithmetic S-matrix:
    println(describe_matrix(c.S_Zsqrtd, c.scale_d))
end
```

`skip_FR = true` is required here because rank-5 pentagon systems have
~238 F-variables and HC is not currently feasible at that scale
(see "Complexity frontier" below).

### 3. Verify a candidate `(F, R)` pair

```julia
report = verify_mtc(F_values, R_values, Nijk; T = T_complex)
# report.pentagon_max, report.hexagon_max, report.ribbon_max
```

## Conceptual framework

### Why conductor-first

The moduli space of rank-`r` modular data of conductor dividing `N`
sits inside the variety of SL(2, ℤ/N)-representations of rank `r`:

    M_{r,N}^{MTC}  ⊂  M_{r,N}^{SL_2(ℤ/N)-rep}

Fixing `N` restricts modular-data entries to `ℤ[ζ_N] ∩ {|z| ≤ D}`, a
finite algebraic set, and aligns the search with Ng–Rowell–Wang–Wen's
(NRWW, arXiv:2203.14829) stratification by SL(2, ℤ/N)-irrep content.

Fixing rank first — as in the classical NRW rank-by-rank tables —
forces unnecessary field extensions when the MTC actually lives in a
smaller cyclotomic field. Conductor-first avoids this.

### Stratification by SL(2, ℤ/N)-irrep decomposition

Each isomorphism class of representation

    ρ  =  ⊕_λ  ρ_λ^{m_λ}

(a partition of `r` with multiplicities `m_λ`) defines a locally closed
stratum of modular data. Phase 1 enumerates these strata; for small `N`
and small rank this is a small combinatorial set.

### Continuous moduli inside a stratum

Even when all `m_λ = 1` — so the stratum's automorphism group is just
a torus acting trivially on S — a continuous moduli appears whenever
different irreps share a T-eigenvalue `θ`. The T-eigenspace

    V_θ  =  ⊕_λ V_θ^{(ρ_λ)}

may have `n_θ ≥ 2`, and rotating its basis by an element of O(n_θ) gives
a new S-matrix with the same T and the same representation-theoretic
structure. The continuous moduli has dimension

    parameter_dim  =  Σ_θ  (n_θ choose 2).

### Verlinde integrality: from continuous to discrete

On the continuous moduli, the Verlinde coefficients
`N_{ij}^k(φ) = Σ_m S_{im}S_{jm}S̄_{km}/S_{1m}` are algebraic functions
of the block-U parameter. The MTC locus is the subset where every
`N_{ij}^k` is a non-negative integer — generically a 0-dimensional
algebraic set inside the continuous family.

For SU(2)_4 at rank 5, `N = 24`, one T-eigenspace has `n_θ = 2`, the
continuous family is `O(2) = P^1`, and integrality picks out two
points corresponding to a rotation-reflection pair
(φ = π/4; derivable analytically from Kac–Peterson).

### Galois-aware F_p validation

Computing over `ℤ[ζ_N]` directly is expensive. Instead ACMG reduces to
F_p at primes `p` with `N | p−1`, runs the block-U sweep and
Verlinde integrality check in F_p exactly, and reconstructs algebraic
entries via multi-prime CRT.

Two primes picking up the same abstract MTC may return different Galois
conjugates (the F_p enumeration order and the choice of primitive root
varies across primes). Fusion tensors are Galois-invariant, but
S-matrices are not. Phase 3's
`group_mtcs_galois_aware(results, anchor_prime)` aligns candidates
across primes into a single Galois sector using 2-prime trial
reconstruction as the matching criterion — a construction not present
in NRWW but necessary once you compute in F_p rather than a global
number field.

### Pentagon / Hexagon / Ribbon

Phase 4 takes `(S, T)` in ℂ and a fusion tensor `Nijk`, and returns
`(F, R)` symbols satisfying:

- **Pentagon**: `TensorCategories.pentagon_equations(Nijk)` →
  homotopy continuation with a random linear slice to break gauge
  symmetry → damped Newton refinement.
- **Hexagon**: `hexagon_equations(Nijk, one_vec, F)` (ACMG's own
  implementation with F baked in and the `R·S = I` constraint) →
  HC again.
- **Ribbon match**: the multiplicity-free relation
  `(R^{ij}_k)² = θ_i · θ_j / θ_k`
  is a necessary condition pinning down which `(F, R)` realises the
  target T. Different pentagon F-classes on the same fusion ring
  give different MTCs (e.g. Fibonacci vs Yang-Lee); the ribbon check
  selects the one matching the input T.

## Pipeline anatomy

### Phase 0 — Atomic SL(2, ℤ/N) irrep catalog

`SL2Reps.jl`. Builds `Vector{AtomicIrrep}` via Oscar's GAP interface to
the SL2Reps package. Each irrep records its dimension, level, parity,
T-eigenvalues (`T_powers` ∈ `Vector{Int}` — exponents of `ζ_N`), and
S-matrix in `ℤ[ζ_N]`. Restricted by `max_rank` for efficiency.

### Phase 1 — Stratum enumeration

`StratumEnum.jl`. `enumerate_strata(catalog, r)` returns all
`Stratum` structs `(m_λ)` with `Σ m_λ · dim(ρ_λ) = r`.

Default is `require_unit_summand = false` — in non-pointed MTCs the
unit object sits *inside* a larger irrep wherever T-eigenvalue 1 appears,
not as a separate 1d_1 summand.

### Phase 2 — Block-U sweep at a single prime

`BlockU.jl`. `find_mtcs_at_prime(catalog, stratum, p)` returns
`Vector{MTCCandidate}`. Steps:

1. Assemble block-diagonal atomic `(S, T)` in `ℤ[ζ_N]` from
   `stratum × catalog`.
2. Reduce to F_p at a chosen primitive N-th root of unity
   `find_zeta_in_Fp(N, p)`.
3. Decompose T into eigenspaces; compute `parameter_dim`.
4. For each degenerate eigenspace `V_θ` with `n_θ ≥ 2`, sweep over
   `O(n_θ)(F_p)` via the Cayley parametrisation.
5. For each block-U, apply to S, check Verlinde integrality with
   `verlinde_find_unit`, build `MTCCandidate`.

Feasibility: `|O(n)(F_p)| ~ p^{n(n-1)/2}`. At p ≈ 100, `n = 3` finishes
in ~20 s; `n ≥ 4` requires an algebraic solver and is not yet
implemented.

### Phase 3 — CRT reconstruction and Galois alignment

`CRT.jl`.

- `group_mtcs_galois_aware(results_by_prime, anchor)`: aligns
  Galois-conjugate candidates across primes into one group per
  sector using 2-prime trial reconstruction.
- `reconstruct_S_matrix(group; scale_d, sqrtd_fn, bound)`: recovers
  `scale · √d · S` entry-wise in `ℤ[√d]` via multi-prime CRT and
  rational reconstruction, returning `Matrix{Tuple{Int, Int}}` with
  entries `(a, b) = a + b·√d`.
- `compute_sqrt3_cyclotomic_mod_p(p)` /
  `compute_sqrt2_cyclotomic_mod_p(p)`: use `ζ_N^k + ζ_N^{-k}` for a
  Galois-consistent branch of `√d` across primes.
- `verify_reconstruction(recon, candidate, d)`: validates the lift
  against a fresh prime.

### Phase 4 — Pentagon + Hexagon + Ribbon verify

Six files, all at `src/` top level:

- `PentagonEquations.jl`: thin wrapper over
  `TensorCategories.pentagon_equations`.
- `PentagonSolver.jl`: HC (`solve_pentagon_homotopy`) and damped Newton
  (`solve_pentagon_newton` / `refine_solution_newton`).
- `HexagonEquations.jl`: custom generator that bakes numerical F-values
  into the coefficient ring (`AcbField`), adds the `R · S = I` constraint.
- `HexagonSolver.jl`: HC for the R-system.
- `ModularDataLift.jl`: discrete-log table `DiscreteLogTable(N, p, ζ_Fp)`
  to lift T from F_p back to ℂ; `ℤ[√d] → ℂ` for S.
- `Verify.jl`:
    - `pentagon_residuals(F, Nijk)`
    - `hexagon_residuals(F, R, Nijk)`
    - `ribbon_residuals(R, T, Nijk)` checking `(R^{ij}_k)² = θ_i θ_j / θ_k`
    - `verify_mtc(F, R, Nijk; T) -> VerifyReport`

All residuals are ∞-norm at `Float64` precision; a clean solution
scores < 1e-10 across all three.

### Phase 5 — End-to-end pipeline

`Phase5.jl`.

- `compute_FR_from_ST(Nijk, T; ...)`: solve pentagon via HC, solve
  hexagon for each F, select `(F, R)` whose ribbon residual vs. `T`
  is below `ribbon_atol`. Returns a NamedTuple with the matched pair
  and a `VerifyReport`.
- `classify_from_group(group, N, stratum, primes; ...)`: CRT + ℂ-lift
  + Phase 4, producing one `ClassifiedMTC`.
- `classify_mtcs_at_conductor(N; max_rank, primes, ...)`: full driver,
  iterating over strata and Galois sectors.

## Output type

```julia
struct ClassifiedMTC
    N::Int                              # conductor
    rank::Int
    stratum::Stratum                    # (m_λ) decomposition
    Nijk::Array{Int, 3}                 # fusion tensor
    S_Zsqrtd::Matrix{Tuple{Int, Int}}   # (a, b) = a + b·√d
    scale_d::Int                        # d for ℤ[√d]
    scale_factor::Int                   # factor so S_ℂ = S_Zsqrtd/(scale·√d)
    used_primes::Vector{Int}
    fresh_primes::Vector{Int}
    verify_fresh::Bool                  # ✓ if all fresh primes cross-check
    S_complex::Matrix{ComplexF64}
    T_complex::Vector{ComplexF64}       # θ_i, N-th roots of unity
    F_values::Union{Vector{ComplexF64}, Nothing}
    R_values::Union{Vector{ComplexF64}, Nothing}
    verify_report::Union{VerifyReport, Nothing}
    galois_sector::Int
end
```

`F_values`, `R_values`, `verify_report` are `nothing` iff `skip_FR=true`
was passed (or Phase 4 couldn't find a match).

## Complexity frontier

| rank | conductor | pentagon F-vars | pentagon HC | notes                      |
|------|-----------|-----------------|-------------|----------------------------|
| 2    | 5         | 5               | ~30 s       | Fibonacci ✓                |
| 3    | 16        | 14              | ~hours      | Ising, mixed volume ~5·10⁵ |
| 5    | 24        | 238             | infeasible  | SU(2)_4 — use `skip_FR`    |

| n_θ | O(n_θ)(F_p) at p ≈ 100 | wall clock   | status      |
|-----|-------------------------|--------------|-------------|
| 2   | ~150                   | instant      | ✓           |
| 3   | ~7·10⁵                 | ~20 s        | ✓           |
| 4   | ~10¹¹                  | infeasible   | needs Gröbner over F_p |

Both frontiers motivate future work on a gauge-fixed pentagon
formulation and an algebraic `O(n)` solver.

## Validation

Test suite (`julia --project=. -e 'using Pkg; Pkg.test()'`) covers:

- F_p arithmetic primitives, Fibonacci (N=5, p=41) and Ising (N=16, p=17)
  as hand-verified cases.
- Atomic catalog construction at N ∈ {5, 8, 16}.
- Stratum enumeration combinatorics.
- O(n) Cayley parametrisation and single-prime block-U sweep.
- CRT primitives, Galois-aware grouping, ℤ[√d] reconstruction with
  fresh-prime cross-validation.
- Pentagon HC + hexagon HC on Fibonacci (end-to-end ribbon match).
- `DiscreteLogTable`, T and S complex lifts.
- **Phase 5**: `compute_FR_from_ST` on Fibonacci (including a negative
  test with wrong T), and `classify_mtcs_at_conductor` on SU(2)_4 at
  `N=24` in `skip_FR` mode.

## Prime selection

For conductor `N`, primes `p` with `N | p-1`:

| N  | Example primes                        |
|----|---------------------------------------|
| 5  | 11, 31, 41, 61, 71                    |
| 8  | 17, 41, 73, 89, 97                    |
| 12 | 13, 37, 61, 73, 97                    |
| 24 | 73, 97, 193, 241, 313, 337, 409       |

A prime is admissible for a particular candidate if additionally
`√d ∈ F_p` for every quadratic irrationality appearing in the MTC
(Chebotarev density guarantees positive-density sets). In practice 4–7
good primes are enough for CRT + cross-validation.

## Design principles

1. **Conductor-first.** The outer parameter is `N`, and all arithmetic
   stays inside `ℚ(ζ_N)`. Rank emerges from stratum enumeration.

2. **F_p as a filter, `ℤ[ζ_N]` as ground truth.** F_p enumeration is
   fast, exact, and parallel; CRT reconstructs the algebraic answer.

3. **ℂ as the workspace for `(F, R)`.** Over ℂ there is topology (enabling
   damping and line search) and TensorCategories' polynomial rings are
   available. Algebraic lift to `ℚ(ζ_N)` is deferred.

4. **Galois awareness.** F_p computation mixes Galois conjugates;
   `group_mtcs_galois_aware` restores sector coherence before CRT.

5. **Necessary conditions first.** Verlinde integrality cuts the
   continuous block-U moduli to a 0-dimensional locus; ribbon
   selects a specific `(F, R)` class.

## Repository layout

```
src/
  ACMG.jl                 — module root, exports
  FpArith.jl              — F_p primitives
  Types.jl                — ModularDatumFp, FusionRule
  ModularData.jl          — (S, T) axiom validation
  FusionExtract.jl        — F_p Verlinde → ℤ lift
  Dimensions.jl           — (stub) quantum dimension enumeration
  Enumerator.jl           — (stub) legacy top-level driver
  SL2Reps.jl              — Phase 0: Oscar + GAP/SL2Reps catalog
  StratumEnum.jl          — Phase 1: (m_λ) enumeration
  BlockU.jl               — Phase 2: O(n) Cayley + single-prime driver
  CRT.jl                  — Phase 3: CRT + Galois-aware grouping
  PentagonEquations.jl    — Phase 4: TensorCategories wrapper
  PentagonSolver.jl       — Phase 4: HC + damped Newton
  HexagonEquations.jl     — Phase 4: custom generator with F fixed
  HexagonSolver.jl        — Phase 4: HC for R-system
  ModularDataLift.jl      — Phase 4: F_p / ℤ[√d] → ℂ
  Verify.jl               — Phase 4: residuals, VerifyReport
  Phase5.jl               — Phase 5: classify_mtcs_at_conductor and friends

test/
  runtests.jl
  test_fparith.jl
  test_fibonacci.jl
  test_ising.jl
  test_sl2reps.jl
  test_stratum_enum.jl
  test_blocku.jl            — pure F_p helpers
  test_block_u_general.jl   — general O(n) Cayley
  test_crt.jl               — CRT + ℤ[√d] reconstruction
  test_phase4_fibonacci.jl  — pentagon + hexagon HC
  test_phase4_lift.jl       — DiscreteLogTable, T and S lifts
  test_phase4_verify.jl     — residuals, ribbon match, VerifyReport
  test_phase5.jl            — Phase 5 pipeline: compute_FR_from_ST,
                              classify_mtcs_at_conductor (SU(2)_4, skip_FR)

scripts/
  su24_integration.jl     — Phase 2 end-to-end for SU(2)_4
  su24_crt.jl             — Phase 3 end-to-end (multi-prime CRT)
  phase5_demo.jl          — Phase 5 end-to-end: Fibonacci (F,R)
                            + SU(2)_4 full pipeline. Pass
                            --full-enumerate to sweep all strata.
  diagnose_crt.jl         — per-prime candidate inspection
```

## Roadmap

**Near-term**

- Rank 5 classification at N ∈ {5, 7, 11, 15, 20} to verify the expected
  empty set.
- T-spectrum pre-filter `enumerate_strata_by_T(catalog, target_spins)`
  to collapse the ~25k rank-5 N=24 strata to O(1) candidates.
- BNRW admissibility layer (Cauchy / Galois norm of D² and
  Frobenius-Schur indicators) as a pre-Verlinde filter.

**Medium-term**

- Gauge-fixed pentagon formulation to bring rank-3 (Ising, 14 F-vars)
  and ideally rank-5 (SU(2)_4, 238 F-vars) into reach.
- Multiple degenerate eigenspaces handled simultaneously (Cartesian
  product of O(n_θ) sweeps).
- NsdGOL6 parsing and the first `n_θ = 3` target from the NRW catalogue.
- Rank 6+ D(ℤ_n) Drinfeld centres: the `n_θ = n` at rank `n²` case.

**Long-term**

- Algebraic solver for `O(n)` sweeps at `n ≥ 4` (Cayley + Gröbner over F_p).
- Algebraic recognition of `ComplexF64` F-symbols as elements of ℚ(ζ_N)
  (PSLQ / LLL).
- Haagerup `Z(H_3)` at rank 12, N = 39.

## References

- Bruillard–Ng–Rowell–Wang, *Rank-finiteness for modular categories*,
  arXiv:1507.05139 (admissibility).
- Ng–Rowell–Wang–Wen, *Reconstruction of modular data from SL(2, ℤ)
  representations*, arXiv:2203.14829 (the irrep-sum + block-U
  framework).
- Ng–Rowell–Wen, *Classification of modular data up to rank 12*,
  arXiv:2308.09670 (ground-truth catalog, `NsdGOL*.g` ancillary files).
- Gannon–Morrison, arXiv:1606.07165 (algorithmic template).
- Evans–Gannon, arXiv:1006.1326 (Haagerup family).
- Kitaev, unpublished notes on topological phases (tangent cohomology,
  Eq. 244–252).

## Status

| Phase | What it does                           | Status | Notes                                  |
|-------|-----------------------------------------|--------|----------------------------------------|
| 0     | SL(2, ℤ/N) atomic catalog               | ✅     | N=24: 178 irreps, max_rank≤5           |
| 1     | Stratum enumeration                     | ✅     | 25k strata at rank 5 N=24              |
| 2     | Block-U sweep + Verlinde (F_p)          | ✅     | O(n) for n ≤ 3; SU(2)_4 verified       |
| 3     | CRT + Galois-aware reconstruction       | ✅     | 7/7 primes on SU(2)_4                  |
| 4     | Pentagon + Hexagon + Ribbon (ℂ)         | ✅     | Fibonacci verified; rank-5 skip        |
| 5     | End-to-end `classify_mtcs_at_conductor` | ✅     | N=24 SU(2)_4 end-to-end (skip_FR)      |

Version: **v0.4-prototype** (April 2026).
