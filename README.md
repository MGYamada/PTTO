# ACMG — Arithmetic Condensed Matter Geometry

Modular Tensor Category (MTC) classification by **conductor**, as a single
automated pipeline from an integer `N` to a list of fully verified MTCs
including their pentagon / hexagon `(F, R)` symbols.

```julia
using ACMG

# Fibonacci fusion ring
Nijk = zeros(Int, 2, 2, 2)
Nijk[1,1,1] = Nijk[1,2,2] = Nijk[2,1,2] = Nijk[2,2,1] = Nijk[2,2,2] = 1

# Fibonacci T
T = ComplexF64[1.0, exp(4π * im / 5)]

result = compute_FR_from_ST(Nijk, T)
# result.F :: Vector{ComplexF64}  pentagon F-symbols
# result.R :: Vector{ComplexF64}  hexagon R-symbols
# result.report :: VerifyReport   pentagon/hexagon/ribbon residuals < 1e-14
```

At the top level, the recommended entry point is
`classify_mtcs_auto(N; ...)`, which runs all five phases and
auto-selects `conductor_mode`, `scale_d`, `primes`, and `max_rank`
from built-in candidate lists:

```julia
auto = classify_mtcs_auto(24; skip_FR = true)
classified = auto.classified

# reproducibility metadata chosen by the auto loop
auto.N_effective
auto.scale_d
auto.primes
```

Under the hood this calls `classify_mtcs_at_conductor(N; ...)`:

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
  Vector{ClassifiedMTC}
```

Each `ClassifiedMTC` carries the stratum decomposition, the exact
arithmetic S-matrix in ℤ[√d], its complex-lifted `(S, T)`, the
`(F, R)` symbols solving pentagon/hexagon over ℂ, and a `VerifyReport`
with pentagon, hexagon and ribbon residuals.

## `conductor_mode` (v0.5.0)

- `conductor_mode = :full_mtc` is the only supported mode.
- Internal search uses `N_effective = lcm(N, 4 * scale_d)`.
- `conductor_mode = :T_only` was removed in v0.5.0.

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

Compatibility policy: `[compat]` uses Julia caret semantics (e.g. `1.3`, `0.8`) to keep narrow ranges with explicit major/minor upper bounds for v0.5.0.

## Worked example: Fibonacci

The Fibonacci category has rank 2, fusion ring `τ ⊗ τ = 1 ⊕ τ`, and
twists `θ_1 = 1`, `θ_τ = exp(4πi/5)`.

> **Conductor note (for this section):**
> Here `N = 5` means the **conductor of the input T-spectrum** used by
> `compute_FR_from_ST`. This function does not run a conductor search;
> it only solves `(F, R)` for the given `(Nijk, T)`.

For Fibonacci, the pentagon system has 5 F-variables and 12 equations;
the hexagon system has 2 variables once F is fixed.

`compute_FR_from_ST` takes a fusion tensor and complex T-eigenvalues
and performs:

1. Pentagon HC on `Nijk` — returns 4 F-solutions, split into 2
   gauge classes.
2. For each F: hexagon HC — each pentagon solution yields 2 R-solutions.
3. Ribbon match `(R^{ij}_k)² = θ_i θ_j / θ_k` against the input `T`
   picks the `(F, R)` pair(s) realising that specific modular datum.

```julia
result = compute_FR_from_ST(Nijk, T_fib; ribbon_atol = 1e-8, verbose = true)

# Result
#   pentagon HC: 4 solutions
#   ribbon matches: 2 of 8 (F,R) pairs
#   chosen pair: F[2] × R[2]
#   VerifyReport(rank=2,
#                pentagon_max = 3.96e-16 over 12 eqs,
#                hexagon_max  = 1.37e-15 over 29 eqs,
#                ribbon_max   = 2.25e-15)
```

The returned `(F, R)` pair can be handed to `verify_mtc(F, R, Nijk; T)`
independently for re-verification, or fed into downstream algebraic
recognition.

As a sanity check, running the same fusion ring against a bogus
`T = (1, i)` — a semion-like twist not consistent with the Fibonacci
pentagon class — returns zero ribbon matches:

```julia
compute_FR_from_ST(Nijk, ComplexF64[1.0, im]).n_matches    # == 0
```

## Difference between `N` in `compute_FR_from_ST` and in `classify_mtcs_at_conductor`

- On the `compute_FR_from_ST(Nijk, T)` side:
  - You are matching `(F, R)` on a fixed fusion ring using the **conductor of the input T-spectrum** (root-of-unity phases in `T`).
  - In the Fibonacci example, `T = (1, exp(4πi/5))`, so `N_T = 5`.
- On the `classify_mtcs_at_conductor(N; ...)` side:
  - You are searching candidates that satisfy consistency conditions for the **full MTC (including the S-field)** starting from `N`.
  - Under `conductor_mode = :full_mtc`, the internal search uses
    `N_effective = lcm(N, 4 * scale_d)`, so the search conductor may differ
    from the conductor seen from `T` alone.

If you want explicit control, a typical call that reaches Fibonacci from
`N = 5` is:

```julia
mtcs = classify_mtcs_at_conductor(
    5;
    max_rank = 2,
    primes = [41, 61],
    conductor_mode = :full_mtc,
    scale_d = 5,
)
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
  give different MTCs; the ribbon check selects the one matching the
  input T.

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

`Pipeline.jl`.

- `compute_FR_from_ST(Nijk, T; ...)`: solve pentagon via HC, solve
  hexagon for each F, select `(F, R)` whose ribbon residual vs. `T`
  is below `ribbon_atol`. Returns a NamedTuple with the matched pair
  and a `VerifyReport`.
- `classify_from_group(group, N, stratum, primes; ...)`: CRT + ℂ-lift
  + Phase 4, producing one `ClassifiedMTC`.
- `classify_mtcs_at_conductor(N; max_rank, primes, ...)`: full driver,
  iterating over strata and Galois sectors.
  - Note: using an incorrect `scale_d` can produce **zero candidates**.
    For Fibonacci (`N = 5`), we recommend `scale_d = 5`.

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

Pentagon HC scales with the number of F-variables via the mixed volume
of the polynomial system; this is the current bottleneck of the
pipeline.

At rank 2 (5 F-vars) pentagon HC finishes in tens of seconds and has
been validated end-to-end on Fibonacci and its sign-conjugate
Yang-Lee. Rank-3 systems (~14 F-vars) have mixed volume above `10⁵`
and push HC to the practical edge. At rank 5 (~238 F-vars) the
mixed volume blows up further; running Phase 4 at this scale requires
a gauge-fixed pentagon formulation that is on the roadmap but not yet
implemented. When this is the case, `classify_mtcs_at_conductor` can
be called with `skip_FR = true` to produce modular data (Phases 0–3
plus the ℂ-lift) without attempting the pentagon/hexagon stage.

The F_p block-U sweep has its own frontier: `|O(n_θ)(F_p)| ~ p^{n_θ(n_θ-1)/2}`,
so `n_θ = 3` is the practical limit of the current naive enumeration.
`n_θ = 4` requires a polynomial solver over F_p (Cayley parametrisation
+ Gröbner).

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
- **Phase 5**: `compute_FR_from_ST` on Fibonacci — including a negative
  test with wrong T — plus a multi-prime integration smoke test of the
  Phase 0 → 3 driver at larger conductor with `skip_FR = true`.

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
  Types.jl                — core layer: ModularDatumFp/FusionRule + F_p primitives +
                            (S, T) validation + Verlinde extraction/lift
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
  Pipeline.jl             — end-to-end pipeline: classify_mtcs_at_conductor and friends

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
  test_fr_pentagon_hexagon_fibonacci.jl  — pentagon + hexagon HC
  test_modular_data_lift.jl               — DiscreteLogTable, T and S lifts
  test_verify_residuals.jl                — residuals, ribbon match, VerifyReport
  test_pipeline_end_to_end.jl             — pipeline end-to-end and FR integration
  test_pipeline_galois_grouping.jl        — CRT branch consistency + sector grouping
  test_pipeline_primes.jl                 — prime selection and conductor mode checks

scripts/
  phase5_demo.jl          — Phase 5 end-to-end demo (Fibonacci)
  su24_integration.jl     — Phase 2 multi-prime driver (legacy)
  su24_crt.jl             — Phase 3 end-to-end CRT (legacy)
  diagnose_crt.jl         — per-prime candidate inspection
```

## Roadmap

**Near-term**

- T-spectrum pre-filter `enumerate_strata_by_T(catalog, target_spins)`
  to collapse large strata enumerations to O(1) candidates.
- BNRW admissibility layer (Cauchy / Galois norm of D² and
  Frobenius-Schur indicators) as a pre-Verlinde filter.

**Medium-term**

- Gauge-fixed pentagon formulation to bring rank-3 and rank-5 systems
  into HC range.
- Multiple degenerate eigenspaces handled simultaneously (Cartesian
  product of O(n_θ) sweeps).

**Long-term**

- Algebraic solver for `O(n)` sweeps at `n ≥ 4` (Cayley + Gröbner over F_p).
- Algebraic recognition of `ComplexF64` F-symbols as elements of ℚ(ζ_N)
  (PSLQ / LLL).

## References

- Bruillard–Ng–Rowell–Wang, *On classification of modular categories by rank*,
  arXiv:1507.05139 (admissibility).
- Ng–Rowell–Wang–Wen, *Reconstruction of modular data from SL_2(ℤ)
  representations*, arXiv:2203.14829 (the irrep-sum + block-U
  framework).
- Ng–Rowell–Wen, *Classification of modular data up to rank 12*,
  arXiv:2308.09670 (ground-truth catalog, `NsdGOL*.g` ancillary files).

## Status

| Phase | What it does                           | Status |
|-------|-----------------------------------------|--------|
| 0     | SL(2, ℤ/N) atomic catalog               | ✅     |
| 1     | Stratum enumeration                     | ✅     |
| 2     | Block-U sweep + Verlinde (F_p)          | ✅     |
| 3     | CRT + Galois-aware reconstruction       | ✅     |
| 4     | Pentagon + Hexagon + Ribbon (ℂ)         | ✅     |
| 5     | End-to-end `classify_mtcs_at_conductor` | ✅     |

Version: **v0.4-prototype** (April 2026).
