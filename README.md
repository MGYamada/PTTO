# ACMG — Arithmetic Condensed Matter Geometry

Modular Tensor Category (MTC) classification by **conductor**, as a single
automated pipeline from an integer `N` to a list of fully verified MTCs
including their pentagon / hexagon `(F, R)` symbols.

```julia
using ACMG

# Fibonacci fusion ring
Nijk = zeros(Int, 2, 2, 2)
Nijk[1,1,1] = Nijk[1,2,2] = Nijk[2,1,2] = Nijk[2,2,1] = Nijk[2,2,2] = 1

result = compute_FR_from_ST(Nijk)
# result.F :: Vector{ComplexF64}  pentagon F-symbols
# result.R :: Vector{ComplexF64}  hexagon R-symbols
# result.report :: VerifyReport   pentagon/hexagon residuals < 1e-14
```

At the top level, the recommended entry point is now
`classify_mtcs_at_conductor(N; ...)`. You can call it with only `N`,
and the pipeline auto-selects `N_effective`, `primes`, and default
search parameters:

```julia
mtcs = classify_mtcs_at_conductor(24; skip_FR = true)
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
with pentagon/hexagon residuals.

## `conductor_mode`

- `conductor_mode = :full_mtc` is the only supported mode.
- Internal search uses
  `N_effective = lcm(N, cyclotomic_requirement(scale_d))`, where
  `cyclotomic_requirement(2|3)=24`, `cyclotomic_requirement(5)=5`,
  else `1` (so in many cases `N_effective = N`).
- `conductor_mode = :T_only` was removed.

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

Compatibility policy: `[compat]` uses Julia caret semantics (e.g. `1.3`, `0.5`) to keep narrow ranges with explicit major/minor upper bounds. The package version is defined in `Project.toml`.

## Worked example: Fibonacci

The Fibonacci category has rank 2, fusion ring `τ ⊗ τ = 1 ⊕ τ`, and
twists `θ_1 = 1`, `θ_τ = exp(4πi/5)`.

> **Conductor note (for this section):**
> Here `N = 5` means the conductor used when deriving the example twists.
> `compute_FR_from_ST` itself does not run a conductor search; it solves
> `(F, R)` from fusion data `Nijk`.

For Fibonacci, the pentagon system has 5 F-variables and 12 equations;
the hexagon system has 2 variables once F is fixed.

`compute_FR_from_ST` takes only a fusion tensor and performs:

1. Pentagon HC on `Nijk` — returns 4 F-solutions, split into 2
   gauge classes.
2. For each F: hexagon HC — each pentagon solution yields 2 R-solutions.
3. Returns numerically valid pentagon/hexagon `(F, R)` candidates.
   Modular-data branch selection against `T` is done later in the
   pipeline roundtrip step.

```julia
result = compute_FR_from_ST(Nijk; verbose = true)

# Result
#   pentagon HC: 4 solutions
#   n_matches: number of numerically valid pentagon/hexagon pairs
#   chosen pair: F[2] × R[2]
#   VerifyReport(rank=2,
#                pentagon_max = 3.96e-16 over 12 eqs,
#                hexagon_max  = 1.37e-15 over 29 eqs)
```

The returned `(F, R)` pair can be handed to `verify_mtc(F, R, Nijk)`
for residual verification, or fed into downstream algebraic recognition.

For diagnostics, `n_matches` reports how many pentagon/hexagon-valid
candidates were found for the given fusion ring:

```julia
compute_FR_from_ST(Nijk).n_matches    # == 0
```

## Difference between `N` in `compute_FR_from_ST` and in `classify_mtcs_at_conductor`

- On the `compute_FR_from_ST(Nijk)` side:
  - You solve pentagon/hexagon `(F, R)` on a fixed fusion ring.
  - Twist/conductor interpretation enters later when checking a chosen branch against modular data.
- On the `classify_mtcs_at_conductor(N; ...)` side:
  - You are searching candidates that satisfy consistency conditions for the **full MTC (including the S-field)** starting from `N`.
  - Under `conductor_mode = :full_mtc`, the internal search uses
    `N_effective = lcm(N, cyclotomic_requirement(scale_d))`, so the search conductor may differ
    from the conductor seen from `T` alone.
  - In practice for Fibonacci, the pipeline-side base conductor should be taken as `N = 20`
    (while the input twist conductor is still `N_T = 5`).

If you want explicit control, a typical call that reaches Fibonacci from
`N = 20` is:

```julia
# Fully automatic from base conductor N
mtcs = classify_mtcs_at_conductor(20; scale_d = 5)

# Or explicit control (optional):
mtcs_explicit = classify_mtcs_at_conductor(
    20;
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

### Pentagon / Hexagon / modular-data consistency

Phase 4 takes `(S, T)` in ℂ and a fusion tensor `Nijk`, and returns
`(F, R)` symbols satisfying:

- **Pentagon**: `TensorCategories.pentagon_equations(Nijk)` →
  homotopy continuation with a random linear slice to break gauge
  symmetry → damped Newton refinement.
- **Hexagon**: `hexagon_equations(Nijk, one_vec, F)` (ACMG's own
  implementation with F baked in and the `R·S = I` constraint) →
  HC again.
- **Ribbon residual diagnostic**: the multiplicity-free relation
  `(R^{ij}_k)² = θ_i · θ_j / θ_k`
  is evaluated for all hexagon solutions, and the minimum-residual pair
  is selected.
- **Modular-data roundtrip check**: from `(F, R, N)` we evaluate
  consistency of reproduced modular data `(S, T)` against the Phase 3
  lifted target `(S, T)` up to cyclotomic Galois action (diagnostic log).

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

`BlockU.jl`. `find_mtcs_at_prime(catalog, stratum, p; ...)` returns
`Vector{MTCCandidate}`. Steps:

1. Assemble block-diagonal atomic `(S, T)` in `ℤ[ζ_N]` from
   `stratum × catalog`.
2. Reduce to F_p at a chosen primitive N-th root of unity
   `find_zeta_in_Fp(N, p)`.
3. Decompose T into eigenspaces; compute `parameter_dim`.
4. For each degenerate eigenspace `V_θ` with `n_θ ≥ 2`, run Phase-2
   block search backend:
   - `search_mode = :groebner` (default): best-effort Gröbner system
     build + point extraction (`O(n)` orthogonality and fixed-unit
     Verlinde/Cayley-linked systems), then Cayley unit-axiom filtered
     search, with optional final fallback to enumeration.
   - `search_mode = :exhaustive`: direct Cayley/reflection sweep over
     `O(n_θ)(F_p)`.
5. For each block-U, apply to S, check Verlinde integrality with
   `verlinde_find_unit`, build `MTCCandidate`.

Useful knobs:
- `max_units_for_groebner`: cap number of fixed-unit systems attempted.
- `groebner_allow_fallback=true`: if set `false`, solver mode does not
  fall back to exhaustive enumeration when extraction is empty.
- `precheck_unit_axiom=true`: run a fast unit-axiom prefilter before
  full Verlinde tensor construction.

Practical presets:
- **Strict solver-only trial**:
  `search_mode=:groebner, groebner_allow_fallback=false, precheck_unit_axiom=true`
- **Hybrid (recommended default)**:
  `search_mode=:groebner, groebner_allow_fallback=true, precheck_unit_axiom=true`
- **Legacy exhaustive baseline**:
  `search_mode=:exhaustive`

Feasibility: `|O(n)(F_p)| ~ p^{n(n-1)/2}`. At p ≈ 100, exhaustive
`n = 3` is still expensive; `n ≥ 4` generally needs solver-assisted
search to be practical.

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

### Phase 4 — Pentagon + Hexagon + modular-data verify

Six files, all at `src/` top level:

- `PentagonEquations.jl`: thin wrapper over
  `TensorCategories.pentagon_equations`.
- `PentagonSolver.jl`: HC (`solve_pentagon_homotopy`) and damped Newton
  (`solve_pentagon_newton` / `refine_solution_newton`).
- `HexagonEquations.jl`: custom generator that bakes numerical F-values
  into the coefficient ring (`AcbField`), adds the `R · S = I` constraint.
  Public helper functions are `coerce_complex` and
  `number_of_variables_in_hexagon_equations` (internal `_`-prefixed
  variants are no longer exported).
- `HexagonSolver.jl`: HC for the R-system.
- `ModularDataLift.jl`: discrete-log table `DiscreteLogTable(N, p, ζ_Fp)`
  to lift T from F_p back to ℂ; `ℤ[√d] → ℂ` for S.
- `Verify.jl`:
    - `pentagon_residuals(F, Nijk)`
    - `hexagon_residuals(F, R, Nijk)`
    - `verify_mtc(F, R, Nijk) -> VerifyReport`

All residuals are ∞-norm at `Float64` precision; a clean solution
scores < 1e-10 on both pentagon and hexagon.

### Phase 5 — End-to-end pipeline

`Pipeline.jl`.

- `compute_FR_from_ST(Nijk; ...)`: solve pentagon via HC, solve
  hexagon for each F, and return numerically valid `(F, R)` candidates.
  Returns a NamedTuple with the selected pair and a `VerifyReport`.
- `classify_from_group(group, N, stratum, primes; ...)`: CRT + ℂ-lift
  + Phase 4, producing one `ClassifiedMTC`.
- `classify_mtcs_at_conductor(N; max_rank = 5, primes = nothing, ...)`: full driver,
  iterating over strata and Galois sectors.
  - Note: using an incorrect `scale_d` can produce **zero candidates**.
    For Fibonacci conductor search (`N = 20`), we recommend `scale_d = 5`.

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
- Pentagon HC + hexagon HC on Fibonacci.
- `DiscreteLogTable`, T and S complex lifts.
- **Phase 5**: `compute_FR_from_ST` on Fibonacci plus a multi-prime
  integration smoke test of the Phase 0 → 3 driver at larger conductor
  with `skip_FR = true`.

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
   continuous block-U moduli to a 0-dimensional locus before Phase 4.

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
  test_verify_residuals.jl                — residuals and VerifyReport
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
