@doc """
    CyclotomicContext

Cyclotomic ground-field context for computations over `Q(ζ_N)`.

# Fields
- `N`: Conductor used to construct the cyclotomic field.
- `field`: Oscar cyclotomic field parent.
- `zeta`: Chosen primitive `N`-th root of unity.
- `real_subfield`: Reserved slot for future real-subfield data.
- `conductor`: Stored conductor, equal to `N` for the current constructor.

# Notes
Stable public API.  In ACMG, `N` is a ground-field choice, not only a post-hoc
invariant inferred from `S` and `T`.
""" CyclotomicContext

@doc """
    modular_data(name; conductor = nothing)

Construct one of ACMG's built-in exact modular-data examples.

# Arguments
- `name`: `:semion`, `:fibonacci`, `:ising`, or `:toric_code` and accepted
  aliases.

# Keyword arguments
- `conductor`: Optional conductor override.  It must match the conductor
  required by the selected built-in example.

# Returns
A `ModularData` object over the selected `CyclotomicContext`.
""" modular_data

@doc """
    HigherCentralChargeResult

Legacy structured Gauss-sum result returned by `higher_central_charge_result`.
The public HCC moment API is `higher_central_charge(data, n)`, which uses
D² normalization and returns the exact value directly.

# Fields
- `ok`: Whether the requested normalization was available.
- `n`: Power used in the higher Gauss sum.
- `normalization`: Normalization convention.
- `value`: Normalized value when `ok` is true.
- `gauss_sum`: Raw exact Gauss sum.
- `denominator`: Exact denominator used for normalization.
- `D_squared`: Exact global-dimension square.
- `conductor`: Conductor used for the computation.
- `status`: Machine-readable status symbol.
- `message`: Human-readable explanation.
""" HigherCentralChargeResult

@doc """
    HCCGeneratingFunction

Structured HCC generating function with exact `weights` and `twists`.
Calling the object on an integer `n` returns the D²-normalized moment `ξ_n`.
""" HCCGeneratingFunction

@doc """
    HCCLocalFactor

Prime-field local Euler factor for an HCC value.  Coefficients are stored in
ascending powers of `T`; `[1, p-a]` represents `1 - aT` over `F_p`.
""" HCCLocalFactor

@doc """
    FREquationSystem

Backend-neutral F/R equation system for a multiplicity-free fusion rule.

# Fields
- `rules`: Fusion rule.
- `variables`: Symbolic variables recorded by the backend-neutral layer.
- `equations`: Pentagon/hexagon equations.
- `metadata`: Construction metadata.
- `assumptions`: Human-readable assumptions.
- `base_ring`: Symbolic base-ring label.

# Notes
Stable as an equation-container boundary.  General solver backends and
finite-field reconstruction are experimental.
""" FREquationSystem

@doc """
    FRData

Exact F/R data container with accessor-based storage.

# Fields
- `rules`: Fusion rule.
- `F_values`: F-symbol coordinates in TensorCategories order.
- `R_values`: Forward R-symbol coordinates.
- `R_inverse_values`: Inverse R-symbol coordinates.
- `metadata`: Object labels and optional extended symbol tables.

# Notes
Stable accessor API.  Raw vector ordering is a compatibility detail; prefer
`F_symbol`, `R_symbol`, `F_values`, and `R_values`.
""" FRData

@doc """
    BraidRepresentation

Braid group representation matrices built from exact `FRData`.

# Fields
- `fr_data`: Exact F/R data.
- `objects`: Boundary object labels.
- `total`: Total charge label.
- `basis`: Fusion-tree basis.
- `generators`: Matrices for the standard braid generators.
- `metadata`: Construction metadata.

# Notes
Stable for multiplicity-free built-in examples.  Finite-field image and
Zariski diagnostics built on top of these matrices are experimental.
""" BraidRepresentation

@doc """
    ZariskiClosureDiagnostics

Experimental diagnostic summary for finite-field Zariski-style evidence.

# Notes
This record reports computational evidence from finite-field calculations.  It
is not a proof of characteristic-zero Zariski closure, density, or exact
classification.
""" ZariskiClosureDiagnostics

@doc """
    FRStatus

Enum describing exact F/R reconstruction status in `ClassifiedMTC`.
""" FRStatus

@doc "Exact F/R reconstruction was intentionally skipped." FRSkipped
@doc "Exact F/R reconstruction found and verified F/R data." FRSolved
@doc "Exact F/R reconstruction found no solution in the attempted search." FRNoSolutionFound
@doc "Exact F/R reconstruction stopped because of timeout-like resource limits." FRTimeoutLikeFailure
@doc "Exact F/R reconstruction failed." FRReconstructionFailed
@doc "F/R data were reconstructed but exact verification failed." FRVerificationFailed

_brief_docs = Dict{Symbol, String}(
    :field => "    field(x)\n\nReturn the cyclotomic field parent associated with a context or exact data object.",
    :zeta => "    zeta(ctx)\n\nReturn the distinguished primitive root of unity for a cyclotomic context.",
    :conductor => "    conductor(x)\n\nReturn the conductor attached to a cyclotomic context or exact modular data.",
    :cond_S => "    cond_S(data)\n\nReturn the declared conductor of the `S`-matrix entries.",
    :cond_T => "    cond_T(data)\n\nReturn the declared conductor of the `T`-matrix entries.",
    :cond_F => "    cond_F(data)\n\nReturn the declared conductor of the F/R layer when known, or `nothing`.",
    :quantum_dimensions => "    quantum_dimensions(data)\n\nReturn quantum dimensions extracted from exact modular data.",
    :higher_central_charge => "    higher_central_charge(data, n = 1; normalization = :D2)\n\nReturn the D²-normalized higher central charge moment `ξ_n`.",
    :higher_central_charges => "    higher_central_charges(data, ns)\n\nReturn D²-normalized higher central charge moments for the supplied integers.",
    :higher_central_charge_result => "    higher_central_charge_result(data; n = 1, normalization = :galois)\n\nReturn the legacy structured Gauss-sum result used by compatibility code.",
    :higher_central_charge_period => "    higher_central_charge_period(data)\n\nReturn the period used by the HCC sequence, usually the twist conductor.",
    :higher_central_charge_sequence => "    higher_central_charge_sequence(data)\n\nReturn one full period of D²-normalized HCC values.",
    :higher_central_charge_generating_function => "    higher_central_charge_generating_function(data)\n\nReturn exact weights and twists representing the HCC generating function.",
    :higher_central_charge_modp => "    higher_central_charge_modp(data, n, p; embedding = nothing)\n\nReduce a D²-normalized HCC value modulo a good prime.",
    :higher_central_charge_sequence_modp => "    higher_central_charge_sequence_modp(data, p; embedding = nothing)\n\nReturn one full HCC period reduced modulo `p`.",
    :hcc_local_factor => "    hcc_local_factor(data, n, p; embedding = nothing)\n\nReturn the prime-field local factor `1 - ξ_{n,p}T`.",
    :hcc_local_factors => "    hcc_local_factors(data, p; embedding = nothing)\n\nReturn one full period of prime-field HCC local factors.",
    :enumerate_strata => "    enumerate_strata(catalog, r; require_unit_summand = false, max_multiplicity = typemax(Int))\n\nEnumerate semisimple representation strata of total dimension `r` from an atomic `SL(2, ℤ/N)` catalog.",
    :count_strata => "    count_strata(catalog, r; kwargs...)\n\nCount rank-`r` strata for an atomic catalog.",
    :describe_stratum => "    describe_stratum(s, catalog)\n\nReturn a compact direct-sum description of a stratum using atomic-catalog labels.",
    :classify_mtcs_at_conductor => "    classify_mtcs_at_conductor(N; kwargs...)\n\nSearch for exact cyclotomic modular-data candidates whose modular representation factors through `SL(2, ℤ/N)`.",
    :classify_mtcs_auto => "    classify_mtcs_auto(N; kwargs...)\n\nConvenience wrapper around `classify_mtcs_at_conductor` that selects rank cutoffs and admissible primes.",
    :recommend_primes => "    recommend_primes(N; kwargs...)\n\nSuggest admissible finite-field primes for conductor-level searches.",
    :recommend_skip_FR => "    recommend_skip_FR(N, rank; kwargs...)\n\nHeuristic recommendation for whether exact F/R solving is likely to be expensive.",
    :galois_orbit => "    galois_orbit(data)\n\nReturn Galois orbit data using the conductor attached to `data`.",
    :frobenius => "    frobenius(x, p)\n\nApply the Frobenius operation at the prime `p` to a supported exact object.",
    :reduce_mod_p => "    reduce_mod_p(x, p)\n\nReduce a supported ACMG object modulo the prime `p`.",
    :EquationVariable => "    EquationVariable\n\nSymbolic variable descriptor used by backend-neutral F/R equation systems.",
    :EquationTerm => "    EquationTerm\n\nSparse monomial term in an `EquationExpr`.",
    :EquationExpr => "    EquationExpr\n\nSparse polynomial expression used by ACMG's backend-neutral equation layer.",
    :PolynomialEquation => "    PolynomialEquation\n\nEquation with left/right `EquationExpr` sides and metadata.",
    :EquationSystem => "    EquationSystem\n\nGeneric backend-neutral equation-system container.",
    :GaugeVariable => "    GaugeVariable\n\nGauge variable descriptor for a multiplicity-free fusion channel `(a,b,c)`.",
    :FiniteFieldEquationSystem => "    FiniteFieldEquationSystem\n\nReduction of an `FREquationSystem` modulo a prime `p`.",
    :fusion_channels => "    fusion_channels(rules, a, b)\n\nReturn simple objects `c` for which `N_ab^c` is nonzero.",
    :is_admissible => "    is_admissible(rules, a, b, c)\n\nReturn whether the fusion channel `a ⊗ b -> c` is present.",
    :is_multiplicity_free => "    is_multiplicity_free(rules)\n\nReturn whether every fusion multiplicity is either zero or one.",
    :require_multiplicity_free => "    require_multiplicity_free(rules)\n\nValidate that a fusion rule is multiplicity-free, returning `true` or throwing.",
    :hexagon_equations => "    hexagon_equations(rules; context = nothing)\n\nGenerate TensorCategories-backed hexagon equations for a multiplicity-free fusion rule.",
    :validate_fr_system => "    validate_fr_system(system)\n\nValidate structural assumptions of an `FREquationSystem`.",
    :gauge_variables => "    gauge_variables(rules)\n\nReturn channel gauge variables for a multiplicity-free fusion rule.",
    :gauge_fix => "    gauge_fix(system; strategy = :toric_snf)\n\nFix a deterministic F-only toric Smith-normal-form slice by substituting selected F-symbol coordinates to `1`; use strategy=:none to disable equation-level gauge fixing.",
    :toric_gauge_data => "    toric_gauge_data(frdata; include_R = true)\n    toric_gauge_data(fusion; field = nothing, conventions = :tensorcategories)\n\nReturn toric gauge coordinates, character matrix, and Smith-normal-form split for multiplicity-free FRData, or pre-reconstruction toric gauge data for a multiplicity-free fusion rule.",
    :toric_gauge_normal_form => "    toric_gauge_normal_form(frdata; include_R = true)\n\nReport the pre-solver toric F-slice and Smith-normal-form stabilizer metadata recorded on solved FRData.",
    :frobenius_metadata => "    frobenius_metadata(p, conductor)\n\nReturn metadata describing a prime relative to a conductor.",
    :check_modular_data => "    check_modular_data(candidate, known_data)\n\nCompare a candidate modular-data object with known data.",
    :semion_fusion_rules => "    semion_fusion_rules()\n\nReturn the multiplicity-free fusion rule for the semion category.",
    :fibonacci_fusion_rules => "    fibonacci_fusion_rules()\n\nReturn the multiplicity-free Fibonacci fusion rule.",
    :toric_code_fusion_rules => "    toric_code_fusion_rules()\n\nReturn the multiplicity-free toric-code fusion rule.",
    :ising_fusion_rules => "    ising_fusion_rules()\n\nReturn the multiplicity-free Ising fusion rule.",
    :fusion_rule => "    fusion_rule(data)\n\nReturn the `FusionRule` stored in an `FRData` object.",
    :F_values => "    F_values(data)\n\nReturn vector-backed F-symbol coordinates in ACMG/TensorCategories order.",
    :R_values => "    R_values(data)\n\nReturn the forward R-symbol coordinate vector.",
    :R_inverse_values => "    R_inverse_values(data)\n\nReturn inverse R-symbol coordinates.",
    :fr_metadata => "    fr_metadata(data)\n\nReturn the metadata dictionary attached to `FRData`.",
    :fr_scalar_type => "    fr_scalar_type(data)\n\nReturn the scalar element type used by an `FRData` object.",
    :fr_value_one => "    fr_value_one(data)\n\nReturn a multiplicative identity compatible with the scalar type of `data`.",
    :has_F_symbol => "    has_F_symbol(data, args...; kwargs...)\n\nReturn whether `F_symbol` can be evaluated for the requested indices.",
    :has_R_symbol => "    has_R_symbol(data, args...; kwargs...)\n\nReturn whether `R_symbol` can be evaluated for the requested indices.",
    :gauge_basis_indices => "    gauge_basis_indices(data, a, b, c)\n\nReturn Hom-basis indices used by scalar gauge transforms.",
    :frdata_from_vectors => "    frdata_from_vectors(rules, Fvals, Rvals; metadata = Dict())\n\nConstruct `FRData` from vector-backed F/R coordinates.",
    :frdata_from_namedtuple => "    frdata_from_namedtuple(nt)\n\nConstruct `FRData` from a named-tuple payload produced by ACMG helpers.",
    :fr_pentagon_values => "    fr_pentagon_values(data)\n\nReturn F-symbol values used by pentagon-equation infrastructure.",
    :fr_hexagon_values => "    fr_hexagon_values(data)\n\nReturn F/R values used by hexagon-equation infrastructure.",
    :GaugeChoice => "    GaugeChoice\n\nAlias for `GaugeParameters`.",
    :FusionPath => "    FusionPath\n\nOne basis path in a left-associated multiplicity-free fusion tree.",
    :FusionTreeBasis => "    FusionTreeBasis\n\nBasis of fusion paths for a sequence of objects and a total charge.",
    :fusion_paths => "    fusion_paths(rules, objects, total)\n\nEnumerate multiplicity-free fusion paths from `objects` to `total`.",
    :fusion_basis => "    fusion_basis(rules, objects, total)\n\nConstruct a `FusionTreeBasis` for braid-representation computations.",
    :dim => "    dim(basis)\n\nReturn the dimension of a fusion-tree basis.",
    :braid_representation => "    braid_representation(fr_data, objects, total)\n\nConstruct braid generators on a fusion space.",
    :braid_generator => "    braid_generator(br, i)\n\nReturn braid generator `σ_i` from a `BraidRepresentation`.",
    :braid_generators => "    braid_generators(br)\n\nReturn all braid generators stored in a `BraidRepresentation`.",
    :check_braid_relations => "    check_braid_relations(br)\n\nCheck braid relations for the generators in `br`.",
    :FiniteFieldBraidRepresentation => "    FiniteFieldBraidRepresentation\n\nExperimental finite-field reduction of a braid representation.",
    :MatrixAlgebraDiagnostics => "    MatrixAlgebraDiagnostics\n\nExperimental finite-field matrix-algebra diagnostic record.",
    :CommutantDiagnostics => "    CommutantDiagnostics\n\nExperimental finite-field commutant diagnostic record.",
    :reconstruct_cyclotomic_element_from_residues => "    reconstruct_cyclotomic_element_from_residues(residues, primes, N)\n\nExperimental CRT helper for reconstructing a cyclotomic element.",
    :assign_F_to_associator! => "    assign_F_to_associator!(args...)\n\nExperimental/internal helper for associator coordinates.",
    :get_hexagon_system => "    get_hexagon_system(args...)\n\nLow-level hexagon-equation system constructor.",
    :get_hexagon_fr_system => "    get_hexagon_fr_system(args...)\n\nLow-level F/R hexagon-equation system constructor.",
    :number_of_variables_in_hexagon_equations => "    number_of_variables_in_hexagon_equations(args...)\n\nReturn the number of variables used by the low-level hexagon system.",
    :solve_pentagon_modular_crt => "    solve_pentagon_modular_crt(args...)\n\nExperimental pentagon solver using modular CRT reconstruction.",
    :solve_pentagon_homotopy => "    solve_pentagon_homotopy(args...)\n\nExperimental pentagon solver hook.",
    :solve_pentagon_newton => "    solve_pentagon_newton(args...)\n\nExperimental pentagon Newton solver hook.",
    :solve_hexagon_modular_crt => "    solve_hexagon_modular_crt(args...)\n\nExperimental hexagon solver using modular CRT reconstruction.",
    :solve_hexagon_homotopy => "    solve_hexagon_homotopy(args...)\n\nExperimental hexagon solver hook.",
    :fr_status => "    fr_status(m::ClassifiedMTC)\n\nReturn the exact F/R reconstruction status for a classification result.",
    :compute_FR_from_ST => "    compute_FR_from_ST(args...)\n\nCompute exact F/R data from exact modular data when ACMG has a supported route.",
    :is_modular_data_automorphism => "    is_modular_data_automorphism(data, perm)\n\nReturn whether `perm` preserves the exact modular data.",
)

for (sym, text) in _brief_docs
    @eval @doc $text $(Symbol(sym))
end

nothing
