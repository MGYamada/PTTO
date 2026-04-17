using Test
import ACMG
import ACMG: build_modular_datum, validate_modular_data, extract_fusion_rule_Fp, extract_and_lift
import ACMG: sqrt_mod, root_of_unity, matmul_mod, matpow_mod, diagmul_right
import ACMG: fusion_isomorphic, FusionRule

@testset "Fibonacci at p=41" begin
    # --- Setup from hand calculation ---
    p = 41
    N = 5        # conductor
    r = 2        # rank

    # ζ_5 in F_41
    ζ5 = root_of_unity(5, p)
    @test ζ5 == 10

    # √5 in F_41
    sqrt5 = sqrt_mod(5, p)
    @test sqrt5 == 13 || sqrt5 == 41 - 13  # either root

    # φ = (1 + √5) / 2
    # Using √5 = 13: φ = 14 / 2 = 7. Using √5 = 28: φ = 29 / 2 = 29*21 mod 41
    # 2⁻¹ mod 41 = 21 (since 2*21 = 42 = 1)
    inv2 = invmod(2, p)
    @test inv2 == 21
    φ = mod((1 + 13) * inv2, p)
    @test φ == 7
    # Verify φ² = φ + 1
    @test mod(φ^2, p) == mod(φ + 1, p)
    @test mod(φ^2, p) == 8

    # D² = 1 + φ² = 9, D = 3
    D² = mod(1 + φ^2, p)
    @test D² == 9
    D = sqrt_mod(D², p)
    @test D == 3 || D == 38  # two roots

    # --- Construct Fibonacci S-matrix ---
    # S = (1/D) * [1  φ ; φ  -1]
    # With D=3, 1/D = 1/3 = 14 mod 41 (since 3*14 = 42 = 1)
    invD = invmod(3, p)  # use D = 3
    @test invD == 14
    S = [mod(invD * 1, p)   mod(invD * φ, p);
         mod(invD * φ, p)   mod(invD * (p-1), p)]
    @test S == [14 16; 16 27]

    # T = diag(1, ζ_5²) = diag(1, 18)  (since 10² = 100 = 18 mod 41)
    T = [1, mod(ζ5^2, p)]
    @test T == [1, 18]

    # --- Validate modular data ---
    result = validate_modular_data(S, T, p, N)
    @test result.valid
    @test result.α == 21       # hand calculation: (ST)³ = 21 I
    @test result.D_squared == 9
    @test result.D == 3 || result.D == 38
    @test result.C == [1, 2]   # all self-dual (S² = I)

    # α should be ζ_{20}^7 in F_41 (central charge c = 14/5)
    ζ20 = root_of_unity(20, p)
    @test ζ20 == 36
    @test powermod(ζ20, 7, p) == 21  # α = ζ_{20}^7

    # --- Build ModularDatumFp ---
    md = build_modular_datum(S, T, p, N)
    @test md !== nothing
    @test md.rank == 2
    @test md.α == 21

    # --- Verlinde extraction ---
    N_Fp = extract_fusion_rule_Fp(md)
    # Fibonacci: τ ⊗ τ = 1 + τ
    # In 1-indexed: N[2, 2, 1] = 1, N[2, 2, 2] = 1
    # And N[1, i, j] = δ_{ij}, N[2, 1, 2] = 1
    @test N_Fp[1, 1, 1] == 1    # 1 ⊗ 1 = 1
    @test N_Fp[1, 1, 2] == 0
    @test N_Fp[1, 2, 1] == 0
    @test N_Fp[1, 2, 2] == 1    # 1 ⊗ τ = τ
    @test N_Fp[2, 1, 2] == 1
    @test N_Fp[2, 2, 1] == 1    # τ ⊗ τ ⊇ 1
    @test N_Fp[2, 2, 2] == 1    # τ ⊗ τ ⊇ τ

    # --- Lift to Z and build FusionRule ---
    fr = extract_and_lift(md)
    @test fr !== nothing
    @test fr.rank == 2
    @test fr.N[2, 2, 1] == 1
    @test fr.N[2, 2, 2] == 1
    @test fr.dual == [1, 2]  # all self-dual

    # --- Isomorphism with reference Fibonacci ---
    fib_reference_N = zeros(Int, 2, 2, 2)
    fib_reference_N[1, 1, 1] = 1
    fib_reference_N[1, 2, 2] = 1
    fib_reference_N[2, 1, 2] = 1
    fib_reference_N[2, 2, 1] = 1
    fib_reference_N[2, 2, 2] = 1
    fib_reference = FusionRule(fib_reference_N)
    @test fusion_isomorphic(fr, fib_reference)
end
