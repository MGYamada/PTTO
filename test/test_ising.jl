using Test
import ACMG
import ACMG: build_modular_datum, validate_modular_data, extract_fusion_rule_Fp, extract_and_lift
import ACMG: sqrt_mod, root_of_unity, matmul_mod, matpow_mod, diagmul_right
import ACMG: fusion_isomorphic, FusionRule

@testset "Ising at p=17" begin
    # --- Setup from hand calculation ---
    p = 17
    N = 16
    r = 3

    # ζ_16 in F_17
    ζ16 = root_of_unity(16, p)
    @test ζ16 == 3

    # √2 = ζ_8 + ζ_8^{-1} in F_17
    sqrt2 = sqrt_mod(2, p)
    @test sqrt2 == 6 || sqrt2 == 11  # 11 is what hand calculation found

    # --- Construct Ising S-matrix ---
    # S = (1/2) * [1  √2  1; √2 0 -√2; 1 -√2 1]
    # Using √2 = 11, 1/2 = 9 (since 2*9 = 18 = 1 mod 17)
    inv2 = invmod(2, p)
    @test inv2 == 9

    # Use sqrt2 = 11 to match hand calculation
    s2 = 11
    S = [mod(inv2 * 1, p)    mod(inv2 * s2, p)       mod(inv2 * 1, p);
         mod(inv2 * s2, p)   0                        mod(inv2 * (p - s2), p);
         mod(inv2 * 1, p)    mod(inv2 * (p - s2), p)  mod(inv2 * 1, p)]
    @test S == [9 14 9; 14 0 3; 9 3 9]

    # T = diag(1, ζ_16, -1) = diag(1, 3, 16)
    T = [1, ζ16, p - 1]
    @test T == [1, 3, 16]

    # --- Validate ---
    result = validate_modular_data(S, T, p, N)
    @test result.valid
    @test result.α == 3       # hand calc: (ST)³ = 3 I, and α = ζ_16
    @test result.D_squared == 4
    @test result.D == 2 || result.D == 15  # two square roots of 4
    @test result.C == [1, 2, 3]  # all self-dual

    # α = ζ_16 check
    @test result.α == ζ16

    # --- Verlinde ---
    md = build_modular_datum(S, T, p, N)
    @test md !== nothing

    N_Fp = extract_fusion_rule_Fp(md)
    # Ising:
    #   1 ⊗ X = X for all X
    #   σ ⊗ σ = 1 + ψ
    #   σ ⊗ ψ = σ
    #   ψ ⊗ ψ = 1
    # 1-indexed: 1 = unit, 2 = σ, 3 = ψ
    @test N_Fp[2, 2, 1] == 1   # σ ⊗ σ ⊇ 1
    @test N_Fp[2, 2, 2] == 0   # σ ⊗ σ ⊉ σ
    @test N_Fp[2, 2, 3] == 1   # σ ⊗ σ ⊇ ψ
    @test N_Fp[2, 3, 2] == 1   # σ ⊗ ψ = σ
    @test N_Fp[2, 3, 1] == 0
    @test N_Fp[2, 3, 3] == 0
    @test N_Fp[3, 3, 1] == 1   # ψ ⊗ ψ = 1
    @test N_Fp[3, 3, 2] == 0
    @test N_Fp[3, 3, 3] == 0

    # --- Build FusionRule ---
    fr = extract_and_lift(md)
    @test fr !== nothing
    @test fr.rank == 3
    @test fr.dual == [1, 2, 3]

    # --- Isomorphism with reference Ising ---
    ising_ref = zeros(Int, 3, 3, 3)
    ising_ref[1, 1, 1] = 1; ising_ref[1, 2, 2] = 1; ising_ref[1, 3, 3] = 1
    ising_ref[2, 1, 2] = 1; ising_ref[2, 2, 1] = 1; ising_ref[2, 2, 3] = 1; ising_ref[2, 3, 2] = 1
    ising_ref[3, 1, 3] = 1; ising_ref[3, 2, 2] = 1; ising_ref[3, 3, 1] = 1
    ising_reference = FusionRule(ising_ref)
    @test fusion_isomorphic(fr, ising_reference)
end

@testset "Ising with different √2 choice" begin
    # Use √2 = 6 (the other square root) — should still give same fusion rule
    p = 17
    N = 16
    inv2 = invmod(2, p)
    s2 = 6  # other choice
    S = [mod(inv2 * 1, p)    mod(inv2 * s2, p)       mod(inv2 * 1, p);
         mod(inv2 * s2, p)   0                        mod(inv2 * (p - s2), p);
         mod(inv2 * 1, p)    mod(inv2 * (p - s2), p)  mod(inv2 * 1, p)]
    T = [1, 3, 16]

    result = validate_modular_data(S, T, p, N)
    @test result.valid

    md = build_modular_datum(S, T, p, N)
    @test md !== nothing
    fr = extract_and_lift(md)
    @test fr !== nothing
    # Same fusion rule regardless of √2 choice
    @test fr.N[2, 2, 1] == 1
    @test fr.N[2, 2, 3] == 1
end
