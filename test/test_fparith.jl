using Test
import ACMG: is_square, sqrt_mod, primitive_root, root_of_unity, roots_of_unity
import ACMG: matmul_mod, matpow_mod, diagmul_right, diagmul_left, lift_symmetric

@testset "FpArith primitives" begin
    # Euler criterion
    @test is_square(1, 17) == true
    @test is_square(0, 17) == true
    @test is_square(2, 17) == true   # 6² = 36 ≡ 2 mod 17
    @test is_square(3, 17) == false
    @test is_square(4, 41) == true
    @test is_square(9, 41) == true   # computed in Fibonacci
    @test is_square(21, 31) == false # computed in Fibonacci/p=31 case

    # Square roots
    @test mod(sqrt_mod(2, 17)^2, 17) == 2
    @test mod(sqrt_mod(9, 41)^2, 41) == 9
    @test mod(sqrt_mod(5, 41)^2, 41) == 5   # Fibonacci: √5 = 13 mod 41
    @test mod(sqrt_mod(5, 11)^2, 11) == 5   # Fibonacci: √5 = 4 mod 11

    # Primitive root
    @test primitive_root(11) == 2
    @test primitive_root(17) == 3
    @test primitive_root(41) == 6

    # Roots of unity
    @test root_of_unity(5, 11) == 4           # Fibonacci test
    @test root_of_unity(5, 41) == 10          # Fibonacci test (p=41)
    @test root_of_unity(16, 17) == 3          # Ising test

    # ζ_N has order exactly N
    for (n, p) in [(5, 11), (5, 41), (16, 17), (8, 17)]
        ζ = root_of_unity(n, p)
        @test powermod(ζ, n, p) == 1
        for k in 1:(n-1)
            @test powermod(ζ, k, p) != 1
        end
    end

    # All n-th roots
    all5_at_41 = roots_of_unity(5, 41)
    @test length(all5_at_41) == 5
    @test all5_at_41[1] == 1
    @test sort(all5_at_41) == sort([1, 10, 100%41, powermod(10, 3, 41), powermod(10, 4, 41)])

    # Matrix operations
    A = [1 2; 3 4]
    B = [5 6; 7 8]
    AB_expected = mod.(A * B, 17)
    @test matmul_mod(A, B, 17) == AB_expected
    @test matpow_mod(A, 3, 17) == mod.(A^3, 17)
    @test matpow_mod(A, 0, 17) == Matrix{Int}(I, 2, 2)

    @test diagmul_right([1 2 3; 4 5 6], [2, 3, 5], 17) == [2 6 15; 8 15 30%17]
    @test diagmul_left([2, 3], [1 2; 3 4], 17) == [2 4; 9 12]

    # Symmetric lift
    @test lift_symmetric(3, 17) == 3
    @test lift_symmetric(15, 17) == -2
    @test lift_symmetric(8, 17) == 8
    @test lift_symmetric(9, 17) == -8
end
