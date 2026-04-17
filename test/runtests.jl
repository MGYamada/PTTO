using Test
using ACMG

@testset "ACMG Prototype" begin
    include("test_fparith.jl")
    include("test_fibonacci.jl")
    include("test_ising.jl")
end
