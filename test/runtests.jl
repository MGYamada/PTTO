using Test
using ACMG

@testset "ACMG Prototype" begin
    include("test_fparith.jl")
    include("test_fibonacci.jl")
    include("test_ising.jl")
    include("test_sl2reps.jl")
    include("test_stratum_enum.jl")
end
