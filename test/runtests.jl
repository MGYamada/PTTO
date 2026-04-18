using Test
using ACMG

@testset "ACMG Prototype" begin
    include("test_fparith.jl")
    include("test_fibonacci.jl")
    include("test_ising.jl")
    include("test_sl2reps.jl")
    include("test_stratum_enum.jl")
    include("test_blocku.jl")
    include("test_crt.jl")
    include("test_block_u_general.jl")
    include("test_phase4_fibonacci.jl")
    include("test_phase4_lift.jl")
    include("test_phase4_verify.jl")
    include("test_kitaev_complex.jl")
end
