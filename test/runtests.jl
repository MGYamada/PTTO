using Test

@testset "ACMG Prototype" begin
    include("test_fparith.jl")
    include("test_ising.jl")
    include("test_sl2reps.jl")
    include("test_blocku.jl")
    include("test_stratum_enum.jl")
    include("test_block_u_general.jl")
    include("test_crt.jl")
    include("test_modular_data_lift.jl")
    include("test_cyclotomic_context.jl")
    include("test_phase4_roundtrip.jl")
    include("test_pipeline_primes.jl")
    include("test_pipeline_galois_grouping.jl")
    include("test_pipeline_end_to_end.jl")
    include("test_pipeline_rank23_smoke.jl")
end
