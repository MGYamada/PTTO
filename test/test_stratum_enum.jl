using Test
using ACMG

"""
Tests for StratumEnum module.

We use synthetic (hand-built) AtomicIrrep catalogs so these tests don't
require Oscar/GAP at all.
"""

# Helper: build a minimal AtomicIrrep with fake S, T for testing
# (we only care about dim, level, parity for stratum enumeration)
function fake_atomic(dim::Int, level::Int, parity::Int, label::String)
    return AtomicIrrep(dim, level, label, nothing, nothing, parity, nothing, level)
end

@testset "enumerate_strata" begin

    @testset "Rank 1: trivial only" begin
        catalog = [fake_atomic(1, 1, +1, "1d_1")]
        strata = enumerate_strata(catalog, 1)
        @test length(strata) == 1
        @test strata[1].multiplicities == Dict(1 => 1)
        @test strata[1].total_dim == 1
    end

    @testset "Rank 2: trivial + one 1d irrep" begin
        catalog = [
            fake_atomic(1, 1, +1, "1d_1"),
            fake_atomic(1, 2, +1, "1d_2"),
        ]
        # By default, no unit summand requirement → all partitions of 2 in dim-1 blocks
        # (2,0), (1,1), (0,2) → 3 strata
        strata = enumerate_strata(catalog, 2)
        @test length(strata) == 3

        # With require_unit_summand: (0,2) excluded
        strata_req = enumerate_strata(catalog, 2; require_unit_summand = true)
        @test length(strata_req) == 2
        @test !(Dict(2 => 2) in [s.multiplicities for s in strata_req])
    end

    @testset "Rank 3 with 2d irrep" begin
        catalog = [
            fake_atomic(1, 1, +1, "1d_1"),
            fake_atomic(2, 5, -1, "2d_5"),
        ]
        # Partitions of 3 = 3×1 or 1 + 2
        strata = enumerate_strata(catalog, 3)
        @test length(strata) == 2
        @test Dict(1 => 3) in [s.multiplicities for s in strata]
        @test Dict(1 => 1, 2 => 1) in [s.multiplicities for s in strata]
    end

    @testset "Rank 5 (3,2)-type decomposition" begin
        # Synthetic catalog mimicking what N=24 provides for the (3,2)-type case
        catalog = [
            fake_atomic(1, 1, +1, "1d_1"),    # [1] trivial
            fake_atomic(1, 3, +1, "1d_3"),    # [2] non-trivial 1d at level 3
            fake_atomic(2, 3, -1, "2d_3"),    # [3] 2d at level 3 (ρ_3^(2))
            fake_atomic(3, 8, +1, "3d_8"),    # [4] 3d at level 8 (ρ_8^(3))
        ]
        strata = enumerate_strata(catalog, 5)
        # Constraint: m1 + m2 + 2*m3 + 3*m4 = 5
        # m4=0:  m3=0: m1+m2=5: 6 ways (m1=0..5) → 6
        #        m3=1: m1+m2=3: 4 ways → 4
        #        m3=2: m1+m2=1: 2 ways → 2
        # m4=1:  m3=0: m1+m2=2: 3 ways → 3
        #        m3=1: m1+m2=0: 1 way → 1
        # Total: 6+4+2+3+1 = 16
        @test length(strata) == 16

        # The key (3,2)-type stratum: 3d_8 + 2d_3 = 3 + 2 = 5, no explicit 1d_1
        # (m1=0, m2=0, m3=1, m4=1)
        @test Dict(3 => 1, 4 => 1) in [s.multiplicities for s in strata]
    end

    @testset "describe_stratum" begin
        catalog = [
            fake_atomic(1, 1, +1, "1d_1"),
            fake_atomic(2, 5, -1, "2d_5"),
        ]
        s = Stratum(Dict(1 => 2, 2 => 1), 4)
        desc = describe_stratum(s, catalog)
        @test occursin("1d_1", desc)
        @test occursin("2d_5", desc)
        @test occursin("⊕", desc)
    end

    @testset "find_unit_indices" begin
        catalog = [
            fake_atomic(1, 1, +1, "1d_1"),
            fake_atomic(1, 1, -1, "1d_1"),   # same label but parity -1
            fake_atomic(2, 3, +1, "2d_3"),
        ]
        unit_idx = find_unit_indices(catalog)
        @test unit_idx == [1]   # only the parity=+1 one
    end
end

