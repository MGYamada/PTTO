if !isdefined(@__MODULE__, :test_ising_fr_data_mod_p_17)
    function test_ising_fr_data_mod_p_17()
        p = 17
        fp(xs) = FpElem[FpElem(x, p) for x in xs]
        F = fp([1, 12, 16, 11, 4, 7, 4, 13, 5, 1, 12, 2, 1, 14])
        R = fp([16, 1, 1, 1, 10, 6, 1, 1, 1, 1])
        Rinv = fp([16, 1, 1, 1, 12, 3, 1, 1, 1, 1])
        return FRData(ising_fusion_rules(), F, R, Rinv;
                      metadata = Dict{Symbol, Any}(:solver_status => :fixture,
                                                   :p => p,
                                                   :base_field => :F_17,
                                                   :name => :ising,
                                                   :source => :precomputed_test_fixture))
    end
end
