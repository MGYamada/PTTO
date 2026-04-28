using Test
using ACMG

function _io_semion_fusion()
    N = zeros(Int, 2, 2, 2)
    N[1, 1, 1] = 1
    N[1, 2, 2] = 1
    N[2, 1, 2] = 1
    N[2, 2, 1] = 1
    return N
end

@testset "IO exports" begin
    data = semion_modular_data()
    Nijk = _io_semion_fusion()
    mtc = ACMG.ClassifiedMTC(8, 8, 2, ACMG.Stratum(Dict(1 => 2), 2), Nijk,
                             2, [17, 41], Int[], true, nothing,
                             data.S, data.T, nothing, nothing, nothing, 1)

    dir = mktempdir()
    json_path = joinpath(dir, "semion.json")
    md_path = joinpath(dir, "semion.md")

    @test save_classification(json_path, mtc) == json_path
    loaded = load_classification(json_path)
    @test loaded["schema"] == "ACMG.classification"
    @test loaded["count"] == 1
    @test loaded["items"][1]["N"] == 8
    @test loaded["items"][1]["rank"] == 2
    @test loaded["items"][1]["F_values"] === nothing

    modular_payload = export_modular_data(mtc; format = :dict)
    @test modular_payload["N"] == 8
    @test length(modular_payload["S"]) == 2

    fusion_payload = export_fusion_rule(mtc; format = :dict)
    @test fusion_payload["rank"] == 2
    @test fusion_payload["Nijk"][2][2][1] == 1

    fr_payload = export_FR(mtc; format = :dict)
    @test fr_payload["F_values"] === nothing
    @test fr_payload["verify_report"] === nothing

    @test write_report(md_path, [mtc]) == md_path
    report = read(md_path, String)
    @test occursin("# ACMG Classification Report", report)
    @test occursin("F/R attached: false", report)
end
