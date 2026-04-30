# Documentation example: lightweight export helpers.
#
# The IO layer stores exact cyclotomic values in textual JSON payloads.  The
# loader returns a Dict so archived data can be inspected without pretending to
# reconstruct arbitrary Oscar parents from JSON.

using ACMG

function main()
    data = semion_modular_data()
    modular_payload = export_modular_data(data; format = :dict)

    @assert modular_payload["type"] == "ModularData"
    @assert modular_payload["N"] == 8
    @assert modular_payload["cond_T"] == 4

    rules = semion_fusion_rules()
    fusion_payload = export_fusion_rule(rules; format = :dict)

    @assert fusion_payload["type"] == "FusionRule"
    @assert fusion_payload["rank"] == 2

    out = joinpath(mktempdir(), "semion_report.md")
    result = classify_mtcs_at_conductor(8;
        max_rank = 2,
        primes = [17, 41],
        skip_FR = true,
        verbose = false,
    )

    @assert all(m -> fr_status(m) == FRSkipped, result)
    write_report(out, result)
    @assert isfile(out)

    println("exported modular-data keys: ", sort(collect(keys(modular_payload))))
    println("classified records: ", length(result))
    println("report path: ", out)
end

main()
