# Documentation example: F/R accessors, gauge metadata, and exact braid matrices.
#
# This example stays inside the built-in semion data so it remains small.  It
# demonstrates the stable F/R accessor boundary and exact braid relation check.

using ACMG

function main()
    rules = semion_fusion_rules()
    system = fr_equation_system(rules)

    @assert is_multiplicity_free(rules)
    @assert validate_fr_system(system)
    @assert system.metadata[:rank] == 2

    fixed = gauge_fix(system; strategy = :safe)
    @assert validate_fr_system(fixed)
    @assert fixed.metadata[:gauge_fix_strategy] == :safe

    data = semion_fr_data()
    @assert validate_frdata_for_gauge(data)
    @assert fusion_coeff(data, :s, :s, :one) == 1
    @assert has_R_symbol(data, :s, :s, :one)

    plan = gauge_fixing_plan(data)
    @assert plan isa Vector

    basis = fusion_basis(fusion_rule(data), [:s, :s, :s], :s)
    @assert dim(basis) == 1

    br = braid_representation(data, [:s, :s, :s], :s)
    @assert length(braid_generators(br)) == 2
    @assert check_braid_relations(br).ok

    println("semion F/R equations: ", length(system.equations))
    println("fusion basis dimension: ", dim(basis))
    println("braid generators: ", length(braid_generators(br)))
end

main()
