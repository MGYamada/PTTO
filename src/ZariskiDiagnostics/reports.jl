struct ZariskiClosureDiagnostics
    exact_closure_computed::Bool
    diagnostic_only::Bool
    representation_dimension::Int
    field_characteristic::Union{Int, Nothing}
    matrix_algebra_dimension::Int
    commutant_dimension::Int
    appears_irreducible::Bool
    appears_full_matrix_algebra::Bool
    determinant_constraints::Vector{Int}
    projective_orders::Vector{Union{Int, Nothing}}
    notes::Vector{String}
end

"""
    zariski_closure_diagnostics(br; max_words = 200, max_degree = 2)

Experimental API.

Collect finite-field diagnostics that are suggestive of the braid image's
Zariski behavior.  The input is a `FiniteFieldBraidRepresentation`; the output
is a `ZariskiClosureDiagnostics` record containing matrix-algebra,
commutant, determinant, projective-order, and note fields.

Mathematical caveats: this function does not compute a full Zariski closure
or classify algebraic groups.  Treat diagnostics as computational evidence,
not as mathematical theorems.  API inputs and outputs may change.
"""
function zariski_closure_diagnostics(br::FiniteFieldBraidRepresentation; max_words = 200, max_degree = 2)
    warn_experimental("zariski_closure_diagnostics")
    alg = generated_matrix_algebra(br)
    comm = commutant(br)
    gdiag = finite_group_diagnostics(br; max_size = max(1000, max_words))
    notes = ["experimental diagnostic only; no full Zariski closure or algebraic-group classification is computed",
             "finite-field image data are arithmetic evidence and not a proof of characteristic-zero Zariski density",
             "sampled vanishing ideal degree <= $max_degree is not treated as complete"]
    return ZariskiClosureDiagnostics(false, true, dim(br.basis), br.p,
        alg.dimension, comm.dimension, comm.appears_absolutely_irreducible,
        alg.is_full_matrix_algebra, gdiag.determinant_values,
        gdiag.projective_generator_orders, notes)
end
