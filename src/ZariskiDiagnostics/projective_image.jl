projective_orders(br::FiniteFieldBraidRepresentation; max_order = 100000) =
    [_projective_order(g, br.p; max_order = max_order) for g in br.generators]
