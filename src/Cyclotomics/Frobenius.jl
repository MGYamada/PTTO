"""
Frobenius action on exact modular data.

Frobenius is represented as the corresponding Galois action modulo the
cyclotomic conductor.
"""

frobenius(data::ModularData, p::Int) = galois_action(data, mod(p, data.context.N))
