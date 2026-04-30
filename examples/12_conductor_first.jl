# Documentation example: conductor-first exact modular data.
#
# This script starts from a conductor, constructs exact cyclotomic modular
# data, and checks a few lightweight invariants.  It is intended to be safe for
# CI and documentation smoke tests.

using ACMG

function main()
    ctx = CyclotomicContext(20)
    data = fibonacci_modular_data(ctx)

    @assert conductor(ctx) == 20
    @assert conductor(data) == 20
    @assert cond_S(data) == 20
    @assert cond_T(data) == 5
    @assert cond_F(data) === nothing

    exact = validate_exact_modular_data(data)
    @assert exact.ok

    dims = quantum_dimensions(data)
    @assert length(dims) == 2

    hcc = higher_central_charge(data; n = 1)
    @assert hcc.ok
    @assert hcc.conductor == 20

    orbit = galois_orbit(data)
    @assert !isempty(orbit)

    fp = reduce_mod_p(data, 41)
    @assert fp.p == 41
    @assert size(fp.S) == (2, 2)

    println("Fibonacci over Q(zeta_20): rank=", length(data.labels))
    println("higher central charge status: ", hcc.status)
    println("finite-field prime: ", fp.p)
end

main()
