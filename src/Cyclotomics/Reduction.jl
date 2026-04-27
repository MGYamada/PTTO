"""
Finite-field reduction helpers for cyclotomic data.

This layer chooses primitive roots in F_p and reduces Q(zeta_N) elements,
vectors, and matrices to the finite fields used by search and reconstruction.
"""

# ============================================================
#  Step 2: F_p reduction of Q(ζ_N) elements
# ============================================================

"""
    find_zeta_in_Fp(N::Int, p::Int) -> Int

Find a primitive N-th root of unity in F_p. Requires N | p-1.
Returns the F_p representative as an Int.

Uses: pick a primitive root g of F_p, then ζ_N = g^{(p-1)/N}.
"""
function find_zeta_in_Fp(N::Int, p::Int)
    (p - 1) % N == 0 || error("N = $N does not divide p - 1 = $(p-1)")
    g = primitive_root(p)
    exponent = div(p - 1, N)
    return powermod(g, exponent, p)
end

"""
    cyclotomic_to_Fp(x, zeta_N_Fp::Int, p::Int) -> Int

Reduce an Oscar cyclotomic field element x ∈ Q(ζ_N) to F_p, given
the value of ζ_N in F_p.

Works by:
1. Extracting coefficients of x in the basis {1, ζ, ζ², ..., ζ^{φ(N)-1}}
2. Reducing each rational coefficient mod p (requires denominator coprime to p)
3. Summing as linear combination of powers of zeta_N_Fp
"""
function cyclotomic_to_Fp(x, zeta_N_Fp::Int, p::Int)
    # Get the coefficient vector of x in the power basis of Q(ζ_N)
    coeffs = Oscar.coefficients(x)
    result = 0
    for (k, c) in enumerate(coeffs)
        # c is a rational number, k-1 is the power of ζ
        num = Int(numerator(c))
        den = Int(denominator(c))
        den % p != 0 || error("Denominator $den divisible by p = $p")
        coeff_Fp = mod(num * invmod(den, p), p)
        zeta_pow = powermod(zeta_N_Fp, k - 1, p)
        result = mod(result + coeff_Fp * zeta_pow, p)
    end
    return result
end

"""
    reduce_matrix_to_Fp(M, zeta_N_Fp::Int, p::Int) -> Matrix{Int}

Reduce an Oscar matrix over Q(ζ_N) to an Int matrix representing values mod p.
"""
function reduce_matrix_to_Fp(M, zeta_N_Fp::Int, p::Int)
    n, m = size(M)
    result = zeros(Int, n, m)
    for i in 1:n
        for j in 1:m
            result[i, j] = cyclotomic_to_Fp(M[i, j], zeta_N_Fp, p)
        end
    end
    return result
end

"""
    reduce_vector_to_Fp(v::Vector, zeta_N_Fp::Int, p::Int) -> Vector{Int}
"""
function reduce_vector_to_Fp(v::Vector, zeta_N_Fp::Int, p::Int)
    return [cyclotomic_to_Fp(x, zeta_N_Fp, p) for x in v]
end
