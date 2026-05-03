"""
Finite-field arithmetic helpers over F_p.

This file contains scalar, root-of-unity, and matrix operations used by
higher arithmetic layers without carrying modular-category state.
"""

using Oscar

# ===== Finite-field arithmetic helpers =====

"""
    FpElem(value, p)

Small exact scalar wrapper for arithmetic in the prime field `F_p`.
It is intentionally lightweight: values are stored as canonical integer
residues, while arithmetic checks that operands live over the same prime.
"""
struct FpElem
    value::Int
    p::Int
    function FpElem(value::Integer, p::Integer)
        pp = Int(p)
        pp > 1 && isprime(pp) || error("FpElem requires a prime p, got $p")
        return new(mod(Int(value), pp), pp)
    end
end

_check_same_field(a::FpElem, b::FpElem) =
    a.p == b.p || error("cannot mix F_$(a.p) and F_$(b.p) elements")

_mul_mod_int(a::Integer, b::Integer, p::Integer) =
    Int(mod(widemul(Int(a), Int(b)), Int(p)))

Base.show(io::IO, x::FpElem) = print(io, x.value, " mod ", x.p)
Base.broadcastable(x::FpElem) = Ref(x)
Base.zero(x::FpElem) = FpElem(0, x.p)
Base.one(x::FpElem) = FpElem(1, x.p)
Base.iszero(x::FpElem) = x.value == 0
Base.:(==)(a::FpElem, b::FpElem) = a.p == b.p && a.value == b.value
Base.hash(x::FpElem, h::UInt) = hash((x.value, x.p), h)
Base.:-(a::FpElem) = FpElem(-a.value, a.p)
function Base.:+(a::FpElem, b::FpElem)
    _check_same_field(a, b)
    return FpElem(a.value + b.value, a.p)
end
function Base.:-(a::FpElem, b::FpElem)
    _check_same_field(a, b)
    return FpElem(a.value - b.value, a.p)
end
function Base.:*(a::FpElem, b::FpElem)
    _check_same_field(a, b)
    return FpElem(_mul_mod_int(a.value, b.value, a.p), a.p)
end
function Base.inv(a::FpElem)
    iszero(a) && error("zero has no inverse in F_$(a.p)")
    return FpElem(invmod(a.value, a.p), a.p)
end
Base.:/(a::FpElem, b::FpElem) = a * inv(b)
Base.:+(a::FpElem, b::Integer) = a + FpElem(b, a.p)
Base.:+(a::Integer, b::FpElem) = FpElem(a, b.p) + b
Base.:-(a::FpElem, b::Integer) = a - FpElem(b, a.p)
Base.:-(a::Integer, b::FpElem) = FpElem(a, b.p) - b
Base.:*(a::FpElem, b::Integer) = a * FpElem(b, a.p)
Base.:*(a::Integer, b::FpElem) = FpElem(a, b.p) * b
Base.:/(a::FpElem, b::Integer) = a / FpElem(b, a.p)
Base.:^(a::FpElem, n::Integer) =
    n >= 0 ? FpElem(powermod(a.value, Int(n), a.p), a.p) : inv(a)^(-n)

function Base.:*(A::AbstractMatrix{FpElem}, B::AbstractMatrix{FpElem})
    size(A, 2) == size(B, 1) ||
        error("matrix dimensions mismatch: $(size(A)) * $(size(B))")
    isempty(A) && return Matrix{FpElem}(undef, size(A, 1), size(B, 2))
    p = A[1, 1].p
    all(x -> x.p == p, A) && all(x -> x.p == p, B) ||
        error("cannot multiply matrices over different prime fields")
    C = fill(FpElem(0, p), size(A, 1), size(B, 2))
    for i in axes(A, 1), j in axes(B, 2)
        s = FpElem(0, p)
        for k in axes(A, 2)
            s += A[i, k] * B[k, j]
        end
        C[i, j] = s
    end
    return C
end

"""
Finite field arithmetic over F_p.

Uses Julia's native Int/BigInt arithmetic with `mod(·, p)` wrapping
rather than a wrapper struct, for minimal overhead in tight loops.

All functions assume p is prime and 2 < p < typemax(Int).
"""

# ----- Euler criterion -----
"""
    is_square(a, p) -> Bool

Test if a is a quadratic residue mod p (p odd prime).
Uses Euler's criterion: a^((p-1)/2) ≡ ±1 (mod p).
"""
function is_square(a::Integer, p::Integer)
    a = mod(a, p)
    a == 0 && return true
    return powermod(a, (p - 1) ÷ 2, p) == 1
end

# ----- Tonelli-Shanks: square root mod p -----
"""
    sqrt_mod(a, p) -> Int

Return some r with r^2 ≡ a (mod p), or throw if a is NQR.
Returns the "smaller" of the two roots (r ≤ p/2) for determinism.
"""
function sqrt_mod(a::Integer, p::Integer)
    a = mod(a, p)
    a == 0 && return 0
    is_square(a, p) || error("$a is not a quadratic residue mod $p")

    # p ≡ 3 (mod 4) easy case
    if p % 4 == 3
        r = powermod(a, (p + 1) ÷ 4, p)
        return min(r, p - r)
    end

    # Tonelli-Shanks for p ≡ 1 (mod 4)
    # Decompose p - 1 = q * 2^s with q odd
    q, s = p - 1, 0
    while q % 2 == 0
        q ÷= 2
        s += 1
    end

    # Find a non-residue z
    z = 2
    while is_square(z, p)
        z += 1
    end

    m = s
    c = powermod(z, q, p)
    t = powermod(a, q, p)
    r = powermod(a, (q + 1) ÷ 2, p)

    while true
        t == 1 && return min(r, p - r)
        # find least i, 0 < i < m, such that t^(2^i) = 1
        i, temp = 0, t
        while temp != 1
            temp = (temp * temp) % p
            i += 1
            i >= m && error("Tonelli-Shanks failed (shouldn't happen)")
        end
        b = powermod(c, 1 << (m - i - 1), p)
        m = i
        c = (b * b) % p
        t = (t * c) % p
        r = (r * b) % p
    end
end

# ----- Primitive root of F_p^* -----
"""
    primitive_root(p) -> Int

Return the smallest positive primitive root modulo p.
"""
function primitive_root(p::Integer)
    p == 2 && return 1
    order_target = p - 1
    # Factor p - 1 to check orders (qualify as Primes.factor to avoid clash with Oscar)
    factors = [f for (f, _) in Primes.factor(order_target)]
    for g in 2:(p-1)
        is_primitive = true
        for q in factors
            if powermod(g, order_target ÷ q, p) == 1
                is_primitive = false
                break
            end
        end
        is_primitive && return Int(g)
    end
    error("no primitive root found mod $p (shouldn't happen for prime p)")
end

# ----- n-th root of unity -----
"""
    root_of_unity(n, p) -> Int

Return a primitive n-th root of unity in F_p.
Requires n | (p - 1).
"""
function root_of_unity(n::Integer, p::Integer)
    n == 1 && return 1
    (p - 1) % n == 0 || error("ζ_$n does not exist in F_$p (need $n | $(p-1))")
    g = primitive_root(p)
    return powermod(g, (p - 1) ÷ n, p)
end

"""
    roots_of_unity(n, p) -> Vector{Int}

Return all n-th roots of unity in F_p, i.e. [ζ_n^0, ζ_n^1, ..., ζ_n^{n-1}].
"""
function roots_of_unity(n::Integer, p::Integer)
    ζ = root_of_unity(n, p)
    result = Vector{Int}(undef, n)
    result[1] = 1
    for k in 2:n
        result[k] = (result[k-1] * ζ) % p
    end
    return result
end

# ----- Oscar polynomial helpers over F_p -----

function _fp_partial_univariate_coeffs(f, varidx::Int,
                                       assigned::Dict{Int, Any}, F)
    coeffs = Dict{Int, Any}()
    for (c, m) in zip(coefficients(f), monomials(f))
        degs = degrees(m)
        term = c
        pow = degs[varidx]
        for i in eachindex(degs)
            i == varidx && continue
            d = degs[i]
            d == 0 && continue
            haskey(assigned, i) || return nothing
            term *= assigned[i]^d
        end
        coeffs[pow] = get(coeffs, pow, zero(F)) + term
    end
    return Dict(k => v for (k, v) in coeffs if !iszero(v))
end

function _fp_roots_from_coeffs(coeffs::AbstractDict{Int}, F, p::Int)
    isempty(coeffs) && return nothing
    maximum(keys(coeffs)) == 0 && return Any[]

    R, x = polynomial_ring(F, :x)
    g = zero(R)
    for (pow, c) in coeffs
        g += R(c) * x^pow
    end
    iszero(g) && return nothing

    roots = Any[]
    for (fac, _) in collect(Oscar.factor(g))
        degree(fac) == 1 || continue
        root = -coeff(fac, 0) / coeff(fac, 1)
        push!(roots, root)
    end
    return sort(unique(roots); by = a -> _fp_elem_to_int(a, p))
end

function _eval_fp_poly(f, vals::Vector)
    F = base_ring(parent(f))
    v = zero(F)
    for (c, m) in zip(coefficients(f), monomials(f))
        term = c
        for (i, d) in enumerate(degrees(m))
            d > 0 && (term *= vals[i]^d)
        end
        v += term
    end
    return v
end

function _fp_elem_to_int(a::FpElem, p::Int)
    a.p == p || error("finite-field element lies in F_$(a.p), expected F_$p")
    return a.value
end

function _fp_elem_to_int(a, p::Int)
    lifted = try
        Oscar.lift(Oscar.ZZ, a)
    catch
        try
            lift(Oscar.ZZ, a)
        catch
            nothing
        end
    end
    lifted === nothing &&
        error("could not lift finite-field element to an integer representative")
    return mod(Int(lifted), p)
end

# ----- Matrix helpers over F_p -----
"""
    matmul_mod(A, B, p) -> Matrix{Int}

Matrix multiplication mod p.
"""
function matmul_mod(A::AbstractMatrix{<:Integer}, B::AbstractMatrix{<:Integer}, p::Integer)
    m, n = size(A)
    n2, k = size(B)
    n == n2 || error("matrix dimensions mismatch: $(size(A)) * $(size(B))")
    pp = Int(p)
    C = zeros(Int, m, k)
    for i in 1:m, j in 1:k
        s = 0
        for ℓ in 1:n
            s = Int(mod(Int128(s) + widemul(Int(A[i, ℓ]), Int(B[ℓ, j])), pp))
        end
        C[i, j] = s
    end
    return C
end

"""
    matpow_mod(A, n, p) -> Matrix{Int}

Compute A^n mod p for square matrix A, using fast exponentiation.
"""
function matpow_mod(A::AbstractMatrix{<:Integer}, n::Integer, p::Integer)
    n == 0 && return Matrix{Int}(I, size(A, 1), size(A, 1))
    n == 1 && return mod.(Int.(A), p)
    result = Matrix{Int}(I, size(A, 1), size(A, 1))
    base = mod.(Int.(A), p)
    while n > 0
        if n & 1 == 1
            result = matmul_mod(result, base, p)
        end
        n >>= 1
        if n > 0
            base = matmul_mod(base, base, p)
        end
    end
    return result
end

"""
    diagmul_right(A, d, p) -> Matrix{Int}

Compute A * diag(d) mod p. Scales column j of A by d[j].
"""
function diagmul_right(A::AbstractMatrix{<:Integer}, d::AbstractVector{<:Integer}, p::Integer)
    m, n = size(A)
    length(d) == n || error("diagonal length mismatch")
    C = zeros(Int, m, n)
    for i in 1:m, j in 1:n
        C[i, j] = (Int(A[i, j]) * Int(d[j])) % p
    end
    return C
end

"""
    diagmul_left(d, A, p) -> Matrix{Int}

Compute diag(d) * A mod p. Scales row i of A by d[i].
"""
function diagmul_left(d::AbstractVector{<:Integer}, A::AbstractMatrix{<:Integer}, p::Integer)
    m, n = size(A)
    length(d) == m || error("diagonal length mismatch")
    C = zeros(Int, m, n)
    for i in 1:m, j in 1:n
        C[i, j] = (Int(d[i]) * Int(A[i, j])) % p
    end
    return C
end

"""
    lift_symmetric(x, p) -> Int

Lift x ∈ F_p = {0,...,p-1} to the symmetric range [-p/2, p/2].
Useful for recognizing small integers (e.g., fusion coefficients).
"""
function lift_symmetric(x::Integer, p::Integer)
    x = mod(x, p)
    return x <= p ÷ 2 ? Int(x) : Int(x) - Int(p)
end
