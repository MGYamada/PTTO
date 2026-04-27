"""
Validation and construction helpers for modular data over finite fields.

This layer checks S/T axioms and wraps validated data as ModularDatumFp.
"""

"""
Modular data validation and construction over F_p.

A modular datum (S, T) must satisfy:
- S is symmetric: S[i,j] = S[j,i]
- S² = C (charge conjugation permutation matrix), C² = I
- T is diagonal with T[0,0] = 1, entries are N-th roots of unity
- (ST)³ = α · S² for some α ∈ μ_{8N}(F_p) (central charge phase)
- S[i, 0] ≠ 0 for all i (so quantum dimensions are well-defined)

This section provides:
- `validate_modular_data`: check all axioms given (S, T, p, N)
- `build_modular_datum`: construct ModularDatumFp from (S, T, p, N) with validation
- `compute_alpha`: compute α from (ST)³ and S²
- `compute_charge_conjugation`: extract C from S²
"""

@inline _invalid_modular_data(reason::String; α::Int = 0, C::Vector{Int} = Int[]) =
    (valid=false, reason=reason, α=α, C=C, D_squared=0, D=0)

"""
    compute_charge_conjugation(S::Matrix{Int}, p::Int) -> Union{Vector{Int}, Nothing}

Given S mod p, compute S² and extract charge conjugation permutation C,
where S² = C (as permutation matrix). Returns C as a vector where C[i] is
the image of i under the permutation (1-indexed).

Returns nothing if S² is not a permutation matrix.
"""
function compute_charge_conjugation(S::AbstractMatrix{Int}, p::Int)
    r = size(S, 1)
    S2 = matmul_mod(S, S, p)
    C = zeros(Int, r)
    for i in 1:r
        found = 0
        for j in 1:r
            if S2[i, j] == 1
                found == 0 || return nothing  # multiple 1s in row
                found = j
            elseif S2[i, j] != 0
                return nothing  # non-0, non-1 entry
            end
        end
        found == 0 && return nothing  # no 1 in row
        C[i] = found
    end
    # Check it's an involution
    for i in 1:r
        C[C[i]] == i || return nothing
    end
    return C
end

"""
    compute_alpha(S::Matrix{Int}, T::Vector{Int}, p::Int) -> Union{Int, Nothing}

Compute α such that (ST)³ = α · S². Returns α ∈ F_p if (ST)³ is indeed a
scalar multiple of S² (entry-wise), otherwise nothing.
"""
function compute_alpha(S::AbstractMatrix{Int}, T::AbstractVector{Int}, p::Int)
    r = size(S, 1)
    ST = diagmul_right(S, T, p)
    ST3 = matpow_mod(ST, 3, p)
    S2 = matmul_mod(S, S, p)

    # Find α: (ST)³[i,j] = α * S²[i,j] for all i,j
    α = -1  # sentinel
    for i in 1:r, j in 1:r
        lhs = ST3[i, j]
        rhs = S2[i, j]
        if rhs == 0
            lhs == 0 || return nothing
        else
            α_candidate = (lhs * invmod(rhs, p)) % p
            if α == -1
                α = α_candidate
            elseif α != α_candidate
                return nothing  # inconsistent
            end
        end
    end
    return α == -1 ? 1 : α  # if S² ≡ 0 (shouldn't happen for modular data), default 1
end

"""
    validate_modular_data(S, T, p, N) -> NamedTuple

Perform all modular datum axiom checks. Returns a NamedTuple:
- valid::Bool
- reason::String (empty if valid, else describes first failure)
- α::Int (central charge phase, if computable)
- C::Vector{Int} (charge conjugation, if computable)
- D_squared::Int
- D::Int (0 if D² is not a square)
"""
function validate_modular_data(S::AbstractMatrix{Int}, T::AbstractVector{Int}, p::Int, N::Int)
    r = size(S, 1)
    # Shape checks
    size(S) == (r, r) || return _invalid_modular_data("S not square")
    length(T) == r || return _invalid_modular_data("T wrong length")
    # Symmetry
    for i in 1:r, j in (i+1):r
        S[i, j] == S[j, i] || return _invalid_modular_data("S not symmetric at ($i,$j)")
    end
    # T[1] = 1 (unit)
    T[1] == 1 || return _invalid_modular_data("T[1] ≠ 1")
    # T entries are N-th roots of unity
    for i in 1:r
        powermod(T[i], N, p) == 1 || return _invalid_modular_data("T[$i] is not an N-th root of unity")
    end
    # S[1, i] ≠ 0 (quantum dimensions well-defined)
    for i in 1:r
        S[1, i] != 0 || return _invalid_modular_data("S[1,$i] = 0")
    end
    # S² = C (permutation matrix, involution)
    C = compute_charge_conjugation(S, p)
    C === nothing && return _invalid_modular_data("S² is not a permutation matrix")
    # (ST)³ = α · S²
    α = compute_alpha(S, T, p)
    α === nothing && return _invalid_modular_data("(ST)³ not proportional to S²"; C=C)
    # Compute D² from quantum dimensions
    # d_i = S[i, 1] / S[1, 1], D² = Σ d_i²
    S00_inv = invmod(S[1, 1], p)
    D_squared = 0
    for i in 1:r
        d_i = (Int(S[i, 1]) * S00_inv) % p
        D_squared = (D_squared + d_i * d_i) % p
    end
    # D is sqrt of D² if exists
    D = is_square(D_squared, p) ? sqrt_mod(D_squared, p) : 0
    return (valid=true, reason="", α=α, C=C, D_squared=D_squared, D=D)
end

"""
    build_modular_datum(S, T, p, N) -> Union{ModularDatumFp, Nothing}

Construct a ModularDatumFp from raw (S, T, p, N) with full validation.
Returns nothing if any axiom fails.
"""
function build_modular_datum(S::AbstractMatrix{Int}, T::AbstractVector{Int}, p::Int, N::Int)
    r = size(S, 1)
    S_mod = mod.(Int.(S), p)
    T_mod = mod.(Int.(T), p)
    result = validate_modular_data(S_mod, T_mod, p, N)
    result.valid || return nothing
    return ModularDatumFp(r, p, N, S_mod, T_mod, result.α, result.D_squared, result.D, result.C)
end
