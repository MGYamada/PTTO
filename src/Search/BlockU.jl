"""
Block-U parametrisation and finite-field MTC reconstruction.

This module implements the core fixed-stratum finite-field search:

1. `build_block_diagonal`: given a stratum {m_λ} and atomic catalog,
   build the block-diagonal atomic (S, T) on V = ⊕_λ V_λ^{m_λ}.

2. `reduce_to_Fp`: reduce an Oscar Q(ζ_N)-matrix to F_p (requires N | p-1
   and a chosen primitive N-th root of unity in F_p).

3. `t_eigenspace_decomposition`: given T (diagonal over F_p), group
   indices by eigenvalue. Returns a Dict{Int, Vector{Int}} mapping
   each T-eigenvalue to the list of basis indices.

4. `sweep_O2`: for a 2-dimensional degenerate eigenspace with indices
   (i, j), enumerate the O(2)(F_p)-rational circle points (u, v) with
   u² + v² = 1, and for each compute S' = U^T S U.

5. `verlinde_check`: compute fusion coefficients N_{ij}^k from an F_p
   modular datum (S', T'), check all are small non-negative integers.

The top-level driver is `find_mtcs_at_prime(catalog, stratum, p)`,
which runs the whole pipeline at a single prime.

Follow-up functions (in a separate CRT module, TODO) will combine
results across multiple primes.
"""

using Oscar

include("BlockU/Types.jl")
include("BlockU/FiniteField.jl")
include("BlockU/Canonicalization.jl")
include("BlockU/Blocks.jl")
include("BlockU/Constraints.jl")
include("BlockU/Enumeration.jl")
include("BlockU/Search.jl")
