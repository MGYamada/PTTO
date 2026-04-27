"""
    select_admissible_primes(N_eff; min_count = 4, window = 2000, start_from = 29)
        -> Vector{Int}

Collect admissible primes `p` such that `(p - 1) % N_eff == 0`, scanning with
`Primes.nextprime` over `(start_from, start_from + window]`.

Notes:
- Enforces `min_count ≥ 2` (CRT needs at least two primes).
- Typical recommendation is `min_count ≥ 4` for robust used/fresh splitting.
- Throws a detailed error if not enough admissible primes are found in the
  attempted range.
"""
function select_admissible_primes(N_eff::Int;
                                  min_count::Int = 4,
                                  window::Int = 2000,
                                  start_from::Int = 29)
    N_eff >= 1 || error("N_eff must be positive, got $N_eff")
    min_count >= 2 || error("min_count must be ≥ 2, got $min_count")
    window >= 1 || error("window must be ≥ 1, got $window")

    upper = start_from + window
    chosen = Int[]

    p = nextprime(start_from)
    while p <= upper && length(chosen) < min_count
        if (p - 1) % N_eff == 0
            push!(chosen, p)
        end
        p = nextprime(p + 1)
    end

    if length(chosen) < min_count
        found = length(chosen)
        error("insufficient admissible primes for N_eff=$N_eff: " *
              "searched range ($(start_from), $upper], found $found < $min_count. " *
              "Reason: condition (p-1) % N_eff == 0 is too restrictive in this window; " *
              "increase `window` or lower `start_from`.")
    end

    return chosen
end
