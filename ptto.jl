using Oscar

function main(N)
    R, = residue_ring(ZZ, N)
    G = special_linear_group(2, R)
    println(order(G))
end

main(5)
