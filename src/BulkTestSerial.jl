using Pairwise, Multiagent

for tag = 3:11
    println("bulk test serial no. ", tag)
    tic()
    alphafile = string("../data/alpha-tula", tag, ".jld")

    d = get_pomdp()
    g = get_grid(d.pomdp.states)

    nuavs = 2:10
    nsim = 1000
    nbulk = 1
    bulk_test(nuavs, alphafile, g, string("cheonan", tag), nsim, nbulk)
    println("bulk test serial no. ", tag, "done in ", toq())
end # for tag