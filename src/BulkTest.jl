using Pairwise, Multiagent

alphafile = string("../data/alpha-tula", ARGS[1], ".jld")
tag = int(ARGS[2])

d = get_pomdp()
g = get_grid(d.pomdp.states)

nuavs = 2:10
nsim = 1000

bulk_test(nuavs, alphafile, g, string("tula", tag), nsim)