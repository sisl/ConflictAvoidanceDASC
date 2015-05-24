using Pairwise, Multiagent

tag = int(ARGS[1])
alphafile = string("../data/alpha-tula", tag, ".jld")

d = get_pomdp()
g = get_grid(d.pomdp.states)

nuavs = 2:10
nsim = 1000
nbulk = 10
bulk_test(nuavs, alphafile, g, string("tula", tag), nsim, nbulk)