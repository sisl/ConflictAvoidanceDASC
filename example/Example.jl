push!(LOAD_PATH, "../src")
using Pairwise, Multiagent

# visualize pairwise encounter
d = get_pomdp()
gen_pairwise_policy(d, 1.0, "example")
# viz_pairwise_policy(d)

# load Q(s, a) table and grid
# alpha = read_alpha()
# g = get_grid(d.pomdp.states)

# visualize policy for simple 3-uavs scenario
# viz_policy(alpha, g)

# simple stress test
# uavs = randuavs(3)
# trajs, nlms, nlmsb, time = stress!(utm, maxmin, true, NT, uavs, alpha, g)
# trajs, nlms, nlmsb, time = stress!(uncrd, maxmin, true, NT, uavs, alpha, g)
# trajs, nlms, nlmsb, time = stress!(coord, maxsum, true, NT, uavs, alpha, g)
# trajs, nlms, nlmsb, time = stress!(naive, maxmin, true, NT, uavs, alpha, g)
# trajs, nlms, nlmsb, time = stress!(distr, maxsum, true, NT, uavs, alpha, g)
# plot_trajs(trajs)

# bulk stress test
# nuavs = 2:10
# nsim = 100
# bulk_test(nuavs, nsim, alpha, g)