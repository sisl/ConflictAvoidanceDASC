using Pairwise, Multiagent

d = get_pomdp()
lambda = float(ARGS[1])
tag = int(ARGS[2])

println("Starting QMDPs for lambda = ", lambda)
gen_pairwise_policy(d, lambda, string("tula", tag))