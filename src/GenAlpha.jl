using Pairwise, Multiagent

d = get_pomdp()
lambda = float(ARGS[1])

println("Starting QMDPs for lambda = ", lambda)
gen_pairwise_policy(d, string("tula", lambdas[tag]), tag)