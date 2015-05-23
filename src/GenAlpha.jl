using Pairwise, Multiagent

d = get_pomdp()
lambdas = logspace(-3, 2, 21)

println("Starting QMDPs for lambda: ", lambdas)
for tag = 1:length(lambdas)
    println("lamba: ", lamdbas[tag])
    gen_pairwise_policy(d, lambdas[tag], tag)
end # tag

println("We're done generating all policy files!")