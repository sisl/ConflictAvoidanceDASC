using Pairwise, Multiagent

d = get_pomdp()
lambdas = logspace(-3, 2, 21)

for tag = 1:length(lambdas)
    gen_pairwise_policy(d, lambdas[tag], tag)
end # tag

println("We're done generating all policy files!")