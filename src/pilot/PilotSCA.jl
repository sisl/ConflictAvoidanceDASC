# load dependencies
push!(LOAD_PATH, "../dvi")

addprocs(int(CPU_CORES / 2))

using DiscreteValueIteration, JLD, PilotSCAs, PilotSCAViz

mdp = SCA()

# check size of mdp
function getBytes(x)
   total = 0;
   fieldNames = typeof(x).names;
   if fieldNames == ()
      return sizeof(x);
   else
     for fieldName in fieldNames
        total += getBytes(getfield(x,fieldName));
     end
     return total;
   end
end

println("mdp of type ", typeof(mdp), " takes up ", getBytes(mdp) / 1000.0, " kB")

# informal validation of transition function
nextStateIndices, probs = nextStates(mdp, 1, 15)
println("next state indices:\n", nextStateIndices, "\n")
println("probabilities:\n", probs, "\n")
println("probabilities sum to ", sum(probs))

# parallel solution
numProcs = int(CPU_CORES / 2)
solver = ParallelSolver(
    numProcs,
    maxIterations = 100,
    tolerance = 1e-2,
    gaussSiedel = true,
    includeV = true,
    includeQ = true,
    includeA = true)

policy = solve(solver, mdp, verbose = true)

# save solution
solQ = sharray2array(policy.Q')
save("../../data/pilot-alpha.jld", "solQ", solQ)

println("End of PilotSCA.jl script.")