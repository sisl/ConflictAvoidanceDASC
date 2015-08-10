# load dependencies
push!(LOAD_PATH, "../dvi")

addprocs(int(CPU_CORES / 2))

using DiscreteValueIteration, JLD, PilotSCAs

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
    maxIterations = 1000,
    tolerance = 1e-4,
    gaussSiedel = true,
    includeV = true,
    includeQ = true,
    includeA = true)

println("\nStarting parallel solver...")
policy = solve(solver, mdp, verbose = true)
println("\nParallel solution generated!")

# save solution
function sharray2array(sharray::SharedArray{Float64, 2})
    result = zeros(sharray.dims)
    for i = 1:sharray.dims[1]
        for j = 1:sharray.dims[2]
            result[i, j] = sharray[i, j]
        end # for j
    end # for i
    return result
end # function sharray2array

solQ = sharray2array(policy.Q')
save("../../data/pilot-alpha.jld", "solQ", solQ)
println("Parallel solution saved...exiting PilotSCA.jl script.")