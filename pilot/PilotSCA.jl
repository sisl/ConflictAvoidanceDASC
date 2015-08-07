# load modules
push!(LOAD_PATH, "./DiscreteValueIteration")

addprocs(int(CPU_CORES / 2))

using DiscreteValueIteration
using PilotSCAs, PilotSCAIterators, PilotSCAConst
using GridInterpolations

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

# parallel solution
numProcs = int(CPU_CORES / 2)
solver = ParallelSolver(
    numProcs,
    maxIterations = 1000,
    tolerance = 1e-6,
    gaussSiedel = false,
    includeV = true,
    includeQ = true,
    includeA = true)

policy = solve(solver, mdp, verbose = true)

# save solution
using JLD

solQ = policy.Q'
save("../data/par-alpha.jld", "solQ", solQ)