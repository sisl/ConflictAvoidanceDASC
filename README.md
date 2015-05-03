# Short-term Conflict Avoidance for Unmanned Aircraft Traffic Management

This directory contains the supplement code for the paper titled "Short-term Conflict Avoidance for Unmanned Aircraft Traffic Management" by Hao Yi Ong and Mykel J. Kochenderfer, submitted to the 2015 Digital Avionics Systems Conference. 

Here you will find implementations of the following:
* stuff
* more stuff
* even more stuff

### Dependencies
The software is implemented entirely in Julia. For the best results, the user should use IJulia. The following packages are required for running all code. 
* PGFPlots
* another package
* yet another package

### Generate policy
To generate the MDP policy for the pairwise encounter subproblem, run the following command. Warning: This can take a while. 
```
julia get_policy
```

To use a pre-generated version of the pairwise encounter, run the following command in the Julia interpreter.
```
julia
policy = read_policy()
```

### Visualize policy as heatmap
To visualize the policy for the pairwise encounter generated in the paper, run the following command.
```
julia viz_heatmap
```

To visualize the policy for the three-aircraft encounter generated in the paper, run the following command.

### Run simulations
To run a single instance encounter and visualize the trajectories, run the following command.
```
julia sim_encounter arg1 arg2
```

To run a batch of simulations and plot the results, run the following command.
```
julia sim_encounter arg1
```

### Layout
```
src/
    stuff1.jl
    stuff2.hs

data/
    policy.jld
    
```