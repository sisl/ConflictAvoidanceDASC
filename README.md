# Short-Term Conflict Avoidance for Unmanned Aircraft Traffic Management

This directory contains the supplement code for the paper titled "Short-Term Conflict Avoidance for Unmanned Aircraft Traffic Management" by Hao Yi Ong and Mykel J. Kochenderfer, to appear in the 2015 Digital Avionics Systems Conference. 

Here you will find implementations of the following:
* Pairwise encounter policy generation and visualization (serial and parallel)
* Multithreat encounter example encounter policy visualization
* Multithreat encounter trajectory simulation and bulk tests

### Dependencies
The software is implemented entirely in Julia. For the best results, the user should use a notebook. An example notebook is provided for the reader's convenience in the example subdirectory. The following Julia packages are required for running all code. 
* [PGFPlots](https://github.com/sisl/PGFPlots.jl)
* [GridInterpolations](https://github.com/sisl/GridInterpolations.jl)
* [ODE](https://github.com/JuliaLang/ODE.jl)
* [Interact](https://github.com/JuliaLang/Interact.jl)
* [JLD](https://github.com/JuliaLang/JLD.jl)
* [DiscreteMDPs](https://github.com/sisl/DiscreteMDPs.jl)

### Layout
```
src/
    DoubleUAVs.jl
    Multiagent.jl
    Pairwise.jl
    POMDPs.jl
    QMDP.jl
    UAVs.jl

data/
    alpha.jld
    results.jld

parallel/
    ParSCA.jl
    ParSCA.ipynb
    SCAs.jl
    SCAConst.jl
    SCAIterators.jl
    README.md

    DiscreteValueIteration/
        DiscreteValueIteration.jl
        helpers.jl
        parallel.jl
        policy.jl
        serial.jl

example/
    Example.jl
    Example.ipynb

README.md
```