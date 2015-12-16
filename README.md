# Short-Term Conflict Avoidance for Unmanned Aircraft Traffic Management

This directory contains the supplement code for the paper titled "Short-Term Conflict Avoidance for Unmanned Aircraft Traffic Management" by Hao Yi Ong and Mykel J. Kochenderfer, published at the 2015 Digital Avionics Systems Conference. 

Here you will find implementations of the following:
* Pairwise encounter policy generation and visualization
* Multithreat encounter example encounter policy visualization
* Multithreat encounter trajectory simulation and bulk tests

### Pilot response model

In order to better capture the response of unmanned aircraft, be it remotely or autonomously piloted, we incorporated a pilot response model for the pairwise encounter policy generation. Details can be found in the `src/pilot` directory.

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
    dvi/
        DiscreteValueIteration.jl
        helpers.jl
        parallel.jl
        policy.jl
        serial.jl

    serial/
        DoubleUAVs.jl
        Multiagent.jl
        Pairwise.jl
        POMDPs.jl
        QMDP.jl
        UAVs.jl

    parallel/
        ParSCA.jl
        ParSCA.ipynb
        SCAs.jl
        SCAConst.jl
        SCAIterators.jl
        README.md

    pilot/
        PilotSCA.jl
        PilotSCA.ipynb
        PilotSCAConst.jl
        PilotSCAViz.jl
        PilotSCAs.jl
        README.md

data/
    alpha.jld
    results.jld

example/
    Example.jl
    Example.ipynb

README.md
```
