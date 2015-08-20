module PilotSCAViz

export viz_pairwise_policy, DoubleUAV, sharray2array, print_pairwise_policy

import PilotSCAConst: Xmin, Xmax, Ymin, Ymax, Bearingmin, Bearingmax, Speedmin, Speedmax
import PilotSCAConst: Xdim, Ydim, Bearingdim, Speeddim, Responsedim, NStates

using GridInterpolations, Interact, PGFPlots


const STATE_DIM = 7
const ACTION_DIM = 2

const X = 1                 # [m] relative x-position
const Y = 2                 # [m] relative y-position
const P = 3                 # [rad] relative heading
const V1 = 4                # [m/s] ac1 speed
const V2 = 5                # [m/s] ac2 speed
const R1 = 6                # [] ac1 response
const R2 = 7                # [] ac2 response

const XMIN = Xmin
const XMAX = Xmax

const YMIN = Ymin
const YMAX = Ymax

const PMIN = Bearingmin
const PMAX = Bearingmax

const VMIN = Speedmin
const VMAX = Speedmax

const XDIM = Xdim
const YDIM = Ydim
const PDIM = Bearingdim
const VDIM = Speeddim
const RDIM = Responsedim

xs = linspace(XMIN, XMAX, XDIM)
ys = linspace(YMIN, YMAX, YDIM)
ps = linspace(PMIN, PMAX, PDIM)
vs = linspace(VMIN, VMAX, VDIM)
rs = 0:1

const NSTATES = NStates
tmp_states = zeros(STATE_DIM, NSTATES)
istate = 1
for ir2 = 1:RDIM
    for ir1 = 1:RDIM
        for iv2 = 1:VDIM
            for iv1 = 1:VDIM
                for ip = 1:PDIM
                    for iy = 1:YDIM
                        for ix = 1:XDIM
                            tmp_states[X, istate] = xs[ix]
                            tmp_states[Y, istate] = ys[iy]
                            tmp_states[P, istate] = ps[ip]
                            tmp_states[V1, istate] = vs[iv1]
                            tmp_states[V2, istate] = vs[iv2]
                            tmp_states[R1, istate] = rs[ir1]
                            tmp_states[R2, istate] = rs[ir2]
                            istate += 1
                        end # for ix
                    end # for iy
                end # for ip
            end # for iv1
        end # for iv2
    end # for ir1
end # for ir2

const TERM_STATE_VAR = 1e5
const TERM_STATE = TERM_STATE_VAR * ones(STATE_DIM)
tmp_states[:, end] = TERM_STATE
const STATES_OBSERVATIONS = tmp_states

const COC = -1
indivActions = deg2rad([-20, -10, 0, 10, 20, rad2deg(COC)])
nIndivActions = length(indivActions)
actions = zeros(2, nIndivActions^2)
iaction = 1
for ia2 = 1:nIndivActions
    for ia1 = 1:nIndivActions
        actions[1, iaction] = indivActions[ia1]
        actions[2, iaction] = indivActions[ia2]
        iaction += 1
    end # for ia1
end # for ia2
const ACTIONS = actions

const SIGMA_V = 2.0                 # [m/s]
const SIGMA_B = deg2rad(4.0)        # [rad]
const SIGMA_B_COC = deg2rad(10.0)   # [rad]
const SIGMAS = [0 SIGMA_B -SIGMA_B 0 0 0 0 0 0;
                0 0 0 SIGMA_B -SIGMA_B 0 0 0 0;
                0 0 0 0 0 SIGMA_V -SIGMA_V 0 0;
                0 0 0 0 0 0 0 SIGMA_V -SIGMA_V]
const WEIGHTS = [1 / 3, 2 / (3 * (size(SIGMAS, 2) - 1)) * 
                             ones(size(SIGMAS, 2) - 1)]
const SIGMA_DEF = zeros(size(SIGMAS, 1))

const NBIN = 500


type POMDP

    states          :: Matrix{Float64}
    actions         :: Array{Float64}
    observations    :: Matrix{Float64}
    term_state      :: Vector{Float64}
    
    function POMDP(states::Matrix, actions::Array, observations::Matrix, term_state::Vector{Float64})
        return new(states, actions, observations, term_state)
    end # function POMDP
    
end # type POMDP


type DoubleUAV

    pomdp::     POMDP
    sigmas::    Matrix{Float64}
    weights::   Vector{Float64}
    
    function DoubleUAV()
        p = POMDP(STATES_OBSERVATIONS, ACTIONS, STATES_OBSERVATIONS, TERM_STATE)
        return new(p, SIGMAS, WEIGHTS)
    end # function DoubleUAV

end # DoubleUAV


type Policy
    alpha       :: Matrix{Float64}
    actions     :: Matrix{Float64}
    nactions    :: Int64
    qvals       :: Vector{Float64}

    function Policy(alpha::Matrix{Float64}, actions::Matrix{Float64})
        return new(alpha, actions, size(actions, 2), zeros(size(actions, 2)))
    end # function Policy
end


function read_policy(actions::Matrix{Float64}, alpha::Matrix{Float64})
    return Policy(alpha, actions)
end # function read_policy


function evaluate(policy::Policy, belief::SparseMatrixCSC{Float64,Int64})
    fill!(policy.qvals, 0.0)
    get_qval!(policy, belief)
    ibest = indmax(policy.qvals)
    return policy.actions[:, ibest], ibest
end # function evaluate


function get_qval!(policy::Policy, belief::SparseMatrixCSC{Float64, Int64})
    fill!(policy.qvals, 0.0)
    for iaction in 1:policy.nactions
        for ib in 1:length(belief.rowval)
            policy.qvals[iaction] += belief.nzval[ib] * policy.alpha[belief.rowval[ib], iaction]
        end # for b
    end # for iaction
end # function get_qval!


function get_belief(pstate::Vector{Float64}, grid::RectangleGrid)
    belief = spzeros(NSTATES, 1)
    if pstate[1] == TERM_STATE_VAR
        belief[end] = 1.0
    else
        indices, weights = interpolants(grid, pstate)
        for i = 1:length(indices)
            belief[indices[i]] = weights[i]
        end # for i
    end # if
    return belief
end # function get_belief


function viz_pairwise_policy(
        d::DoubleUAV,
        alpha::Matrix{Float64},
        nbin::Int64=250)

    states = d.pomdp.states
    actions = d.pomdp.actions

    xs = sort(unique(vec(states[X, 1:end - 1])))
    ys = sort(unique(vec(states[Y, 1:end - 1])))
    ps = sort(unique(vec(states[P, 1:end - 1])))
    vs = sort(unique(vec(states[V1, 1:end - 1])))
    rs = sort(unique(vec(states[R1, 1:end - 1])))

    grid = RectangleGrid(xs, ys, ps, vs, vs, rs, rs)

    ps_deg = rad2deg(ps)
    pstart = ps_deg[1]
    pend = ps_deg[end]
    pdiv = ps_deg[2] - ps_deg[1]
    
    vstart = vs[1]
    vend = vs[end]
    vdiv = vs[2] - vs[1]

    policy = read_policy(actions, alpha)

    @manipulate for p = pstart:pdiv:pend, 
                    v0 = vstart:vdiv:vend, 
                    v1 = vstart:vdiv:vend,
                    r0 = 0:1,
                    r1 = 0:1
        
        # ownship uav
        function get_heat1(x::Float64, y::Float64)
            
            action, _ = evaluate(policy, get_belief(
                [x, y, deg2rad(p), v0, v1, r0, r1], grid))
            
            if abs(action[1]) > 0.5
                return 2.0
            end # if
            
            return rad2deg(action[1])

        end # function get_heat1
        
        # intruder uav
        function get_heat2(x::Float64, y::Float64)
            
            action, _ = evaluate(policy, get_belief(
                [x, y, deg2rad(p), v0, v1, r0, r1], grid))
            
            if abs(action[2]) > 0.5
                return 2.0
            end # if
            
            return rad2deg(action[2])

        end # function get_heat2
        
        g = GroupPlot(2, 1, groupStyle = "horizontal sep=3cm")
        push!(g, Axis([
            Plots.Image(get_heat1, (int(XMIN), int(XMAX)), 
                        (int(YMIN), int(YMAX)), 
                        zmin = -20, zmax = 20,
                        xbins = nbin, ybins = nbin,
                        colormap = ColorMaps.Named("RdBu"), colorbar = false),
            Plots.Node(L">", 0, 0, style="rotate=0,font=\\Huge"),
            Plots.Node(L">", 2600, 2600, style=string("rotate=", p, ",font=\\Huge"))
            ], width="10cm", height="10cm", xlabel="x (m)", ylabel="y (m)", title="Ownship action"))
        push!(g, Axis([
            Plots.Image(get_heat2, (int(XMIN), int(XMAX)), 
                        (int(YMIN), int(YMAX)), 
                        zmin = -20, zmax = 20,
                        xbins = nbin, ybins = nbin,
                        colormap = ColorMaps.Named("RdBu")),
            Plots.Node(L">", 0, 0, style="rotate=0,font=\\Huge"),
            Plots.Node(L">", 2600, 2600, style=string("rotate=", p, ",font=\\Huge"))],
            width="10cm", height="10cm", xlabel="x (m)", title="Intruder action"))
        g
    end # for p, v0, v1, r0, r1

end # function viz_pairwise_policy


function print_pairwise_policy(
        d::DoubleUAV,
        alpha::Matrix{Float64},
        pval::Float64 = 180.0,
        v0val::Float64 = 10.0,
        v1val::Float64 = 10.0)

    nbin = NBIN

    states = d.pomdp.states
    actions = d.pomdp.actions

    xs = sort(unique(vec(states[X, 1:end - 1])))
    ys = sort(unique(vec(states[Y, 1:end - 1])))
    ps = sort(unique(vec(states[P, 1:end - 1])))
    vs = sort(unique(vec(states[V1, 1:end - 1])))
    rs = sort(unique(vec(states[R1, 1:end - 1])))

    grid = RectangleGrid(xs, ys, ps, vs, vs, rs, rs)

    ps_deg = rad2deg(ps)
    pstart = ps_deg[1]
    pend = ps_deg[end]
    pdiv = ps_deg[2] - ps_deg[1]
    
    vstart = vs[1]
    vend = vs[end]
    vdiv = vs[2] - vs[1]

    policy = read_policy(actions, alpha)

    for r0 = 0:1, r1 = 0:1
        
        p = pval
        v0 = v0val
        v1 = v1val
        
        # ownship uav
        function get_heat1(x::Float64, y::Float64)
            
            action, _ = evaluate(policy, get_belief(
                [x, y, deg2rad(p), v0, v1, r0, r1], grid))
            
            if abs(action[1]) > 0.5
                return 2.0
            end # if
            
            return rad2deg(action[1])

        end # function get_heat1
        
        # intruder uav
        function get_heat2(x::Float64, y::Float64)
            
            action, _ = evaluate(policy, get_belief(
                [x, y, deg2rad(p), v0, v1, r0, r1], grid))
            
            if abs(action[2]) > 0.5
                return 2.0
            end # if
            
            return rad2deg(action[2])

        end # function get_heat2
        
        g = GroupPlot(2, 1, groupStyle = "horizontal sep=3cm")
        push!(g, Axis([
            Plots.Image(get_heat1, (int(XMIN), int(XMAX)), 
                        (int(YMIN), int(YMAX)), 
                        zmin = -20, zmax = 20,
                        xbins = nbin, ybins = nbin,
                        colormap = ColorMaps.Named("RdBu"), colorbar = false),
            Plots.Node(L">", 0, 0, style="rotate=0,font=\\Huge"),
            Plots.Node(L">", 1800, 1800, style=string("rotate=", p, ",font=\\Huge"))
            ], width="10cm", height="10cm", xlabel="x (m)", ylabel="y (m)", title="Ownship action"))
        push!(g, Axis([
            Plots.Image(get_heat2, (int(XMIN), int(XMAX)), 
                        (int(YMIN), int(YMAX)), 
                        zmin = -20, zmax = 20,
                        xbins = nbin, ybins = nbin,
                        colormap = ColorMaps.Named("RdBu")),
            Plots.Node(L">", 0, 0, style="rotate=0,font=\\Huge"),
            Plots.Node(L">", 1800, 1800, style=string("rotate=", p, ",font=\\Huge"))],
            width="10cm", height="10cm", xlabel="x (m)", title="Intruder action"))
        PGFPlots.save(string("pair-", r0, r1, ".tex"), g)
    end # for p, v0, v1, r0, r1

end # function print_pairwise_policy


function sharray2array(sharray::SharedArray{Float64, 2})
    result = zeros(sharray.dims)
    for i = 1:sharray.dims[1]
        for j = 1:sharray.dims[2]
            result[i, j] = sharray[i, j]
        end # for j
    end # for i
    return result
end # function sharray2array

end # module PilotSCAViz