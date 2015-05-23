module Pairwise

using DoubleUAVs, QMDP, POMDPs, GridInterpolations, Interact, 
      PGFPlots, HDF5, JLD

export get_pomdp, gen_pairwise_policy, viz_pairwise_policy


const ALPHA_FILE = "../data/alpha.jld"
const ALPHA_VARIABLE = "alpha"

const DIM_POMDP_STATES = 5
const DIM_STATES = 8
const TERM_STATE_VAR = 1e5
const TERM_STATE = TERM_STATE_VAR * ones(DIM_POMDP_STATES)
const TERM_STATE_REAL = TERM_STATE_VAR * ones(DIM_STATES)

const XMIN = -2e3
const XMAX = 2e3
const YMIN = -2e3
const YMAX = 2e3
const PMIN = 0.0
const PMAX = 2 * pi


type Policy
    alpha       :: Matrix{Float64}
    actions     :: Matrix{Float64}
    nactions    :: Int64
    qvals       :: Vector{Float64}

    function Policy(alpha::Matrix{Float64}, actions::Matrix{Float64})
        return new(alpha, actions, size(actions, 2), zeros(size(actions, 2)))
    end # function Policy
end


function get_pomdp()
    return DoubleUAV()
end # function get_pomdp


function gen_pairwise_policy(d::DoubleUAV, lambda::Float64, 
                             filetag::UTF8String, saveAlpha::Bool=true)
    alpha = qmdp(d, lambda)
    if saveAlpha
        @printf("Writing alpha vector to %s...", ALPHA_FILE)
        tic()
        filename = string("../data/alpha-", filetag, ".jld")
        save(filename, ALPHA_VARIABLE, alpha)
        @printf("done in %.2e.\n\n", toq())
    end # if
    return Policy(alpha, d.pomdp.actions)
end # function gen_policy


function get_pomdp_state(state1::Vector{Float64}, state2::Vector{Float64})
    x1 = state1[X]
    x2 = state2[X]
    y1 = state1[Y]
    y2 = state2[Y]
    p1 = state1[P]
    p2 = state2[P]

    v1 = sqrt(state1[XDOT]^2 + state1[YDOT]^2)
    v2 = sqrt(state2[XDOT]^2 + state2[YDOT]^2)

    dx = x2 - x1
    dy = y2 - y1
    
    xr = dx * cos(p1) + dy * sin(p1)
    yr = -dx * sin(p1) + dy * cos(p1)
    pr = p2 - p1
    
    if xr < XMIN || xr > XMAX ||
       yr < YMIN || yr > YMAX ||
       state1[1] == TERM_STATE_VAR ||
       state2[1] == TERM_STATE_VAR
        return TERM_STATE
    else
        return [xr, yr, pr, v1, v2]
    end # if
end # function get_pomdp_state


function get_belief(pstate::Vector{Float64}, grid::RectangleGrid)
    belief = spzeros(NSTATES, 1)
    if pstate[1] == TERM_STATE_VAR
        belief[end] = 1.0
    else
        indices, weights = interpolants(grid, pstate[1:3])
        for i = 1:length(indices)
            belief[indices[i]] = weights[i]
        end # for i
    end # if
    return belief
end # function get_belief


function get_qval!(policy::Policy, belief::SparseMatrixCSC{Float64, Int64})
    fill!(policy.qvals, 0.0)
    for iaction in 1:policy.nactions
        for ib in 1:length(belief.rowval)
            policy.qvals[iaction] += belief.nzval[ib] * policy.alpha[belief.rowval[ib], iaction]
        end # for b
    end # for iaction
end # function get_qval!


function evaluate(policy::Policy, belief::SparseMatrixCSC{Float64,Int64})
    fill!(policy.qvals, 0.0)
    get_qval!(policy, belief)
    ibest = indmax(policy.qvals)
    return policy.actions[:, ibest], ibest
end # function evaluate


function read_policy(actions::Matrix{Float64})
    alpha = load(ALPHA_FILE, ALPHA_VARIABLE)
    return Policy(alpha, actions)
end # function read_policy


function viz_pairwise_policy(d::DoubleUAV)
    states = d.pomdp.states
    actions = d.pomdp.actions

    xs = sort(unique(vec(states[1, 1:end - 1])))
    ys = sort(unique(vec(states[2, 1:end - 1])))
    ps = sort(unique(vec(states[3, 1:end - 1])))
    vs = sort(unique(vec(states[4, 1:end - 1])))
    ps_deg = rad2deg(ps)
    grid = RectangleGrid(xs, ys, ps, vs, vs)

    pstart = ps_deg[1]
    pend = ps_deg[end]
    pdiv = ps_deg[2] - ps_deg[1]
    
    vstart = vs[1]
    vend = vs[end]
    vdiv = vs[2] - vs[1]

    policy = read_policy(actions)

    @manipulate for p = pstart:pdiv:pend, 
                    v0 = vstart:vdiv:vend, 
                    v1 = vstart:vdiv:vend
        # ownship uav
        function get_heat1(x::Float64, y::Float64)
            action, _ = evaluate(policy, get_belief(
                [x, y, deg2rad(p), v0, v1], grid
            ))
            if abs(action[1]) > 0.5
                return -2.0
            end # if
            return rad2deg(action[1])
        end # function get_heat1
        
        # intruder uav
        function get_heat2(x::Float64, y::Float64)
            action, _ = evaluate(policy, get_belief(
                [x, y, deg2rad(p), v0, v1], grid
            ))
            if abs(action[2]) > 0.5
                return -2.0
            end # if
            return rad2deg(action[2])
        end # function get_heat2
        
        g = GroupPlot(2, 1, groupStyle = "horizontal sep=3cm")
        push!(g, Axis([
            Plots.Image(get_heat1, (int(XMIN), int(XMAX)), 
                        (int(YMIN), int(YMAX)), 
                        zmin = -20, zmax = 20,
                        xbins = 150, ybins = 150,
                        colormap = ColorMaps.Named("RdBu"), colorbar = false),
            Plots.Node(L">", 0, 0, style="rotate=0,font=\\Huge"),
            Plots.Node(L">", 1800, 1800, style=string("rotate=", p, ",font=\\Huge"))
            ], width="10cm", height="10cm", xlabel="x (m)", ylabel="y (m)", title="Ownship action"))
        push!(g, Axis([
            Plots.Image(get_heat2, (int(XMIN), int(XMAX)), 
                        (int(YMIN), int(YMAX)), 
                        zmin = -20, zmax = 20,
                        xbins = 150, ybins = 150,
                        colormap = ColorMaps.Named("RdBu")),
            Plots.Node(L">", 0, 0, style="rotate=0,font=\\Huge"),
            Plots.Node(L">", 1800, 1800, style=string("rotate=", p, ",font=\\Huge"))],
            width="10cm", height="10cm", xlabel="x (m)", title="Intruder action"))
        g
    end # for p, v0, v1
end # function viz_pairwise_policy

end # module Pairwise