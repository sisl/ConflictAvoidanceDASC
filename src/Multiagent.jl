#=
Description
    Module for multi-agent UAV algorithm and simulations.
=#

module Multiagent

using DoubleUAVs, POMDPs, UAVs, GridInterpolations, QMDP, Interact,
      PGFPlots, HDF5, JLD

export read_alpha, get_grid, viz_policy, randuavs, stress!, plot_trajs,
       maxsum, maxmin, utm, uncrd, coord, naive, distr, bulk_test


const ALPHA_FILE = "../data/alpha.jld"
const ALPHA_VARIABLE = "alpha"

const X = 1; const XDOT = 2     # x-direction
const Y = 3; const YDOT = 4     # y-direction
const P = 5; const PDOT = 6     # heading
const B = 7; const BDOT = 8     # bank angle

const PROBLEM = DoubleUAV()

const PX = 1         # state variable indexing
const PY = 2
const PP = 3
const PV1 = 4
const PV2 = 5
const NPSTATES = size(PROBLEM.pomdp.states, 1)
const NSTATES = size(PROBLEM.pomdp.states, 2)

const COC = -1
const INDIV_ACTIONS = deg2rad([-20, -10, 0, 10, 20, rad2deg(COC)])
const NINDIV_ACTIONS = length(INDIV_ACTIONS)
actions = zeros(2, NINDIV_ACTIONS^2)
iaction = 1
for ia2 = 1:NINDIV_ACTIONS
    for ia1 = 1:NINDIV_ACTIONS
        actions[1, iaction] = INDIV_ACTIONS[ia1]
        actions[2, iaction] = INDIV_ACTIONS[ia2]
        iaction += 1
    end # for ia1
end # for ia2
const ACTIONS = actions
const NACTIONS = size(ACTIONS, 2)

const DT = 5.0              # [s]
const DTI = 1.0             # [s]
const G = 9.8               # [m/s^2]

const DIM_ACTIONS = 2
const DIM_POMDP_STATES = 5
const DIM_STATES = 8        # TYPE_SIMPLE state dimensionality
const TERM_STATE = 1e5 * ones(DIM_POMDP_STATES)
const TERM_STATE_REAL = 1e5 * ones(DIM_STATES)

const XMIN = -2e3
const XMAX = 2e3
const YMIN = -2e3
const YMAX = 2e3
const PMIN = 0.0
const PMAX = 2 * pi
const VMIN = 10
const VMAX = 20

const WN = 0.15             # [rad/s]
const MIN_SEP = 500.0       # [m]
const MIN_SEP_INIT = 600.0  # [m]

const FEAS_PERCENTILE = 0.8
const CRASH_REWARD = -5e2

const NT = 100              # []
const MAX_ITER_JESP = 100

const MIN_RAD = 2000        # [m]
const MAX_RAD = 3000        # [m]
const RANGE_RAD = MAX_RAD - MIN_RAD

const MU = 0.0              # [rad] gaussian noise mean
const SIGMA = deg2rad(2)    # [rad] gaussian noise variance

const CONSENSUS_REWARD = 0.5

const SIGMA_V = 1.0         # [m/s]
const SIGMA_B = deg2rad(2)  # [rad]
const SIGMA_P = deg2rad(2)  # [rad]
const SIGMA_XY = 50         # [m]

const LARGE_NUMBER = int64(1e9)


type Policy
    alpha :: Matrix{Float64}
    actions :: Matrix{Float64}
    nactions :: Int64
    qvals :: Vector{Float64}
    function Policy(alpha::Matrix{Float64}, actions::Matrix{Float64})
        return new(alpha, actions, size(actions, 2), zeros(size(actions, 2)))
    end # function Policy
end


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


function read_policy(d::DoubleUAV, filename::ASCIIString)
    alpha = load(ALPHA_FILE, ALPHA_VARIABLE)
    return Policy(alpha, d.pomdp.actions)
end # function read_pol


function read_alpha()
    return load(ALPHA_FILE, ALPHA_VARIABLE)
end # function read_alpha


function spdot(v1::SparseMatrixCSC{Float64, Int64},
               v2::SparseMatrixCSC{Float64, Int64})
    sum = 0.0
    if length(v1.rowval) < length(v2.rowval)
        for iv1 in 1:length(v1.rowval)
            sum += v1.nzval[iv1] * v2[v1.rowval[iv1]]
        end # for iv1
    else
        for iv2 in 1:length(v2.rowval)
            sum += v2.nzval[iv2] * v1[v2.rowval[iv2]]
        end # for iv2
    end # if
    return sum
end # function spdot


function spdot(v1::Vector{Float64},
               v2::SparseMatrixCSC{Float64, Int64})
    sum = 0.0
    for iv2 in 1:length(v2.rowval)
        sum += v2.nzval[iv2] * v1[v2.rowval[iv2]]
    end # for iv2
    return sum
end # function spdot


function spdot(m1::Matrix{Float64}, icol::Int64,
               v2::SparseMatrixCSC{Float64, Int64})
    sum = 0.0
    for iv2 in 1:length(v2.rowval)
        sum += v2.nzval[iv2] * m1[v2.rowval[iv2], icol]
    end # for iv2
    return sum
end # function spdot


function spdot(v1::SparseMatrixCSC{Float64, Int64},
               v2::Vector{Float64})
    sum = 0.0
    for iv1 in 1:length(v1.rowval)
        sum += v1.nzval[iv1] * v2[v1.rowval[iv1]]
    end # for iv1
    return sum
end # function spdot


function spdot(v1::SparseMatrixCSC{Float64, Int64},
               m2::Matrix{Float64}, icol::Int64)
    sum = 0.0
    for iv1 in 1:length(v1.rowval)
        sum += v1.nzval[iv1] * m2[v1.rowval[iv1], icol]
    end # for iv1
    return sum
end # function spdot


function randuavs(nuav::Int64)
    return [UAV(VMIN + (VMAX - VMIN) * rand()) for iuav = 1:nuav]
end # function randuavs


function randactions(nuav::Int64)
    return [0.0 for iuav = 1:nuav] # for standardization
end # function randactions


function get_grid(states::Matrix{Float64})
    g = RectangleGrid(sort(unique(vec(states[PX, 1:end - 1]))),
                      sort(unique(vec(states[PY, 1:end - 1]))),
                      sort(unique(vec(states[PP, 1:end - 1]))),
                      sort(unique(vec(states[PV1, 1:end - 1]))),
                      sort(unique(vec(states[PV2, 1:end - 1]))))
    return g
end # function get_grid


function get_pomdp_state(s1::Vector{Float64}, s2::Vector{Float64})
    dx = s2[X] - s1[X]
    dy = s2[Y] - s1[Y]
    
    xr = dx * cos(s1[P]) + dy * sin(s1[P])
    yr = -dx * sin(s1[P]) + dy * cos(s1[P])
    pr = s2[P] - s1[P]
    v1 = sqrt(s1[XDOT]^2 + s1[YDOT]^2)
    v2 = sqrt(s2[XDOT]^2 + s2[YDOT]^2)

    if xr < XMIN || xr > XMAX ||
       yr < YMIN || yr > YMAX ||
       s1 == TERM_STATE_REAL || s2 == TERM_STATE_REAL
        return TERM_STATE
    else
        return [xr, yr, pr, v1, v2]
    end # if
end # function get_pomdp_state


function get_belief(pstate::Vector{Float64}, grid::RectangleGrid)
    belief = spzeros(NSTATES, 1)
    if pstate == TERM_STATE
        belief[end] = 1.0
    else
        indices, weights = interpolants(grid, pstate)
        belief[indices] = weights
    end # if
    return belief
end # function get_belief


function get_belief!(pstate::Vector{Float64}, grid::RectangleGrid, 
                     beliefs::SparseMatrixCSC{Float64, Int64},
                     uavIdx::Int64)
    if pstate == TERM_STATE
        beliefs[end, uavIdx] = 1.0
    else
        indices, weights = interpolants(grid, pstate)
        beliefs[indices, uavIdx] = weights
    end # if
end # function get_belief


function action_idx(a1::Float64, a2::Float64)
    idx1 = 0
    idx2 = 0
    for iaction in 1:NINDIV_ACTIONS
        if INDIV_ACTIONS[iaction] == a1
            idx1 = iaction
        end # if
        if INDIV_ACTIONS[iaction] == a2
            idx2 = iaction
        end # if
        if idx1 != 0 && idx2 != 0
            break
        end # if
    end # for iaction
    return sub2ind([NINDIV_ACTIONS, NINDIV_ACTIONS], idx1, idx2)
end # function action_idx


function and(a1::BitArray{1}, a2::BitArray{1})
    @assert length(a1) == length(a2)
    for i = 1:length(a1)
        if a2[i] == false
            a1[i] = false
        end # if
    end # for i
    return a1
end # function and


function consensus(action::Float64, suggest::Symbol)
    if (suggest == :left && action > 0) ||
       (suggest == :right && action < 0) ||
       (suggest == :straight && action == 0)
        return true
    else # suggest == :anything or no consensus
        return false
    end
end # function consensus


function maxsum(ownship::Int64, uavs::Vector{UAV}, 
                actions::Vector{Float64}, alpha::Matrix{Float64}, 
                grid::RectangleGrid, suggest::Symbol=:anything)
    #=
    Returns the max-sum ownship action and utility value based 
    on utility fusion. 
    =#
    nuav = length(uavs)
    utils = zeros(NINDIV_ACTIONS)

    # compute belief states for each uav pair
    beliefs = spzeros(NSTATES, nuav)
    for iu = 1:nuav
        if iu == ownship
            continue
        end # if
        state = get_pomdp_state(uavs[ownship].state, uavs[iu].state)
        get_belief!(state, grid, beliefs, iu)
    end # for iu

    # find action corresponding to max-sum utility
    for iaction = 1:NINDIV_ACTIONS
        # compute sum of utility corresponding to ownship action
        for iu = 1:nuav
            if iu == ownship
                continue
            end # if
            ia = action_idx(INDIV_ACTIONS[iaction], actions[iu])
            if beliefs[end, iu] == 1
                if INDIV_ACTIONS[iaction] == 0
                    utils[iaction] += 0.0
                else
                    utils[iaction] += -0.5
                end # if
            else
                utils[iaction] += spdot(alpha, ia, beliefs[:, iu])
            end # if
        end # for iu

        # add online cost due to suggestion
        if consensus(INDIV_ACTIONS[iaction], suggest)
            utils[iaction] += CONSENSUS_REWARD * abs(utils[iaction])
        end # if
    end # for iaction

    # compute feasible range of actions
    bestAction = INDIV_ACTIONS[indmax(utils)]
    bestUtil = maximum(utils)
    
    return bestAction, bestUtil
end # function maxsum


function maxmin(ownship::Int64, uavs::Vector{UAV}, 
                actions::Vector{Float64}, alpha::Matrix{Float64}, 
                grid::RectangleGrid, suggest::Symbol=:anything)
    #=
    Returns the max-min ownship action and utility value based 
    on utility fusion. 
    =#
    nuav = length(uavs)
    utils = Inf * ones(NINDIV_ACTIONS)

    # compute belief states for each uav pair
    beliefs = spzeros(NSTATES, nuav)
    for iu = 1:nuav
        if iu == ownship
            continue
        end # if
        state = get_pomdp_state(uavs[ownship].state, uavs[iu].state)
        get_belief!(state, grid, beliefs, iu)
    end # for iu

    # find action corresponding to max-min utility
    for iaction = 1:NINDIV_ACTIONS
        # compute minimum utility corresponding to ownship action
        for iu = 1:nuav
            if iu == ownship
                continue
            end # if
            aIdx = action_idx(INDIV_ACTIONS[iaction], actions[iu])
            if beliefs[end, iu] == 1
                if INDIV_ACTIONS[iaction] == 0
                    util = 0.0
                else
                    util = -0.5
                end # if
            else
                util = spdot(alpha, aIdx, beliefs[:, iu])
            end # if
            if utils[iaction] > util
                utils[iaction] = util
            end # if
        end # for iu

        # add online cost due to suggestion
        if consensus(INDIV_ACTIONS[iaction], suggest)
            utils[iaction] += CONSENSUS_REWARD * abs(utils[iaction])
        end # if
    end # for iaction

    # compute feasible range of actions
    bestAction = INDIV_ACTIONS[indmax(utils)]
    bestUtil = maximum(utils)
    
    return bestAction, bestUtil
end # function maxmin


function maxsum(iaction::Int64, ownship::Int64, uavs::Vector{UAV}, 
                actions::Vector{Float64}, alpha::Matrix{Float64}, 
                grid::RectangleGrid, suggest::Symbol=:anything)
    #=
    Returns the max-sum ownship utility value based on utility 
    fusion. 
    =#
    nuav = length(uavs)
    util = 0
    
    # compute belief states for each uav pair
    beliefs = spzeros(NSTATES, nuav)
    for iu = 1:nuav
        if iu == ownship
            continue
        end # if
        state = get_pomdp_state(uavs[ownship].state, uavs[iu].state)
        get_belief!(state, grid, beliefs, iu)
    end # for iu
    
    # compute sum of utility corresponding to ownship action
    for iu = 1:nuav
        if iu == ownship
            continue
        end # if
        ia = action_idx(INDIV_ACTIONS[iaction], actions[iu])
        if beliefs[end, iu] == 1
            if INDIV_ACTIONS[iaction] == 0
                util += 0.0
            else
                util += -0.5
            end # if
        else
            util += spdot(alpha, ia, beliefs[:, iu])
        end # if
    end # for iu

    # add online cost due to suggestion
    if consensus(INDIV_ACTIONS[iaction], suggest)
        util += CONSENSUS_REWARD * abs(util)
    end # if
    
    return util
end # function maxsum


function maxmin(iaction::Int64, ownship::Int64, uavs::Vector{UAV}, 
                actions::Vector{Float64}, alpha::Matrix{Float64}, 
                grid::RectangleGrid, suggest::Symbol=:anything)
    #=
    Returns the max-min ownship utility value based on utility 
    fusion. 
    =#
    nuav = length(uavs)
    min_util = Inf

    # compute belief states for each uav pair
    beliefs = spzeros(NSTATES, nuav)
    for iu = 1:nuav
        if iu == ownship
            continue
        end # if
        state = get_pomdp_state(uavs[ownship].state, uavs[iu].state)
        get_belief!(state, grid, beliefs, iu)
    end # for iu
   
    # compute minimum utility corresponding to ownship action
    for iu = 1:nuav
        if iu == ownship
            continue
        end # if
        aIdx = action_idx(INDIV_ACTIONS[iaction], actions[iu])
        if beliefs[end, iu] == 1
            if INDIV_ACTIONS[iaction] == 0
                util = 0.0
            else
                util = -0.5
            end # if
        else
            util = spdot(alpha, aIdx, beliefs[:, iu])
        end # if
        if min_util > util
            min_util = util
        end # if
    end # for iu

    # add online cost due to suggestion
    if consensus(INDIV_ACTIONS[iaction], suggest)
        utils[iaction] += CONSENSUS_REWARD * abs(utils[iaction])
    end # if
    
    return min_util
end # function maxmin


function get_distance(u1::Int64, u2::Int64, uavs::Vector{UAV})
    return norm([uavs[u1].state[X] - uavs[u2].state[X], 
                 uavs[u1].state[Y] - uavs[u2].state[Y]])
end # function get_distance


function get_distance(u1::Int64, u2::Int64, uavs::Vector{UAV},
                      quavs::Vector{UAV})
    return norm([uavs[u1].state[X] - quavs[u2].state[X], 
                 uavs[u1].state[Y] - quavs[u2].state[Y]])
end # function get_distance


function greedy(uavs::Vector{UAV}, alpha::Matrix{Float64}, 
                grid::RectangleGrid, utilFn::Function)
    #=
    Returns the uav actions solution via a greedy search
    for the best policy against the closest uav, assuming
    all other intruders are white noise intruders.
    =#
    nuav = length(uavs)
    closest = zeros(Int64, nuav)
    actions = zeros(nuav)
    # ranges = Array(Vector{Float64}, nuav)
    utils = zeros(nuav)

    # find closest other uav
    for iu = nuav:-1:2
        distances = Inf * ones(iu - 1)
        for iu2 = 1:iu - 1
            distances[iu2] = get_distance(iu, iu2, uavs)
        end # for iu2
        closest[iu] = indmin(distances)
        closest[closest[iu]] = iu
    end # for iu

    # compute best actions
    for iu = 1:nuav
        action, util = 
            utilFn(1, [uavs[iu], uavs[closest[iu]]],
                   zeros(2), alpha, grid)
        actions[iu] = action
        utils[iu] = util
    end # for iu
    return actions, utils
end # function greedy


function greedy(uavs::Vector{UAV}, alpha::Matrix{Float64}, 
                grid::RectangleGrid, utilFn::Function,
                quavs::Vector{UAV})
    #=
    Returns the uav actions solution via a greedy search
    for the best policy against the closest uav, assuming
    all other intruders are white noise intruders.
    =#
    nuav = length(uavs)
    closest = zeros(Int64, nuav)
    actions = zeros(nuav)
    # ranges = Array(Vector{Float64}, nuav)
    utils = zeros(nuav)

    # find closest other uav
    for iu = nuav:-1:2
        distances = Inf * ones(iu - 1)
        for iu2 = 1:iu - 1
            distances[iu2] = get_distance(iu, iu2, uavs, quavs)
        end # for iu2
        closest[iu] = indmin(distances)
        closest[closest[iu]] = iu
    end # for iu

    # compute best actions
    for iu = 1:nuav
        action, util = 
            utilFn(1, [uavs[iu], quavs[closest[iu]]], 
                   zeros(2), alpha, grid)
        actions[iu] = action
        utils[iu] = util
    end # for iu
    return actions, utils
end # function greedy


function jesp(uavs::Vector{UAV}, alpha::Matrix{Float64}, 
                  grid::RectangleGrid, utilFn::Function)
    #=
    Returns the utility fusion uav actions solution via
    the joint equilibrium based search for policy method.
    =#
    nuav = length(uavs)
    
    actions = zeros(nuav)
    
    prevActions = Inf * ones(nuav)
    utils = -Inf * ones(nuav)
    
    bestUtil = -Inf
    bestActions = copy(actions)
    
    niter = 0
    solutionImproving = false

    while niter < MAX_ITER_JESP
        netUtil = 0.0
        prevActions = copy(actions)
        
        jespOrder = randperm(nuav)
        for iu = jespOrder
            action, util = 
                utilFn(iu, uavs, actions, alpha, grid)
            utils[iu] = util
            actions[iu] = action
            
            if bestUtil < sum(utils)
                bestUtil = sum(utils)
                bestActions = copy(actions)
                if iu == jespOrder[end]
                    solutionImproving = true
                end # if
            else
                solutionImproving = false
            end # if
        end

        if !solutionImproving
            break
        end # if

        niter += 1
        if niter == MAX_ITER_JESP
            print("Max iter reached, check solution\n")
        end # if
    end # while
    return bestActions, bestUtil, niter
end # function jesp


function fuse_uncrd(uavs::Vector{UAV}, alpha::Matrix{Float64}, 
                    grid::RectangleGrid, utilFn::Function)
    #=
    Returns the utility fusion uav actions solution via
    a greedy search for policy method, assuming all other
    intruders are white noise intruders.
    =#
    nuav = length(uavs)
    dummyActions = zeros(nuav)
    actions = zeros(nuav)
    utils = zeros(nuav)

    for iu = 1:nuav
        action, util = 
            utilFn(iu, uavs, dummyActions, alpha, grid)
        actions[iu] = action
        utils[iu] = util
    end # for iu
    return actions, utils
end # function fuse_uncrd


function fuse_uncrd(quavs::Vector{UAV}, alpha::Matrix{Float64}, 
                    grid::RectangleGrid, utilFn::Function)
    #=
    Returns the utility fusion uav actions solution via
    a greedy search for policy method, assuming all other
    intruders are white noise intruders.
    =#
    nuav = length(quavs)
    dummyActions = zeros(nuav)
    actions = zeros(nuav)
    utils = zeros(nuav)

    for iu = 1:nuav
        iuavs = [quavs[iu], quavs[1:nuav .!= iu]]
        action, util = 
            utilFn(1, iuavs, dummyActions, alpha, grid)
        actions[iu] = action
        utils[iu] = util
    end # for iu
    return actions, utils
end # function fuse_uncrd


function fuse_coord(quavs::Vector{UAV}, alpha::Matrix{Float64}, 
                    grid::RectangleGrid, utilFn::Function)
    #=
    Returns the utility fusion uav actions solution via
    a search for the best policy based on coordination messages
    from other uavs. Note that each uav still assumes that
    all other uavs are white noise intruders.
    =#
    nuav = length(quavs)
    actions = zeros(nuav)
    utils = zeros(nuav)

    # compute what each uav thinks is best for all uavs
    sactions = zeros(nuav, nuav)
    for iu = 1:nuav
        jactions, _ = jesp(quavs, alpha, grid, utilFn)
        sactions[:, iu] = jactions
    end # for iu

    # generate what each uav guesses it should do
    gactions = diag(sactions)

    # compute best action with suggested turning sense
    for iu = 1:nuav
        # average coordination messages to get suggested sense
        jaction = sum(sactions[iu, 1:nuav .!= iu])
        if jaction > 0
            suggest = :left
        elseif jaction < 0
            suggest = :right
        else # jaction == 0
            suggest = :straight
        end # if

        # compute best action
        action, util = 
            utilFn(iu, quavs, gactions, alpha, grid, suggest)
        actions[iu] = action
        utils[iu] = util
    end # for iu
    return actions, utils
end # function fuse_coord


function fuse_distr(quavs::Vector{UAV}, alpha::Matrix{Float64}, 
                    grid::RectangleGrid, utilFn::Function)
    nuav = length(quavs)
    actions = zeros(nuav)
    pactions = Inf * ones(nuav)
    utils = zeros(nuav)

    niter = 0
    while pactions != actions && niter < MAX_ITER_JESP
        niter += 1
        pactions = actions
        # compute each uav's best greedy action in parallel
        for iu = 1:nuav
            action, util = 
                utilFn(iu, quavs, pactions, alpha, grid)
            actions[iu] = action
            utils[iu] = util
        end # for iu
    end # while
    return actions, utils
end # function fuse_distr


function noisify(uavs::Vector{UAV})
    #=
    Adds noise to the state variables of uavs.
    =#
    quavs = deepcopy(uavs)
    nuav = length(quavs)
    for iu = 1:nuav
        quavs[iu].state[X] += SIGMA_XY * randn()
        quavs[iu].state[Y] += SIGMA_XY * randn()
        quavs[iu].state[P] += SIGMA_P * randn()
        speed = SIGMA_V * randn() + 
                sqrt(quavs[iu].state[XDOT]^2 + 
                     quavs[iu].state[YDOT]^2)
        quavs[iu].state[XDOT] = speed * cos(quavs[iu].state[P])
        quavs[iu].state[YDOT] = speed * sin(quavs[iu].state[P])
    end # for iu
    return quavs
end # function noisify


function utm(uavs::Vector{UAV}, alpha::Matrix{Float64}, 
             grid::RectangleGrid, utilFn::Function, noise::Bool)
    #=
    Centralized controller computes the best action for all planes
    based on received measurements from each UAV and sends advisories.
    For UTM-like centralized controller, or ACAS-like coordinated
    controller with perfect information (measurements + messages). 
    =#
    quavs = noisify(uavs)
    tic()
    actions, _ = jesp(quavs, alpha, grid, utilFn)
    if noise
        addnoise!(actions)
    end # if
    return actions
end # function utm


function uncrd(uavs::Vector{UAV}, alpha::Matrix{Float64}, 
               grid::RectangleGrid, utilFn::Function, noise::Bool)
    #=
    Each UAV plans its best action based on white-noise assumption
    for all other uavs using utility fusion (no coordination). For
    ACAS-like uncoordinated controller with quantized state messages.
    =#
    quavs = noisify(uavs)
    actions, _ = fuse_uncrd(quavs, alpha, grid, utilFn)
    if noise
        addnoise!(actions)
    end # if
    return actions
end # function uncrd


function coord(uavs::Vector{UAV}, alpha::Matrix{Float64}, 
               grid::RectangleGrid, utilFn::Function, noise::Bool)
    #=
    Each UAV plans its best action based on coordination messages
    from all other uavs using utility fusion. For ACAS-like
    coordinated controller with quantized state messages.
    =#
    quavs = noisify(uavs)
    actions, _ = fuse_coord(quavs, alpha, grid, utilFn)
    if noise
        addnoise!(actions)
    end # if
    return actions
end # function uncrd


function naive(uavs::Vector{UAV}, alpha::Matrix{Float64}, 
               grid::RectangleGrid, utilFn::Function, noise::Bool)
    #=
    Each UAV plans its best action based on white-noise assumption
    for the closest other uav (no utility fusion; greedy takes
    care of this, even though the api is the same as acas and utm).
    For either UTM or ACAS-like baseline controller.
    =#
    actions, _ = greedy(uavs, alpha, grid, utilFn)
    if noise
        addnoise!(actions)
    end # if
    return actions
end # function naive


function distr(uavs::Vector{UAV}, alpha::Matrix{Float64}, 
               grid::RectangleGrid, utilFn::Function, noise::Bool)
    #=
    Each UAV greedily selects its own action in a distributed
    fashion at each iteration and passes its solution to server.
    =#
    quavs = noisify(uavs)
    actions, _ = fuse_distr(quavs, alpha, grid, utilFn)
    if noise
        addnoise!(actions)
    end # if
    return actions
end # function naive


function set_scenario!(iu::Int64, uavs::Vector{UAV}, 
                       scenario::Vector{Float64})
    uavs[iu].state[X] = scenario[1]
    uavs[iu].state[Y] = scenario[2]
    uavs[iu].state[P] = scenario[3]
    uavs[iu].state[B] = 0.0
    uavs[iu].state[XDOT] = uavs[iu].speed * cos(uavs[iu].state[P])
    uavs[iu].state[YDOT] = uavs[iu].speed * sin(uavs[iu].state[P])
    uavs[iu].state[PDOT] = 0.0
    uavs[iu].state[BDOT] = 0.0
end # function set_scenario


function viz_policy(alpha::Matrix{Float64}, grid::RectangleGrid)
    uavs = randuavs(3)
    set_scenario!(2, uavs, [1250, 600, deg2rad(180)])
    set_scenario!(3, uavs, [1250, -600, deg2rad(180)])

    utilFnRange = [maxsum, maxmin]
    prange = 0:30:360    

    @manipulate for iutilFn = 1:length(utilFnRange), p = prange
        utilFn = utilFnRange[iutilFn]
        function getmap1(x::Float64, y::Float64)
            set_scenario!(1, uavs, [x, y, deg2rad(p)])
            actions, _ = jesp(uavs, alpha, grid, utilFn)
            if abs(actions[1]) > 0.5
                return 0.0
            end # if
            return rad2deg(actions[1])
        end # function getmap1
        
        function getmap2(x::Float64, y::Float64)
            set_scenario!(1, uavs, [x, y, deg2rad(p)])
            actions, _ = jesp(uavs, alpha, grid, utilFn)
            if abs(actions[2]) > 0.5
                return 0.0
            end # if
            return rad2deg(actions[2])
        end # function getmap2

        function getmap3(x::Float64, y::Float64)
            set_scenario!(1, uavs, [x, y, deg2rad(p)])
            actions, _ = jesp(uavs, alpha, grid, utilFn)
            if abs(actions[3]) > 0.5
                return 0.0
            end # if
            return rad2deg(actions[3])
        end # function getmap3

        g = GroupPlot(3, 1, groupStyle="horizontal sep=1.5cm")
        push!(g, Axis([
            Plots.Image(getmap1, (-3000, 3000), (-3000, 3000), 
                        zmin = -20, zmax = 20, xbins = 150, ybins = 150,
                        colormap = ColorMaps.Named("RdBu"), 
                        colorbar = false),
    Plots.Node(L">", 1250, 600, style="rotate=180,font=\\huge"),
    Plots.Node(L">", 1250, -600, style="rotate=180,font=\\huge")], 
                      width="9cm", height="8cm", 
                      xlabel="x (m)", ylabel="y (m)", 
                      title="First aircraft action",
                      style="xtick={-3000,-2000,...,3000}"))
        push!(g, Axis([
            Plots.Image(getmap2, (-3000, 3000), (-3000, 3000), 
                        zmin = -20, zmax = 20, xbins = 150, ybins = 150,
                        colormap = ColorMaps.Named("RdBu"), 
                        colorbar = false),
    Plots.Node(L">", 1250, 600, style="rotate=180,font=\\huge"),
    Plots.Node(L">", 1250, -600, style="rotate=180,font=\\huge")], 
                      width="9cm", height="8cm", 
                      xlabel="x (m)" ,
                      title="Second aircraft action",
                      style="xtick={-3000,-2000,...,3000}"))
        push!(g, Axis([
            Plots.Image(getmap3, (-3000, 3000), (-3000, 3000), 
                        zmin = -20, zmax = 20, xbins = 150, ybins = 150,
                        colormap = ColorMaps.Named("RdBu")),
    Plots.Node(L">", 1250, 600, style="rotate=180,font=\\huge"),
    Plots.Node(L">", 1250, -600, style="rotate=180,font=\\huge")], 
            width="9cm", height="8cm", 
                      xlabel="x (m)",
                      title="Third aircraft action",
                      style="xtick={-3000,-2000,...,3000}"))
        g
    end # for utilFn, p
end # function viz_policy


function norm_angle(angle::Float64)
    return ((angle % (2 * pi)) + 2 * pi) % (2 * pi)
end # function norm_angle


function rand_uavpos()
    # randomly generate position in polar coordinates
    r = MIN_RAD + RANGE_RAD * rand()
    p = norm_angle(2 * pi * rand())

    # convert into cartesian coordinates
    x = r * cos(p)
    y = r * sin(p)
    head = norm_angle(p + pi)
    return [x, y, head]
end # function rand_uavpos


function initconflict(x::Float64, y::Float64, uavs::Vector{UAV})
    for iu = 1:length(uavs)
        if (x - uavs[iu].state[X])^2 + 
           (y - uavs[iu].state[Y])^2 <
           (MIN_SEP_INIT)^2 
            return true
        end # if
    end # for iu
    return false
end # function inconflict


function inconflict(uav::UAV, uavs::Vector{UAV})
    for iu = 1:length(uavs)
        if MIN_SEP^2 > (uav.state[X] - uavs[iu].state[X])^2 + 
                       (uav.state[Y] - uavs[iu].state[Y])^2
            return true
        end # if
    end # for iu
    return false
end # function inconflict
 

function inconflict(uav1::UAV, uav2::UAV)
    if MIN_SEP^2 > (uav1.state[X] - uav2.state[X])^2 + 
                   (uav1.state[Y] - uav2.state[Y])^2
        return true
    else
        return false
    end # if
end # function inconflict


function init_test!(uavs::Vector{UAV})
    for iu = 1:length(uavs)
        newpos = rand_uavpos()
        while initconflict(newpos[1], newpos[2], uavs[1:iu - 1])
            newpos = rand_uavpos()
        end # while
        set_scenario!(iu, uavs, newpos)
    end # for iu
end # function init_test!


function plot_trajs(alg::Function, util::Function, alpha::Matrix{Float64}, g::RectangleGrid)
    trajs, n, c, t = stress!(alg, util, true, NT, randuavs(4), alpha, g)
    print("number of collisions = ", sum(c), ", average decision time = ")
    @printf("%.3e ms", t * 1000)
    PGFPlots.Axis([Plots.Linear(vec(trajs[1, X, :]), vec(trajs[1, Y, :]), mark="none", style="blue,->,smooth,thick"), 
                   Plots.Linear(vec(trajs[2, X, :]), vec(trajs[2, Y, :]), mark="none", style="red,->,smooth,thick"), 
                   Plots.Linear(vec(trajs[3, X, :]), vec(trajs[3, Y, :]), mark="none", style="green,->,smooth,thick"),
                   Plots.Linear(vec(trajs[4, X, :]), vec(trajs[4, Y, :]), mark="none", style="black,->,smooth,thick"),
                   Plots.Linear([trajs[1, X, 1]], [trajs[1, Y, 1]], mark="*", style="blue"),
                   Plots.Linear([trajs[2, X, 1]], [trajs[2, Y, 1]], mark="*", style="red"),
                   Plots.Linear([trajs[3, X, 1]], [trajs[3, Y, 1]], mark="*", style="green"),
                   Plots.Linear([trajs[4, X, 1]], [trajs[4, Y, 1]], mark="*", style="black")],
                  xlabel="x (m)", ylabel="y (m)", width="12cm", height="12cm", 
                  style="axis equal,xtick={-4000,-2000,...,4000}, ytick={-4000,-2000,...,4000}")
end # function plot_traj


function pid_policy(bankAngle::Float64)
    function policy(t::Float64, state::Vector{Float64})
        return 2 * WN * (-state[BDOT]) + WN^2 * (bankAngle - state[B])
    end # function policy
    return policy
end # function pid_policy


function init_trajs(uavs::Vector{UAV}, nt::Int64)
    nuavs = length(uavs)
    trajs = zeros(nuavs, DIM_STATES, nt * int64(DT / DTI) + 1)
    for iu = 1:nuavs
        trajs[iu, :, 1] = uavs[iu].state
    end # for iu
    return trajs
end # function init_trajs


function gen_traj(uav::UAV, policy::Function)
    _, traj = generate_trajectory(uav, policy, 1, DTI)
    return traj
end # function gen_traj


function sim_trajs!(uavs::Vector{UAV}, actions::Vector{Float64},
                    dtcounts::Vector{Int64}, trajs::Array{Float64,3},
                    conflicts::BitArray{2})
    nuavs = length(uavs)
    nlms = 0

    # generate pid policies
    pids = Array(Function, nuavs)
    for iu = 1:nuavs
        # correct COC to 0 degrees nominal turn
        if abs(actions[iu]) > 0.5
            actions[iu] = 0.0
        end # if
        # add noise to bank angle commanded
        pids[iu] = pid_policy(actions[iu] + SIGMA_B * randn())
    end # for iu

    # simulate trajectories
    for idt = 1:int64(DT / DTI)
        # update state for each uav
        for iu = 1:nuavs
            dtcounts[iu] += 1            
            traj = gen_traj(uavs[iu], pids[iu])
            trajs[iu, :, dtcounts[iu]] = traj[:, end]
            set_state!(uavs[iu], traj[:, end])
        end # for iu

        # check conflict for each uav pair
        for iu1 = 1:nuavs
            for iu2 = iu1 + 1:nuavs
                if inconflict(uavs[iu1], uavs[iu2])
                    nlms += 1
                    conflicts[iu1, iu2] = true
                end # if
            end # for iu2
        end # for iu1
    end # for idt
    
    return nlms
end # function sim_trajs!


function gaussnoise(mu::Float64=MU, sigma::Float64=SIGMA)
    return mu + sigma * randn()
end # function gaussnoise


function addnoise!(actions::Vector{Float64}, 
                   noise::Function=gaussnoise)
    for iu = 1:length(actions)
        actions[iu] += noise()
    end # for iu
end # function addnoise!


function stress!(policy::Function, utilFn::Function, noise::Bool, nt::Int64, 
                 uavs::Vector{UAV}, alpha::Matrix{Float64}, grid::RectangleGrid)
    init_test!(uavs)
    nuav = length(uavs)
    trajs = init_trajs(uavs, nt)
    nlms = 0
    dtcounts = ones(Int64, nuav)
    conflicts = falses(nuav, nuav)
    avgtime = 0

    for it = 1:nt
        tic()
        actions = policy(uavs, alpha, grid, utilFn, noise)
        avgtime += toq()
        nlms += sim_trajs!(uavs, actions, dtcounts, trajs, conflicts)
    end # for it
    avgtime /= nt

    return trajs, nlms, conflicts, avgtime
end # function stress!

function bulk_test(nuavs::UnitRange{Int64}, nsim::Int64, 
                   alpha::Matrix{Float64}, g::RectangleGrid,
                   nbulk::Int64=1)
    for ibulk = 1:nbulk
        print("beginning stress tests...\n")
        maxminNLMS = zeros(Int64, 5, length(nuavs))
        maxsumNLMS = zeros(Int64, 5, length(nuavs))
        maxminNLMSBool = zeros(Int64, 5, length(nuavs))
        maxsumNLMSBool = zeros(Int64, 5, length(nuavs))
        maxminTimes = zeros(Float64, 5, length(nuavs))
        maxsumTimes = zeros(Float64, 5, length(nuavs))

        for inu = 1:length(nuavs)
            tic()
            nu = nuavs[inu]
            uavs = randuavs(nu)
            
            for isim = 1:nsim
                _, n1, c1, t1 = stress!(utm, maxmin, true, NT, uavs, alpha, g)
                _, n2, c2, t2 = stress!(uncrd, maxmin, true, NT, uavs, alpha, g)
                _, n3, c3, t3 = stress!(coord, maxmin, true, NT, uavs, alpha, g)
                _, n4, c4, t4 = stress!(naive, maxmin, true, NT, uavs, alpha, g)
                _, n5, c5, t5 = stress!(distr, maxmin, true, NT, uavs, alpha, g)
                maxminNLMS[:, inu] += [n1, n2, n3, n4, n5]
                maxminNLMSBool[:, inu] += [sum(c1), sum(c2), sum(c3), sum(c4), sum(c5)]
                maxminTimes[:, inu] += [t1, t2, t3, t4, t5]
            end # for isim
            
            for isim = 1:nsim
                _, n1, c1, t1 = stress!(utm, maxsum, true, NT, uavs, alpha, g)
                _, n2, c2, t2 = stress!(uncrd, maxsum, true, NT, uavs, alpha, g)
                _, n3, c3, t3 = stress!(coord, maxsum, true, NT, uavs, alpha, g)
                _, n4, c4, t4 = stress!(naive, maxsum, true, NT, uavs, alpha, g)
                _, n5, c5, t5 = stress!(distr, maxsum, true, NT, uavs, alpha, g)
                maxsumNLMS[:, inu] += [n1, n2, n3, n4, n5]
                maxsumNLMSBool[:, inu] += [sum(c1), sum(c2), sum(c3), sum(c4), sum(c5)]
                maxsumTimes[:, inu] += [t1, t2, t3, t4, t5]
            end # for isim

            maxminTimes /= nsim
            maxsumTimes /= nsim
            @printf("nuavs = %d: cputime = %.2e sec\n", nu, toq())
            print("maxmin nlms: ", vec(maxminNLMS[:, inu]), "\n")
            print("maxmin nlms bool: ", vec(maxminNLMSBool[:, inu]), "\n")
            print("average maxmin times: ", vec(maxminTimes[:, inu]), "\n")
            print("maxsum nlms: ", vec(maxsumNLMS[:, inu]), "\n")
            print("maxsum nlms bool: ", vec(maxsumNLMSBool[:, inu]), "\n")
            print("average maxsum times: ", vec(maxsumTimes[:, inu]), "\n")
        end # for inu

        outfile = string("../data/results-", ibulk, ".jld")
        save(outfile, "maxminNLMS", maxminNLMS, "maxsumNLMS", maxsumNLMS, 
             "maxminNLMSBool", maxminNLMSBool, "maxsumNLMSBool", maxsumNLMSBool,
             "maxminTimes", maxminTimes, "maxsumTimes", maxsumTimes)
        println("Results saved to ", outfile)
    end # for ibulk
end # function bulk_test

end # module Multiagent