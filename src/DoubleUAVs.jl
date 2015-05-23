#=
Description
    Module for pairwise UAV encounters. The dynamics implemented are 
    based on the simple kinematic model, as seen in the UAVs module.
    You can modify the discretization scheme by varying XDIM, YDIM, 
    PDIM, and VDIM. The current values are chosen for a reasonable
    compute time for quick visualization purposes.
=#

module DoubleUAVs

using POMDPs

export NSTATES, COC, DoubleUAV, get_next_state, gen_reward


const DT = 5.0              # [s]
const DTI = 1.0             # [s]
const G = 9.8               # [m/s^2]

const MIN_SEP_SQ = 500.0^2  # [m]
const RBF_STD_DIST = 500.0  # [m]
const RBF_INV_VAR = 1 / RBF_STD_DIST^2   # [1/m^2]

const PEN_ACTION = 0.02
const PEN_CLOSENESS = 10
const PEN_MIN_SEP = 1000
const PEN_CONFLICT = 1.0

const STATE_DIM = 5
const ACTION_DIM = 2

const X = 1                 # [m] relative x-position
const Y = 2                 # [m] relative y-position
const P = 3                 # [rad] relative heading
const V1 = 4                # [m/s] ac1 speed
const V2 = 5                # [m/s] ac2 speed

const XDIM = 51
const XMIN = -2e3
const XMAX = 2e3

const YDIM = 51
const YMIN = -2e3
const YMAX = 2e3

const PDIM = 37
const PMIN = 0.0
const PMAX = 2 * pi

const VDIM = 3
const VMIN = 10
const VMAX = 20

xs = linspace(XMIN, XMAX, XDIM)
ys = linspace(YMIN, YMAX, YDIM)
ps = linspace(PMIN, PMAX, PDIM)
vs = linspace(VMIN, VMAX, VDIM)

const NSTATES = XDIM * YDIM * PDIM * VDIM^2 + 1
tmp_states = zeros(STATE_DIM, NSTATES)
istate = 1
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
                    istate += 1
                end # for ix
            end # for iy
        end # for ip
    end # for iv1
end # for iv2

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


type DoubleUAV

    pomdp::     POMDP
    sigmas::    Matrix{Float64}
    weights::   Vector{Float64}
    
    function DoubleUAV()
        p = POMDP(STATES_OBSERVATIONS, ACTIONS, STATES_OBSERVATIONS, TERM_STATE)
        return new(p, SIGMAS, WEIGHTS)
    end # function DoubleUAV

end # DoubleUAV


function norm_angle(angle::Float64)
    return ((angle % (2 * pi)) + 2 * pi) % (2 * pi)
end # function norm_angle


function gen_next_state(dt::Float64 = DT)
    function get_next_state(state::Vector{Float64}, iaction::Int64, isigma::Int64=1)
        if state[1] == TERM_STATE_VAR
            return TERM_STATE
        else
            x = state[1]
            y = state[2]
            p2 = state[3]
            v1 = state[4]
            v2 = state[5]

            b1 = ACTIONS[1, iaction]
            b2 = ACTIONS[2, iaction]

            # for sigma-point sampling
            if b1 != COC
                b1 = b1 + SIGMAS[1, isigma]
            else 
                b1 = sign(SIGMAS[1, isigma]) * SIGMA_B_COC
            end # if

            if b2 != COC
                b2 = b2 + SIGMAS[2, isigma]
            else 
                b2 = sign(SIGMAS[2, isigma]) * SIGMA_B_COC
            end # if

            v1 = v1 + SIGMAS[3, isigma]
            v2 = v2 + SIGMAS[4, isigma]
        
            if b1 == 0 || b2 == 0 # straight-line path(s)
                if b2 != 0 # uav1 straight path
                    dp2 = dt * G * tan(b2) / v2
                    r2 = abs(v2^2 / (G * tan(b2)))
                    xr = x + r2 * sign(dp2) * (sin(p2) - sin(p2 - dp2)) - v1 * dt
                    yr = y + r2 * sign(dp2) * (-cos(p2) + cos(p2 - dp2))
                    pr = norm_angle(p2 + dp2)
                    
                    if xr < XMIN || xr > XMAX ||
                       yr < YMIN || yr > YMAX
                        return TERM_STATE
                    else
                        return [xr, yr, pr, v1, v2]
                    end # if
                elseif b1 != 0 # uav2 straight path
                    dp1 = dt * G * tan(b1) / v1
                    r1 = abs(v1^2 / (G * tan(b1)))
                    x = x + v2 * dt * cos(p2) - r1 * sign(dp1) * sin(dp1)
                    y = y + v2 * dt * sin(p2) - r1 * sign(dp1) * (cos(dp1) - 1)
                    xr = x * cos(dp1) + y * sin(dp1)
                    yr = -x * sin(dp1) + y * cos(dp1)
                    pr = norm_angle(p2 - dp1)
                    
                    if xr < XMIN || xr > XMAX ||
                       yr < YMIN || yr > YMAX
                        return TERM_STATE
                    else
                        return [xr, yr, pr, v1, v2]
                    end # if
                else # both straight paths
                    xr = x + v2 * dt * cos(p2) - v1 * dt
                    yr = y + v2 * dt * sin(p2)
                    pr = norm_angle(p2)
                    
                    if xr < XMIN || xr > XMAX ||
                       yr < YMIN || yr > YMAX
                        return TERM_STATE
                    else
                        return [xr, yr, pr, v1, v2]
                    end # if
                end # if
            else # both curved paths
                dp1 = dt * G * tan(b1) / v1
                dp2 = dt * G * tan(b2) / v2
                r1 = abs(v1^2 / (G * tan(b1)))
                r2 = abs(v2^2 / (G * tan(b2)))
                x = x + r2 * sign(dp2) * (sin(p2) - sin(p2 - dp2)) - r1 * sign(dp1) * sin(dp1)
                y = y + r2 * sign(dp2) * (-cos(p2) + cos(p2 - dp2)) - r1 * sign(dp1) * (-1 + cos(dp1))
                xr = x * cos(dp1) + y * sin(dp1)
                yr = -x * sin(dp1) + y * cos(dp1)
                pr = norm_angle(p2 + dp2 - dp1)

                if xr < XMIN || xr > XMAX ||
                   yr < YMIN || yr > YMAX
                    return TERM_STATE
                else
                    return [xr, yr, pr, v1, v2]
                end # if
            end # if
        end # if
    end # function get_next_state
    return get_next_state
end # function gen_next_state


function gen_reward(pen_conflict::Float64=PEN_CONFLICT)
    get_next_state = gen_next_state(DTI)
    function get_reward(state::Vector{Float64}, iaction::Int64)
        action = ACTIONS[:, iaction]

        reward = 0
        
        if action[1] == COC
            action[1] = 0.0
        else
            reward = -pen_conflict
        end # if

        if action[2] == COC
            action[2] = 0.0
        else
            reward = reward - pen_conflict
        end # if

        if state[1] == TERM_STATE_VAR
            return reward - PEN_ACTION * rad2deg(norm(action))^2
        else
            reward = reward - PEN_ACTION * rad2deg(norm(action))^2
            
            minSepSq = Inf
            for ti = 1:DT / DTI
                minSepSq = min(minSepSq, norm(state)^2)
                state = get_next_state(state, iaction)
            end # for ti

            if minSepSq < MIN_SEP_SQ
                reward = reward - PEN_MIN_SEP
            end # if

            reward = reward - PEN_CLOSENESS * exp(-minSepSq * RBF_INV_VAR)
            return reward
        end # if
    end # function get_reward
    return get_reward
end # function gen_reward


get_next_state = gen_next_state(DT)

end # module DoubleUAVs