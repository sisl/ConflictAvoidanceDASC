#=
Description
    Module for UAV type objects.

    state: [ x,   xdot,       // x-direction
             y,   ydot,       // y-direction
             psi, psidot,     // heading
             phi, phidot  ]   // bank angle
    input:   phiddot          // bank angle 'acceleration'
=#

module UAVs

using ODE

export UAV, generate_trajectory, set_state!


const SPEED_DEF = 10.0
const MU = 0.0
const SIGMA = 0.05
const G = 9.8
const NSTATES = 8
const NINPUTS = 1


type UAV
    state       :: Vector{Float64}
    nStates     :: Int64
    nInputs     :: Int64
    speed       :: Float64

    function UAV(speed::Float64 = SPEED_DEF, 
                 nStates::Int64 = NSTATES, 
                 nInputs::Int64 = NINPUTS)
        return new(zeros(nStates), nStates, nInputs, speed)
    end # function UAV
end # type UAV


function dynamics(uav::UAV, policy::Function, nT::Int64, dT::Float64)
    V = uav.speed
    tStart = 0
    tEnd = tStart + dT * nT
    tVec = tStart:dT:tEnd

    function dxdt(t::Float64, state::Vector{Float64})
        dx = zeros(size(state))
        dx[1] = V * cos(state[5])                   # xdot
        dx[3] = V * sin(state[5])                   # ydot
        dx[5] = G * tan(state[7]) / V               # psidot
        dx[7] = state[8]                            # phidot
        dx[8] = policy(t, state)                    # phiddot
        dx[2] = -V * sin(state[5]) * dx[7]          # xddot
        dx[4] = V * cos(state[5]) * dx[7]           # yddot
        dx[6] = G / V * sec(state[8])^2 * dx[8]     # psiddot
        return dx
    end # function dxdt

    X = zeros(uav.nStates, nT + 1)
    state = uav.state
    X[:, 1] = state
    for ti = 1:length(tVec) - 1
        tSpan = [tVec[ti], tVec[ti + 1]]
        t, xRaw = ode45(dxdt, state, tSpan) # non-stiff
        # t, xRaw = ode23s(dxdt, state, tSpan) # stiff
        state = xRaw[end]
        X[:, ti + 1] = state
    end # for ti
    return tVec, X
end # function discrete_dynamics


#=
Description
    Returns the time-state trajectory pair for the given input
    policy |policy|, number of time steps |nT| and discrete time 
    period |dT|.
=#
function generate_trajectory(uav::UAV, policy::Function, nT::Int64, dT::Float64)
    if uav.nInputs != size(policy(0.0, uav.state), 1)
        error("Dimension mismatch between uav.nInputs and policy output")
    end # if
    T, X = dynamics(uav, policy, nT, dT)
    return T, X
end # function generate_trajectory


function set_state!(uav::UAV, newState::Vector{Float64})
    if uav.nStates == size(newState, 1)
        uav.state = newState
    else
        error("Dimension mismatch between uav.nStates and newState")
    end # if
end # function set_state!

end # module UAVs