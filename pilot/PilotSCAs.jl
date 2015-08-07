module PilotSCAs

export SCA
export states, actions
export numStates, numActions
export reward, nextStates

export State, Action

using GridInterpolations, DiscreteMDPs

import DiscreteMDPs.DiscreteMDP
import DiscreteMDPs.reward
import DiscreteMDPs.nextStates
import DiscreteMDPs.states
import DiscreteMDPs.actions
import DiscreteMDPs.numStates
import DiscreteMDPs.numActions

using SCAConst, SCAIterators

type SCA <: DiscreteMDP
    
    nStates::Int64
    nActions::Int64
    states::StateIterator
    actions::ActionIterator
    grid::RectangleGrid
    
    function SCA()
        
        states = StateIterator(Xs, Ys, Bearings, Speeds, Speeds)
        actions = ActionIterator(Actions)
        grid = RectangleGrid(Xs, Ys, Bearings, Speeds, Speeds)
        
        return new(NStates, NActions, states, actions, grid)

    end # function SCA
    
end # type SCA


type State
    
    x::Float64
    y::Float64
    bearing::Float64
    speedOwnship::Float64
    speedIntruder::Float64
    clearOfConflict::Bool
    
end # type State


type Action
    
    ownship::Symbol
    intruder::Symbol
    
end # type Action


# Returns an interator over the states.
function states(mdp::SCA)

    return mdp.states
    
end # function states


# Returns an iterator over the actions.
function actions(mdp::SCA)

    return mdp.actions
    
end # function actions


function numStates(mdp::SCA)
    
    return mdp.nStates
    
end # function numStates


function numActions(mdp::SCA)

    return mdp.nActions
    
end # function numActions


function reward(mdp::SCA, istate::Int64, iaction::Int64)

    state = State(0.0, 0.0, 0.0, 10.0, 10.0, true)
    if istate < mdp.nStates
        state = gridState2state(ind2x(mdp.grid, istate))
    end # if

    action = ind2a(mdp.actions.actions, iaction)

    return reward(mdp, state, action)

end # function reward


function reward(mdp::SCA, state::State, action::Action)
    
    reward = 0.0
    
    if action.ownship != :clearOfConflict
        reward -= PenConflict
    end # if
    
    if action.intruder != :clearOfConflict
        reward -= PenConflict
    end # if
    
    turnOwnship = getTurnAngle(action.ownship)
    turnIntruder = getTurnAngle(action.intruder)
    reward -= PenAction * (turnOwnship^2 + turnIntruder^2)
    
    if !state.clearOfConflict
        minSepSq = Inf
        for ti = 1:DT / DTI
            minSepSq = min(minSepSq, getSepSq(state))
            state = getNextState(state, action, DTI)
        end # for ti
        
        if minSepSq < MinSepSq
            reward -= PenMinSep
        end # if
        
        reward -= PenCloseness * exp(-minSepSq * InvVar)
    end # if
        
    return reward
    
end # function reward


function ind2a(actions::Vector{Symbol}, iaction::Int64)

    iOwnship = iaction % length(actions)
    if iOwnship == 0
        iOwnship = 6
    end # if

    iIntruder = (iaction - iOwnship) / length(actions) + 1

    return Action(actions[iOwnship], actions[iIntruder])

end # function ind2a


# Returns turn angle corresponding to action in degrees.
function getTurnAngle(action::Symbol)
    
    if action == :clearOfConflict
        return 0.0
    elseif action == :straight
        return 0.0
    elseif action == :left10
        return 10.0
    elseif action == :right10
        return -10.0
    elseif action == :left20
        return 20.0
    elseif action == :right20
        return -20.0
    else
        throw(ArgumentError("illegal action symbol"))
    end # if
    
end # function getTurnAngle


function getSepSq(state::State)
    
    if state.clearOfConflict
        return Inf
    else
        return state.x^2 + state.y^2
    end # if
    
end # function getSepSq


function getNextState(
        state::State,
        action::Action,
        sigmaTurnOwnship::Float64 = 0.0,
        sigmaTurnIntruder::Float64 = 0.0,
        dt::Float64 = DT)
    
    newState = deepcopy(state)
    
    if !state.clearOfConflict
        
        turnOwnship = deg2rad(getTurnAngle(action.ownship)) + sigmaTurnOwnship
        turnIntruder = deg2rad(getTurnAngle(action.intruder)) + sigmaTurnIntruder

        if turnOwnship == 0.0 || turnIntruder == 0.0  # straight line path(s)
            
            if turnIntruder != 0.0  # ownship straight path
                
                gtan = G * tan(turnIntruder)
                bearingChange = dt * gtan / state.speedIntruder
                radiusIntruder = abs(state.speedIntruder^2 / gtan)
                
                newX = state.x + radiusIntruder * sign(bearingChange) * (sin(state.bearing) - sin(state.bearing - bearingChange)) - state.speedOwnship * dt
                newY = state.y + radiusIntruder * sign(bearingChange) * (-cos(state.bearing) + cos(state.bearing - bearingChange))
                newBearing = norm_angle(state.bearing + bearingChange)

                if newX < Xmin || newX > Xmax || newY < Ymin || newY > Ymax
                    newState.clearOfConflict = true
                else
                    newState.x = newX
                    newState.y = newY
                    newState.bearing = newBearing
                end # if
                
            elseif turnOwnship != 0.0  # intruder straight path

                gtan = G * tan(turnOwnship)
                bearingChange = dt * gtan / state.speedOwnship
                radiusOwnship = abs(state.speedOwnship^2 / gtan)
                
                x = state.x + state.speedIntruder * dt * cos(state.bearing) - radiusOwnship * sign(bearingChange) * sin(bearingChange)
                y = state.y + state.speedIntruder * dt * sin(state.bearing) - radiusOwnship * sign(bearingChange) * (cos(bearingChange) - 1)
                
                newX = x * cos(bearingChange) + y * sin(bearingChange)
                newY = -x * sin(bearingChange) + y * cos(bearingChange)
                newBearing = norm_angle(state.bearing - bearingChange)

                if newX < Xmin || newX > Xmax || newY < Ymin || newY > Ymax
                    newState.clearOfConflict = true
                else
                    newState.x = newX
                    newState.y = newY
                    newState.bearing = newBearing
                end # if
                
            else  # both straight paths

                newX = state.x + state.speedIntruder * dt * cos(state.bearing) - state.speedOwnship * dt
                newY = state.y + state.speedIntruder * dt * sin(state.bearing)
                newBearing = norm_angle(state.bearing)

                if newX < Xmin || newX > Xmax || newY < Ymin || newY > Ymax
                    newState.clearOfConflict = true
                else
                    newState.x = newX
                    newState.y = newY
                    newState.bearing = newBearing
                end # if
                
            end # if
            
        else  # both curved paths
            
            gtanOwnship = G * tan(turnOwnship)
            gtanIntruder = G * tan(turnIntruder)
            
            bearingChangeOwnship = dt * gtanOwnship / state.speedOwnship
            bearingChangeIntruder = dt * gtanIntruder / state.speedIntruder
            
            radiusOwnship = abs(state.speedOwnship^2 / gtanOwnship)
            radiusIntruder = abs(state.speedIntruder^2 / gtanIntruder)
            
            x = state.x + radiusIntruder * sign(bearingChangeIntruder) * (sin(state.bearing) - sin(state.bearing - bearingChangeIntruder))
              - radiusOwnship * sign(bearingChangeOwnship) * sin(bearingChangeOwnship)
            y = state.y + radiusIntruder * sign(bearingChangeIntruder) * (-cos(state.bearing) + cos(state.bearing - bearingChangeIntruder))
              - radiusOwnship * sign(bearingChangeOwnship) * (-1 + cos(bearingChangeOwnship))
            
            newX = x * cos(bearingChangeOwnship) + y * sin(bearingChangeOwnship)
            newY = -x * sin(bearingChangeOwnship) + y * cos(bearingChangeOwnship)
            newBearing = norm_angle(state.bearing + bearingChangeIntruder - bearingChangeOwnship)

            if newX < Xmin || newX > Xmax || newY < Ymin || newY > Ymax
                newState.clearOfConflict = true
            else
                newState.x = newX
                newState.y = newY
                newState.bearing = newBearing
            end # if
            
        end # if
        
    end # if
    
    return newState
    
end # function getNextState


function norm_angle(angle::Float64)
    return ((angle % (2 * pi)) + 2 * pi) % (2 * pi)
end # function norm_angle


# Returns next states and associated transition probabilities.
function nextStates(mdp::SCA, istate::Int64, iaction::Int64)

    if istate == mdp.nStates
        return [mdp.nStates], [1.0]
    end # if

    state = gridState2state(ind2x(mdp.grid, istate))
    action = ind2a(mdp.actions.actions, iaction)

    # sigma sampling
    nominalIndices, nominalProbs = nextStatesSigma(mdp, state, action)
    speedIndices, speedProbs = sigmaSpeed(mdp, state, action)
    bankIndices, bankProbs = sigmaBank(mdp, state, action)
    
    return [
            nominalIndices,
            speedIndices,
            bankIndices], 
        [
            nominalProbs * SigmaWeightNominal,
            speedProbs * SigmaWeightOffNominal,
            bankProbs * SigmaWeightOffNominal]

end # function nextStates


function sigmaSpeed(mdp::SCA, state::State, action::Action)

    # negative sigma
    state.speedOwnship -= SigmaSpeed
    negIndicesOwnship, negProbsOwnship = nextStatesSigma(mdp, state, action)

    # positive sigma
    state.speedOwnship += 2 * SigmaSpeed
    posIndicesOwnship, posProbsOwnship = nextStatesSigma(mdp, state, action)

    # restore original
    state.speedOwnship -= SigmaSpeed

    # negative sigma
    state.speedIntruder -= SigmaSpeed
    negIndicesIntruder, negProbsIntruder = nextStatesSigma(mdp, state, action)

    # positive sigma
    state.speedIntruder += 2 * SigmaSpeed
    posIndicesIntruder, posProbsIntruder = nextStatesSigma(mdp, state, action)

    # restore original
    state.speedIntruder -= SigmaSpeed

    return [negIndicesOwnship, posIndicesOwnship, negIndicesIntruder, posIndicesIntruder],
           [negProbsOwnship, posProbsOwnship, negProbsIntruder, posProbsIntruder]

end # function sigmaSpeed


function sigmaBank(mdp::SCA, state::State, action::Action)

    sigmaBankVal = SigmaBank
    
    if action.ownship == :clearOfConflict
        sigmaBankVal = SigmaBankCOC
    end # if

    # negative sigma
    negIndicesOwnship, negProbsOwnship = nextStatesSigmaAction(mdp, state, action, -sigmaBankVal, 0.0)

    # positive sigma
    posIndicesOwnship, posProbsOwnship = nextStatesSigmaAction(mdp, state, action, sigmaBankVal, 0.0)

    if action.intruder == :clearOfConflict
        sigmaBankVal = SigmaBankCOC
    end # if

    # negative sigma
    negIndicesIntruder, negProbsIntruder = nextStatesSigmaAction(mdp, state, action, 0.0, -sigmaBankVal)

    # positive sigma
    posIndicesIntruder, posProbsIntruder = nextStatesSigmaAction(mdp, state, action, 0.0, sigmaBankVal)

    return [negIndicesOwnship, posIndicesOwnship, negIndicesIntruder, posIndicesIntruder], 
           [negProbsOwnship, posProbsOwnship, negProbsIntruder, posProbsIntruder]

end # function sigmaBank


function nextStatesSigma(mdp::SCA, state::State, action::Action)
    
    trueNextState = getNextState(state, action)
    gridNextState = getGridState(trueNextState)
    
    if trueNextState.clearOfConflict
        return [mdp.nStates], [1.0]
    else
        nextStateIndices, probs = interpolants(mdp.grid, gridNextState)
        return nextStateIndices, probs
    end # if

end # function nextStatesSigma


function nextStatesSigmaAction(
        mdp::SCA,
        state::State,
        action::Action,
        sigmaTurnOwnship::Float64,
        sigmaTurnIntruder::Float64)
    
    trueNextState = getNextState(state, action, sigmaTurnOwnship, sigmaTurnIntruder)
    gridNextState = getGridState(trueNextState)
    
    if trueNextState.clearOfConflict
        return [mdp.nStates], [1.0]
    else
        nextStateIndices, probs = interpolants(mdp.grid, gridNextState)
        return nextStateIndices, probs
    end # if

end # function nextStatesSigmaAction


function getGridState(state::State)
    
    return [
        state.x,
        state.y,
        state.bearing,
        state.speedOwnship,
        state.speedIntruder]

end # function getGridState


function index2state(mdp::SCA, stateIndices::Vector{Int64})
    
    states = Array(State, length(stateIndices))
    
    for index = 1:length(stateIndices)
        states[index] = gridState2state(ind2x(mdp.grid, index))
    end # for index
    
    return states
    
end # function index2state


function gridState2state(gridState::Vector{Float64})
    
    return State(
        gridState[1],  # x
        gridState[2],  # y
        gridState[3],  # bearing
        gridState[4],  # speedOwnship
        gridState[5],  # speedIntruder
        false)  # clearOfConflict

end # function gridState2state

end # module PilotSCAs