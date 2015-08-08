module PilotSCAIterators

export StateIterator, ActionIterator

using PilotSCAs


type StateIterator
    
    xs::Vector{Float64}
    ys::Vector{Float64}
    bearings::Vector{Float64}
    speedsOwnship::Vector{Float64}
    speedsIntruder::Vector{Float64}
    
    xCounter::Int64
    yCounter::Int64
    bearingCounter::Int64
    speedOwnshipCounter::Int64
    speedIntruderCounter::Int64

    nStates::Int64

    function StateIterator(
            xs::Vector{Float64},
            ys::Vector{Float64},
            bearings::Vector{Float64},
            speedsOwnship::Vector{Float64},
            speedsIntruder::Vector{Float64})
        
        nStates = 
            length(xs) * 
            length(ys) * 
            length(bearings) * 
            length(speedsOwnship) * 
            length(speedsIntruder) +
            1  # +1 for clear of conflict state

        return new(
            xs,  # xs
            ys,  # ys
            bearings,  # bearings
            speedsOwnship,  # speedsOwnship
            speedsIntruder,  # speedsIntruder
            1,  # xCounter
            1,  # yCounter
            1,  # bearingCounter
            1,  # speedOwnshipCounter
            1,  # speedIntruderCounter
            nStates)  # nStates
          
    end # function StateIterator
    
end # type StateIterator


Base.start(states::StateIterator) = 1
Base.done(states::StateIterator, stateCounter::Int64) = stateCounter - 1 == states.nStates


function Base.next(states::StateIterator, stateCounter::Int64)
    
    nextState = State(
        states.xs[states.xCounter],  # x
        states.ys[states.yCounter],  # y
        states.bearings[states.bearingCounter],  # bearing
        states.speedsOwnship[states.speedOwnshipCounter],  # speedOwnship
        states.speedsIntruder[states.speedIntruderCounter],  # speedIntruder
        false)  # clearOfConflict

    updateCounters!(states)

    if stateCounter - 1 == states.nStates - 1
        nextState.clearOfConflict = true
    end # if

    return nextState, stateCounter + 1

end # function Base.next


function updateCounters!(states::StateIterator)

    if (states.xCounter) % length(states.xs) == 0
        states.xCounter = 1
        if (states.yCounter) % length(states.ys) == 0
            states.yCounter = 1
            if (states.bearingCounter) % length(states.bearings) == 0
                states.bearingCounter = 1
                if (states.speedOwnshipCounter) % length(states.speedsOwnship) == 0
                    states.speedOwnshipCounter = 1
                    if (states.speedIntruderCounter) % length(states.speedsIntruder) == 0
                        states.speedIntruderCounter = 1  # we are done
                    else
                        states.speedIntruderCounter += 1
                    end # if: speedIntruder
                else
                    states.speedOwnshipCounter += 1
                end # if: speedOwnship
            else
                states.bearingCounter += 1
            end # if: bearing
        else 
            states.yCounter += 1
        end # if: y
    else
        states.xCounter += 1
    end # if: x

end # function updateCounters!


type ActionIterator
    
    actions::Vector{Symbol}
    
    ownshipCounter::Int64
    intruderCounter::Int64

    nActions::Int64

    function ActionIterator(actions::Vector{Symbol})

        return new(
            actions,  # actions
            1,  # ownshipCounter
            1,  # intruderCounter
            length(actions)^2)  # nActions

    end # function ActionIterator

end # type ActionIterator


Base.start(actions::ActionIterator) = 1
Base.done(actions::ActionIterator, actionCounter::Int64) = actions.nActions == actionCounter - 1


function Base.next(actions::ActionIterator, actionCounter::Int64)
    
    nextAction = Action(
        actions.actions[actions.ownshipCounter],  # actionOwnship
        actions.actions[actions.intruderCounter])  # actionIntruder

    updateCounters!(actions)
    
    return nextAction, actionCounter + 1

end # function Base.next


function updateCounters!(actions::ActionIterator)

    if (actions.ownshipCounter) % length(actions.actions) == 0
        actions.ownshipCounter = 1
        if (actions.intruderCounter) % length(actions.actions) == 0
            actions.intruderCounter = 1
        else
            actions.intruderCounter += 1
        end # if: intruder
    else
        actions.ownshipCounter += 1
    end # if: ownship

end # function updateCounters!

end # module PilotSCAIterators