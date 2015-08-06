module SCAIterators

export StateIterator, ActionIterator

using SCAs


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


function Base.done(states::StateIterator, stateCounter::Int64)
    return stateCounter - 1 == states.nStates
end # function Base.done


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

    if (states.speedIntruderCounter) % length(states.speedsIntruder) == 0
        states.speedIntruderCounter = 1
        if (states.speedOwnshipCounter) % length(states.speedsOwnship) == 0
            states.speedOwnshipCounter = 1
            if (states.bearingCounter) % length(states.bearings) == 0
                states.bearingCounter = 1
                if (states.yCounter) % length(states.ys) == 0
                    states.yCounter = 1
                    if (states.xCounter) % length(states.xs) == 0
                        states.xCounter = 1  # we are done
                    else
                        states.xCounter += 1
                    end # if: x
                else
                    states.yCounter += 1
                end # if: y
            else
                states.bearingCounter += 1
            end # if: bearing
        else 
            states.speedOwnshipCounter += 1
        end # if: speedOwnship
    else
        states.speedIntruderCounter += 1
    end # if: speedIntruder

end # function updateCounters!


type ActionIterator
    actions::Vector{Symbol}
end # type ActionIterator


Base.start(actions::ActionIterator) = 1
Base.done(actions::ActionIterator, iaction::Int64) = length(actions.actions) == iaction - 1
Base.next(actions::ActionIterator, iaction::Int64) = actions.actions[iaction], iaction + 1

end # module