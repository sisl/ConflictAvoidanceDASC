#=
Description
    Module for POMDP type objects with state, action, and
    observation space, etc. defined.
=#

module POMDPs

export POMDP


type POMDP

    states          :: Matrix
    actions         :: Array
    observations    :: Matrix
    term_state      :: Vector{Float64}
    
    function POMDP(states::Matrix, actions::Array, observations::Matrix, 
                   term_state::Vector{Float64})
        return new(states, actions, observations, term_state)
    end # function POMDP
    
end # type POMDP

end # module POMDPs