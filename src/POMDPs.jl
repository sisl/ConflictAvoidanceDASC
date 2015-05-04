#=
Description
    Module for POMDP type objects with state, action, and
    observation space, etc. defined.
=#

module POMDPs

export
    POMDP


type POMDP

    states          :: Matrix
    actions         :: Array
    observations    :: Matrix
    get_next_state  :: Function
    get_reward      :: Function
    get_obs_prob    :: Function
    term_state      :: Vector{Float64}
    
    function POMDP(states::Matrix, actions::Array, observations::Matrix, 
                   get_next_state::Function, get_reward::Function, 
                   get_obs_prob::Function, term_state::Vector{Float64})
        return new(states, actions, observations, get_next_state, get_reward, 
                   get_obs_prob, term_state)
    end # function POMDP
    
end # type POMDP

end # module POMDPs