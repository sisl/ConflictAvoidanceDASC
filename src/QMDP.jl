#=
Description
    Module for QMDP solver to generate approximate POMDP solution for
    the DoubleUAV problem.
=#

module QMDP

using POMDPs, DoubleUAVs, GridInterpolations

export qmdp


const MAX_ITER = 1000
const GAMMA = 0.95
const ALPHA_TOL = 1e-6
const STATE_DIM = 5
const TERM_STATE_VAR = 1e5
const TERM_STATE = TERM_STATE_VAR * ones(STATE_DIM)


function get_max_alphas!(V::Vector{Float64}, alpha::Matrix{Float64}, 
                         nStates::Int64)
    nActions = size(alpha, 2)
    for istate = 1:nStates
        V[istate] = alpha[istate, 1]
        for iaction = 2:nActions
            if V[istate] < alpha[istate, iaction]
                V[istate] = alpha[istate, iaction]
            end # if
        end # for iaction
    end # for istate
end # function get_max_alphas!


function qmdp(states::Matrix, actions::Matrix, sigmas::Matrix{Float64}, 
              weights::Vector{Float64}, verbose::Bool, maxIter::Int64, 
              gamma::Float64, alphaTol::Float64)
    nStates = size(states, 2)
    nActions = size(actions, 2)
    nSigmas = size(sigmas, 2)
    g = RectangleGrid(sort(unique(vec(states[1, 1:end - 1]))),
                      sort(unique(vec(states[2, 1:end - 1]))),
                      sort(unique(vec(states[3, 1:end - 1]))),
                      sort(unique(vec(states[4, 1:end - 1]))),
                      sort(unique(vec(states[5, 1:end - 1]))))
    alpha = zeros(nStates, nActions)
    V = zeros(nStates)
    
    @printf("Running QMDP alpha vectors approximation...\n")
    cputime = 0
    for iter = 1:maxIter
        tic()
        residual = 0
        for istate = 1:nStates
            state = TERM_STATE
            if istate != nStates # not terminal state
                state = ind2x(g, istate)
            end # if

            for iaction = 1:nActions
                action = actions[:, iaction]
                vnext = 0

                for isigma = 1:nSigmas
                    snext = get_next_state(state, action, sigmas, isigma)
                    if snext[1] == TERM_STATE_VAR
                        vnext = vnext + V[end] * weights[isigma]
                    else
                        vnext = vnext + interpolate(g, V, snext) * weights[isigma]
                    end # if
                end # for isigma

                alphaPrev = alpha[istate, iaction]
                alpha[istate, iaction] = get_reward(state, action) + 
                                         gamma * vnext
                residual = residual + (alpha[istate, iaction] - alphaPrev)^2
            end # for iaction
        end # for istate
        
        get_max_alphas!(V, alpha, nStates)

        iterTime = toq()
        cputime = cputime + iterTime
        if verbose
            @printf("Iteration %d: residual = %.2e, cputime = %.2e\n",
                    iter, residual, iterTime)
        end # if
        if residual < alphaTol; break; end
        if iter == maxIter
            println("Warning: maximum number of iterations reached;", 
                    "solution may be inaccurate")
        end # if
    end # for iter
    @printf("QMDP computations done!\ncputime = %.2e sec\n\n", cputime)
    return alpha
end # function qmdp


function qmdp(d::DoubleUAV, verbose::Bool = false, maxIter::Int64 = MAX_ITER, 
              gamma::Float64 = GAMMA, alphaTol::Float64 = ALPHA_TOL)
    p = d.pomdp
    alpha = qmdp(p.states, p.actions, d.sigmas, d.weights, verbose, maxIter, 
                 gamma, alphaTol)
    return alpha
end # function qmdp

end # module QMDP