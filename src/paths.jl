# TODO: test all the funcitonality in this file


function _add_neighbours!(stack, path, current_variable, current_time, A, B)
    n = size(A, 1)

    # contemporanous neighbours
    contemporanous_neighbours = findall(x->x!=0, A[current_variable,:])
    contemporanous_neighbours = filter(x->x!=current_variable, contemporanous_neighbours)
    contemporanous_neighbours_time = [(x, current_time) for x in contemporanous_neighbours]
    for cn in contemporanous_neighbours_time
        push!(stack, push!(copy(path), cn))
    end

    # lagged neighbours
    lagged_neighbours = findall(x->x!=0, B[current_variable,:])
    # B contrains all lagged matrices bind along columns. So if 
    # we had four variables and the 6th entry is non-zero, then this 
    # would correspond to the second variable in the system at lag 2.
    lagged_neighbours_time = [
        ((mod(x,n)==0) ? n : mod(x, n), 
        current_time - ceil(Int,x/n))
        for x in lagged_neighbours
    ]
    for ln in lagged_neighbours_time
        push!(stack, push!(copy(path), ln))
    end
end

function find_paths(A::AbstractMatrix, B::AbstractMatrix, from::Tuple{Int, Int}, to::Tuple{Int, Int})
    stack = []
    final_paths = []

    push!(stack, [to])
    while length(stack) > 0
        current_path = pop!(stack)
        current_variable = current_path[end][1]
        current_time = current_path[end][2]
        if current_variable == from[1] && current_time == from[2]
            push!(final_paths, current_path)
            continue
        end
        if current_time < from[2]
            # We are generally either staying in the same period or going 
            # backwards. If the path crossed the from time, then we can 
            # no longer reach the from node
            continue
        end
        add_neighbours!(stack, current_path, current_variable, current_time, A, B)
    end

    return reverse.(final_paths)
end


"""

Find all directed causal paths from one variable to another variable in a SVAR. 

The algorithm goes backwards from the 'to' node to the 'from' node.

## Arguments 

- `svar::SVAR`: an `SVAR` model from the `MacroEconometrics.jl` package.
- `from::Tuple{Int, Int}`: (from variable number, from time period)
- `to::Tuple{Int, Int}`: (to variable number, to time period)

## Returns

An array of paths. 

## Notes

This method is computationally very expensive and therefore not provided for
estimated models. Especially for Bayesian estimated models or Frequentist models
using Bootstrapping, this method is pretty much infeasible since it takes too
much time. Use one of the alternative methods in this package or use the Algebra
of Transmission Effects to derive a custom effect. 

"""
function find_paths(svar::SVAR{E}, from::Tuple{Int, Int}, to::Tuple{Int, Int}) where {E<:FixedEstimated}
    A = svar.A.value
    B = svar.B.value

    return find_paths(A, B, from, to)
end


function calculate_path_effect(A::AbstractMatrix, B::AbstractMatrix, path)
    n = size(A, 1)
    path_effect = 1.0
    previous_node = path[1]
    for node in path[2:end]
        lag = node[2] - previous_node[2]
        if lag == 0
            # Use contemporanous coefficients
            path_effect *= -A[node[1],previous_node[1]]
        else 
            # Use lagged matrix
            col_idx = (lag-1)*n + previous_node[1]
            path_effect *= B[node[1],col_idx]
        end
        previous_node = node
    end
    return path_effect
end
function calculate_path_effect(svar::SVAR{E}, path) where {E<:FixedEstimated}
    A = svar.A.value
    B = svar.B.value

    return calculate_path_effect(A, B, path)
end
"""
We can also do this for Bayesian estimated methods, since, given the path, the calculations are fast.
"""
function calculate_path_effect(svar::SVAR{E}, path) where {E<:BayesianEstimated}
    ndraws = size(svar.A.value, 3)
    nchains = size(svar.A.value, 4)
    path_effects = Array{Float64}(undef, ndraws, nchains)
    for d in 1:ndraws
        for c in 1:nchains
            A = svar.A.value[:,:,d,c]
            B = svar.B.value[:,:,d,c]
            path_effects[d, c] = calculate_path_effect(A, B, path)
        end
    end
    return path_effects
end

function mediation(svar::SVAR{E}, from::Tuple{Int, Int}, to::Tuple{Int, Int}, condition::Function) where {E<:FixedEstimated}
    paths = find_paths(svar, from, to)
    return mediation(svar, paths, condition)
end
"""
This should work for all estimation methods, since, given the paths, the calculations are fast.
"""
function mediation(svar::SVAR{E}, paths, condition::Function) where {E<:Estimated}
    mediating_paths = filter(condition, paths)
    mediating_effect = 0
    for p in mediating_paths
        mediating_effect .+= calculate_path_effect(svar, p)
    end
    return mediating_effect
end