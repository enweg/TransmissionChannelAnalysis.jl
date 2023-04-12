
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

"""
    find_paths(svar::SVAR{<:FixedEstimated}, from::Tuple{Int, Int}, to::Tuple{Int, Int})
    find_paths(A::AbstractMatrix, B::AbstractMatrix, from::Tuple{Int, Int}, to::Tuple{Int, Int})

Find all directed causal paths from one variable to another variable in a SVAR. 

The algorithm goes backwards from the 'to' node to the 'from' node.

## Arguments 

- `A::AbstractMatrix`: Contemporanous matrix of a SVAR of dimensions n×n where n
  is the number of variables in the SVAR
- `B::AbstractMatrix`: Lag matrix of the SVAR. This is of the form [B_1 B_2 B_3
  ... B_p] so that B is of dimension n×np where p is the number of lags in the
  SVAR
- `svar::SVAR`: an `SVAR` model from the `MacroEconometrics.jl` package.
- `from::Tuple{Int, Int}`: (from variable number, from time period)
- `to::Tuple{Int, Int}`: (to variable number, to time period)

## Returns

An array of paths where each path is an array of Tuple{Int, Int} where the
tuples are of the form (Variable, Time period).

## Notes

This method is computationally very expensive and therefore not provided for
estimated models. Especially for Bayesian estimated models or Frequentist models
using Bootstrapping, this method is pretty much infeasible since it takes too
much time. Use one of the alternative methods in this package or use the Algebra
of Transmission Effects to derive a custom effect.

"""
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
        _add_neighbours!(stack, current_path, current_variable, current_time, A, B)
    end

    return reverse.(final_paths)
end
function find_paths(svar::SVAR{E}, from::Tuple{Int, Int}, to::Tuple{Int, Int}) where {E<:FixedEstimated}
    A = svar.A.value
    B = svar.B.value

    return find_paths(A, B, from, to)
end

"""
    calculate_path_effect(A::AbstractMatrix, B::AbstractMatrix, path)
    calculate_path_effect(svar::SVAR{<:FixedEstimated}, path)

Calculate the causal effect along a directed path. 

## Arguments

- `A::AbstractMatrix`: Contemporanous matrix of the SVAR of dimension n×n where
  n is the number of variables in the SVAR. 
- `B::AbstractMatrix`: Lag matrix of the SVAR: B=[B_1 ... B_p] and thus B is of
  dimension n×np where p is the number of lags in the SVAR. 
- `path::Vector{Tuple{Int, Int}}`: The directed path. Each tuple is of the form
  (variable number, time period)

## Returns 

- Returns a `Real`
"""
function calculate_path_effect(A::AbstractMatrix, B::AbstractMatrix, path::Vector{Tuple{Int, Int}})
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
function calculate_path_effect(svar::SVAR{E}, path::Vector{Tuple{Int, Int}}) where {E<:FixedEstimated}
    A = svar.A.value
    B = svar.B.value

    return calculate_path_effect(A, B, path)
end
"""
    calculate_path_effect(svar::SVAR{<:BayesianEstimated}, path)

## Returns

- In case of a Bayesian estimated SVAR, the returned object will be of type
`BayesianEstimated` with the last dimension corresponding to the chains, and the
second to last dimension corresponding to the draws.

"""
function calculate_path_effect(svar::SVAR{E}, path) where {E<:BayesianEstimated}
    f = (A, B) -> calculate_path_effect(A, B, path)
    path_effects = map(f, svar.A, svar.B)
    return  path_effects
end


"""
    mediation(svar::SVAR{<:FixedEstimated}, from::Tuple{Int, Int}, to::Tuple{Int, Int}, condition::Function)
    mediation(svar::SVAR{<:Estimated}, paths, condition::Function) where {E<:Estimated}

Calculate the mediation effect. This is the total path effect along multiple
paths. The first version includes a seach for paths which is computationally
expensive and thus only implemented for fixed SVARs and not for estimated SVARs.  

## Arguments

- `svar::SVAR`: SVAR model coming from `MacroEconometrics.jl`. 
- `from::Tuple{Int, Int}`: Origin node of type (variable number, time period). 
- `to::Tuple{Int, Int}`: Destination node. See `from`. 
- `condition::Function`: A function returning a boolean. `condition` will be
  forwarded to `filter` in order to filter the paths to only those that one is
  interested in.
- `paths::Vector{Vector{Tuple{Int, Int}}}`: A vector of paths. Also see
  [`calculate_path_effect`](@ref)

## Returns

- The total effect along the filtered paths. 
- In case of a Bayesian estimated SVAR, this will be a `BayesianEstimated`
  object with the last dimension being the chains, and the second to last
  dimension corresponding to the draws.
"""
function mediation(svar::SVAR{E}, from::Tuple{Int, Int}, to::Tuple{Int, Int}, condition::Function) where {E<:FixedEstimated}
    paths = find_paths(svar, from, to)
    return mediation(svar, paths, condition)
end
function mediation(svar::SVAR{E}, paths::Vector{Vector{Tuple{Int, Int}}}, condition::Function) where {E<:Estimated}
    mediating_paths = filter(condition, paths)
    if length(mediating_paths) == 0
        throw(ErrorException("All paths have been sorted out by the condition"))
    end
    mediating_effect = 0 .* calculate_path_effect(svar, mediating_paths[1])
    for p in mediating_paths
        mediating_effect += calculate_path_effect(svar, p)
    end
    if E <: BayesianEstimated
        return BayesianEstimated(mediating_effect, nothing)
    end
    return mediating_effect
end