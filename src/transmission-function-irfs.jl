using Combinatorics: combinations

"""
    calculate_Q_and_only(from::Int, 
        irfs::AbstractMatrix{T}, 
        irfs_ortho::AbstractMatrix{T}, 
        vars::AbstractVector{Int}, 
        multiplier::Number
    ) where {T}

Calculate the transmission effect of a transmission condition/query that
involves only ANDs. 

## Arguments

- `from::Int`: ID of shock. 
- `irfs::AbstractMatrix{T}`: IRFs in transmission form. See
  [`to_transmission_irfs`](@ref).
- `irfs_ortho::AbstractMatrix{T}`: Orthogonalised IRFs in transmission form. See
  [`to_transmission_irfs`](@ref).
- `multiplier::Number`: Multiplier.

## Returns 

- Returns a `Vector{T}` with entry `i` corresponding to the transmission effect
  on variable ``x_i``. If ``x_k`` is the variable in the transmission condition
  with the highest subscript, then all entries in the returned vector with index
  less thatn `k` are `NaN` since interpretation of those results is nonsensical.

## Notes

- Only for internal use. 

"""
function calculate_Q_and_only(from::Int, 
    irfs::AbstractMatrix{T}, 
    irfs_ortho::AbstractMatrix{T}, 
    vars::AbstractVector{Int}, 
    multiplier::Number
) where {T}

    if isempty(vars)
        # indicating TRUE
        return irfs[:, from] .* multiplier
    end
    
    vars = sort(vars)
    
    effect = zeros(T, size(irfs, 1))
    effect .= multiplier * irfs[vars[1], from]
    for i in 1:(lastindex(vars)-1)
        effect .*= irfs_ortho[vars[i+1], vars[i]] ./ irfs_ortho[vars[i], vars[i]]
    end
    effect .*= irfs_ortho[:, vars[end]] ./ irfs_ortho[vars[end], vars[end]]

    # effect[1:maximum(vars)] .= NaN
    return effect
end

"""
    transmission(from::Int, 
        irfs::AbstractMatrix{T}, 
        irfs_ortho::AbstractMatrix{T}, 
        q::Q, 
        ::Type{Val{:irfs}}
    ) where {T}

Given a transmission condition `q`, calculate the transmission effect using the
`:irfs` method. 

## Arguments

- `from::Int`: Shock number. 
- `irfs::AbstractMatrix{T}`: Impulse responses. These should be in the form of
  the structural transmission model. See also [`to_transmission_irfs`](@ref). 
- `irfs_ortho::AbstractMatrix{T}`: Orthogonalised IRFs. These should be in the
  form of the structural transmission model. See also [`to_transmission_irfs`](@ref). 
- `q::Q`: A transmission condition. See also [`Q`](@ref). 

## Returns 

- Returns a `Vector{T}` with entry `i` corresponding to the transmission effect
  on variable ``x_i``. If ``x_k`` is the variable in the transmission condition
  with the highest subscript, then all entries in the returned vector with index
  less thatn `k` are `NaN` since interpretation of those results is nonsensical. 

## Examples

```julia
k = 6
h = 3
s = "(x1 | x2) & !x3"
cond = make_condition(s)

irfs = randn(k, k, h+1)
irfs_ortho = randn(k, k, h+1)

irfs = to_transmission_irfs(irfs)
irfs_ortho = to_transmission_irfs(irfs_ortho)

effect = transmission(1, irfs, irfs_ortho, cond; method = :irfs)
```
"""
function transmission(from::Int, 
    irfs::AbstractMatrix{T}, 
    irfs_ortho::AbstractMatrix{T}, 
    q::Q, 
    ::Type{Val{:irfs}}
) where {T}

    @info "Using method :irfs to calculate transmission effect."

    vars_and, vars_not, multiplier = get_varnums_and_multiplier(q)
    effects = Vector{Vector{T}}(undef, length(vars_and))
    Threads.@threads for i = 1:lastindex(vars_and)
        v_and = vars_and[i]
        v_not = vars_not[i]
        m = multiplier[i]
        effects[i] = transmission(from, irfs, irfs_ortho, v_and, v_not, m, Val{:irfs})
    end
    return sum(effects)
end
function transmission(from::Int, 
    irfs::AbstractMatrix{T}, 
    irfs_ortho::AbstractMatrix{T},
    vars_and::AbstractVector{Int}, 
    vars_not::AbstractVector{Int}, 
    multiplier::Number, 
    ::Type{Val{:irfs}}
) where {T}

    effect = calculate_Q_and_only(from, irfs, irfs_ortho, vars_and, one(T))
    combs = combinations(vars_not)
    for c in combs
        v = unique(vcat(vars_and, c))
        effect += (-1)^length(c) * calculate_Q_and_only(from, irfs, irfs_ortho, v, one(T))
    end
    return multiplier * effect
end
