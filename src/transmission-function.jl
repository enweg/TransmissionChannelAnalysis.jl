using LinearAlgebra

"""
    get_varnums_and_multiplier(q::Q)

Obtain the AND and NOT expressions of a transmission condition. Also return the
multiplier for each term. 

A valid transmission condition is a set of terms involving only AND and NOT
expressions. For each term, the AND and NOT expressions are collected and a
vector of the respective variable numbers is returned. 

## Arguments

- `q::Q`: A transmission condition. See also [`Q`](@ref). 

## Returns

1. `and_nums::Vector{Vector{Int}}`: Contains for each term a vector of variable
   numbers that are included via AND in the term. 
2. `and_not_nums::Vector{Vector{Int}}`: Contains for each term a vector of
   variable numbers that are included via NOT in the term. 
3. `multiplier::Vector{Number}`: Contains for each term the multiplier. 

## Examples

```julia
q = Q(["x1", "!x2", "x1 & !x2"], [1, 2, 3])
and_nums, and_not_nums, multiplier = get_varnums_and_multiplier(q)
# output: 
# and_nums = [[1], [], [1]]
# and_not_nums = [[], [2], [2]]
# multiplier = [1, 2, 3]
```
"""
function get_varnums_and_multiplier(q::Q)
    and_nums = [collect([parse(Int, m[1]) for m in eachmatch(r"(?<!!)x(\d+)", q.vars[i])]) for i = 1:lastindex(q.vars)]
    and_not_nums = [collect([parse(Int, m[1]) for m in eachmatch(r"!x(\d+)", q.vars[i])]) for i in 1:lastindex(q.vars)]
    return and_nums, and_not_nums, q.multiplier
end

"""
    to_transmission_irfs(irfs::AbstractArray{T, 3})

Transform a standard three dimensional IRF array into a IRF matrix. The 

## Arguments

- `irfs::AbstractArray{T, 3}`: IRF Array of dimension n_variabels × n_shocks ×
  n_horizons with the first horizons corresponding to horizon 0. 

## Returns

- `Matrix{T}` of dimension (n_variables * n_horizons) × (n_variables *
  n_horizons). This is the same as what would be obtained via
  ``(I-B)^{-1}\\mathbb{Q}`` using the notation of $WEGNER. 
"""
function to_transmission_irfs(irfs::AbstractArray{T, 3}) where {T}
    k = size(irfs, 1)
    max_horizon = size(irfs, 3) - 1
    irfs = reduce(vcat, eachslice(irfs; dims = 3))
    irfs = reduce(hcat, [vcat(zeros(k*h, k), irfs[1:(end-k*h), :]) for h = 0:max_horizon])
    return irfs
end

"""
  transmission(from::Int, arr1::AbstractMatrix{T}, arr2::AbstractMatrix{T}, q::Q; method = :BQbb) where {T}


Given a transmission condition `q`, calculate the transmission effect using the
either the `:BQbb` method (the default), or the `:irfs` method. 

## Arguments

- `from::Int`: Shock number. 
- `arr1::AbstractMatrix{T}`. In case of `:BQbb` this must be `B`, in case of
  `:irfs` this must be `irfs`. See the documentation for the specific methods
  for `transmission(..., ::Type{Val{:BQbb}})` and `transmission(...,::Type{Val{:irfs}})`. 
- `arr2::AbstractMatrix{T}`: In case of `:BQbb` this must be `Qbb`, in case of
  `:irfs` this must be `irfs_ortho`. See the documentation for the specific methods
  for `transmission(..., ::Type{Val{:BQbb}})` and `transmission(...,::Type{Val{:irfs}})`.  
- `q::Q`: A transmission condition. See also [`Q`](@ref).

## Keyword Arguments

- `method::Symbol`: Either `:BQbb` in which case the transmission effect will be
  calculated using the second method in $WEGNER, or `:irfs` in which case the
  transmission effect is calculated using the first method in $WEGNER. 

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

B = randn(k*(h+1), k*(h+1))
Qbb = randn(k*(h+1), k*(h+1))

effect = transmission(1, B, Qbb, cond; method = :BQbb)
effect = transmission(1, B, Qbb, cond)  # same as above; default is :BQbb

irfs = randn(k, k, h+1)
irfs_ortho = randn(k, k, h+1)

irfs = to_transmission_irfs(irfs)
irfs_ortho = to_transmission_irfs(irfs_ortho)

effect = transmission(1, irfs, irfs_ortho, cond; method = :irfs)
```
"""
function transmission(from::Int, arr1::AbstractMatrix{T}, arr2::AbstractMatrix{T}, q::Q; method = :BQbb) where {T}
  return transmission(from, arr1, arr2, q, Val{method})
end