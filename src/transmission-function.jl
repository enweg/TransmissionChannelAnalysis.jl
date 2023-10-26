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
    apply_and!(B::AbstractMatrix{T}, Qbb::AbstractMatrix{T}, from::Int, var::Int)
    
Manipulate `B` and `Qbb` so that `var` lies on all paths. This corresponds to
zeroing out all edges going directly from the shock to any variables ordered
after `var` and zeroing out any edges going from variables ordered before `var`
to any variables ordered after `var`.  

## Arguments

- `B::AbstractMatrix{T}`: Part of the structural transmission representation in
  $WEGNER. See also [`make_structural_B`](@ref). 
- `Qbb::AbstractMatrix{T}`: Part of the structural transmission representation
  in $WEGNER. See also [`make_structural_Qbb`](@ref). 
- `from::Int`: The shock number. 
- `var::Int`: The variable number that must lie on all paths. 

## Notes

- This function is meant for internal use only. 
"""
function apply_and!(B::AbstractMatrix{T}, Qbb::AbstractMatrix{T}, from::Int, var::Int) where {T}
    Qbb[(var+1):end, from] .= 0
    B[(var+1):end, 1:(var-1)] .= 0
    return B, Qbb
end

"""
    apply_not!(B::AbstractMatrix{T}, Qbb::AbstractMatrix{T}, from::Int, var::Int)
    
Manipulate `B` and `Qbb` so that `var` lies on no paths. This corresponds to
zeroing out the edge from the shock to `var`, and zeroing out all edges from
variables ordered before `var` to `var`. The paper mentions also zeroing out
edges leaving `var`, but this is not necessary.  

## Arguments

- `B::AbstractMatrix{T}`: Part of the structural transmission representation in
  $WEGNER. See also [`make_structural_B`](@ref). 
- `Qbb::AbstractMatrix{T}`: Part of the structural transmission representation
  in $WEGNER. See also [`make_structural_Qbb`](@ref). 
- `from::Int`: The shock number. 
- `var::Int`: The variable number that must lie on all paths. 

## Notes

- This function is meant for internal use only. 
"""
function apply_not!(B::AbstractMatrix{T}, Qbb::AbstractMatrix{T}, from::Int, var::Int) where {T}
    Qbb[var, from] = 0
    B[var, :] .= 0
    return B, Qbb
end

"""
    transmission(from::Int, B::Matrix{T}, Qbb::Matrix{T}, q::Q) where {T}

Given a transmission condition `q`, calculate the transmission effect. 

## Arguments

- `from::Int`: Shock number. 
- `B::AbstractMatrix{T}`: Part of the structural transmission representation in
  $WEGNER. See also [`make_structural_B`](@ref). 
- `Qbb::AbstractMatrix{T}`: Part of the structural transmission representation
  in $WEGNER. See also [`make_structural_Qbb`](@ref). 
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

B = randn(k*(h+1), k*(h+1))
Qbb = randn(k*(h+1), k*(h+1))

effect = transmission(1, B, Qbb, cond)
```
"""
function transmission(from::Int, B::AbstractMatrix{T}, Qbb::AbstractMatrix{T}, q::Q) where {T}
    var_and, var_not, multiplier = get_varnums_and_multiplier(q)
    return transmission(from, B, Qbb, var_and, var_not, multiplier)
end
function transmission(from::Int, B::Matrix{T}, Qbb::Matrix{T}, var_and::Vector{Int}, var_not::Vector{Int}) where {T}
    B_tilde = copy(B)
    Qbb_tilde = copy(Qbb)
    for v in var_and
        apply_and!(B_tilde, Qbb_tilde, from, v)
    end
    for v in var_not
        apply_not!(B_tilde, Qbb_tilde, from, v)
    end
    irfs = (I - B_tilde) \ Qbb_tilde[:, from]
    irfs[1:maximum(vcat(var_and, var_not))] .= T(NaN)
    return irfs
end 
function transmission(from::Int, B::Matrix{T}, Qbb::Matrix{T}, var_and::Vector{Vector{Int}}, var_not::Vector{Vector{Int}}, multiplier::Vector{Number}) where {T}
    effects = Vector{Vector{T}}(undef, length(var_and))
    for i in 1:lastindex(var_and)
        v_and = var_and[i]
        v_not = var_not[i]
        m = multiplier[i]
        effects[i] = m .* transmission(from, B, Qbb, v_and, v_not)
    end
    return sum(effects)
end