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
    transmission(from::Int, 
        B::AbstractMatrix{T},
        Qbb::AbstractMatrix{T}, 
        q::Q, 
        ::Type{Val{:BQbb}}
    ) where {T}

Given a transmission condition `q`, calculate the transmission effect using the
`:BQbb` method. 

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
function transmission(from::Int, 
  B::AbstractMatrix{T},
  Qbb::AbstractMatrix{T}, 
  q::Q, 
  ::Type{Val{:BQbb}}
) where {T}

    @info "Using method :BQbb to calculate transmission effect."
    var_and, var_not, multiplier = get_varnums_and_multiplier(q)
    return transmission(from, B, Qbb, var_and, var_not, multiplier, Val{:BQbb})
end
function transmission(
  from::Int, 
  B::Matrix{T}, 
  Qbb::Matrix{T}, 
  var_and::Vector{Int}, 
  var_not::Vector{Int}, 
  ::Type{Val{:BQbb}}
) where {T}

    B_tilde = copy(B)
    Qbb_tilde = copy(Qbb)
    for v in var_and
        apply_and!(B_tilde, Qbb_tilde, from, v)
    end
    for v in var_not
        apply_not!(B_tilde, Qbb_tilde, from, v)
    end
    irfs = (I - B_tilde) \ Qbb_tilde[:, from]
    # var_and and var_not are empty if effect reduced to total effect
    isempty(var_and) && isempty(var_not) && return irfs
    irfs[1:maximum(vcat(var_and, var_not))] .= T(NaN)
    return irfs
end 
function transmission(from::Int, 
  B::Matrix{T}, 
  Qbb::Matrix{T}, 
  var_and::Vector{Vector{Int}}, 
  var_not::Vector{Vector{Int}}, 
  multiplier::Vector{Number}, 
  ::Type{Val{:BQbb}}
) where {T}

    effects = Vector{Vector{T}}(undef, length(var_and))
    for i in 1:lastindex(var_and)
        v_and = var_and[i]
        v_not = var_not[i]
        m = multiplier[i]
        effects[i] = m .* transmission(from, B, Qbb, v_and, v_not, Val{:BQbb})
    end
    return sum(effects)
end