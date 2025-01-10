"""
    apply_and!(B::AbstractMatrix, Omega::AbstractMatrix, from::Int, var::Int)
    
Manipulate `B` and `Omega` so that `var` lies on all paths. This corresponds to
zeroing out all edges going directly from the shock to any variables ordered
after `var` and zeroing out any edges going from variables ordered before `var`
to any variables ordered after `var`.  

## Arguments

- `B::AbstractMatrix`: Part of the systems representation in $WEGNER.
- `Omega::AbstractMatrix`: Part of the systems representation in $WEGNER.
- `from::Int`: The shock number. 
- `var::Int`: The variable number that must lie on all paths. 

## Notes

- This function is meant for internal use only. 
- See also `apply_not!`, `make_B`, `make_Omega`, and `make_systems_form`
"""
function apply_and!(B::AbstractMatrix, Omega::AbstractMatrix, from::Int, var::Int)
    Omega[(var+1):end, from] .= 0
    B[(var+1):end, 1:(var-1)] .= 0
    return B, Omega
end

"""
    apply_not!(B::AbstractMatrix, Omega::AbstractMatrix, from::Int, var::Int)
    
Manipulate `B` and `Omega` so that `var` lies on no paths. This corresponds to
zeroing out the edge from the shock to `var`, and zeroing out all edges from
variables ordered before `var` to `var`. 

## Arguments

- `B::AbstractMatrix`: Part of the systems representation in $WEGNER.
- `Omega::AbstractMatrix`: Part of the systems representation in $WEGNER.
- `from::Int`: The shock number. 
- `var::Int`: The variable number that cannot lie on any paths. 

## Notes

- This function is meant for internal use only. 
- See also `apply_and!`, `make_B`, `make_Omega`, and `make_systems_form`
"""
function apply_not!(B::AbstractMatrix, Omega::AbstractMatrix, from::Int, var::Int)
    Omega[var, from] = 0
    B[var, 1:(var-1)] .= 0
    return B, Omega
end

"""
    transmission(
        from::Int, 
        B::AbstractMatrix{T},
        Omega::AbstractMatrix{T}, 
        q::Q, 
        ::Type{Val{:BOmega}}
    ) where {T}

Given a transmission condition `q`, calculate the transmission effect using the
`:BOmega` method. 

## Arguments

- `from::Int`: Shock number. 
- `B::AbstractMatrix{T}`: Part of the systems representation in $WEGNER. 
- `Omega::AbstractMatrix{T}`: Part of the systems representation in $WEGNER. 
- `q::Q`: A transmission condition. See also `Q`. 

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
Omega = randn(k*(h+1), k*(h+1))

effect = transmission(1, B, Omega, cond)
```
"""
function transmission(
    from::Int, 
    B::AbstractMatrix,
    Omega::AbstractMatrix, 
    q::Q, 
    ::Type{Val{:BOmega}}
)

    @info "Using method :BOmega to calculate transmission effect."
    var_and, var_not, multiplier = get_varnums_and_multiplier(q)
    return transmission(from, B, Omega, var_and, var_not, multiplier, Val{:BOmega})
end
function transmission(
    from::Int, 
    B::Matrix, 
    Omega::Matrix, 
    var_and::Vector{Int}, 
    var_not::Vector{Int}, 
    ::Type{Val{:BOmega}}
)

    B_tilde = copy(B)
    Omega_tilde = copy(Omega)
    for v in var_and
        apply_and!(B_tilde, Omega_tilde, from, v)
    end
    for v in var_not
        apply_not!(B_tilde, Omega_tilde, from, v)
    end
    irfs = (I - B_tilde) \ Omega_tilde[:, from]
    # var_and and var_not are empty if effect reduced to total effect
    isempty(var_and) && isempty(var_not) && return irfs

    # In the theory we assume that all var_and and var_not indices are less 
    # than the index of the outcome variable (otherwise the qustion does not make
    # sense to begin with). We still assumes this here. However, a user might 
    # not be aware of this requirement and may thus interpret the outcome of a 
    # variable that has an index less than the indices in var_and and var_not. 
    # Having an index less than var_not is not a problem, since the path did 
    # not go through any variables in var_not. However, having an index less than 
    # var_and is a problem, since the path cannot have gone through at least 
    # one variable in var_and. Thus, correct is that all effects less than 
    # the maximum(var_and) are zero. 
    irfs[1:maximum(var_and)] .= 0

    return irfs
end 
function transmission(
    from::Int, 
    B::Matrix, 
    Omega::Matrix, 
    var_and::Vector{Vector{Int}}, 
    var_not::Vector{Vector{Int}}, 
    multiplier::Vector{Number}, 
    ::Type{Val{:BOmega}}
)

    effects = Vector{Vector{Float64}}(undef, length(var_and))
    for i in 1:lastindex(var_and)
        v_and = var_and[i]
        v_not = var_not[i]
        m = multiplier[i]
        effects[i] = m .* transmission(from, B, Omega, v_and, v_not, Val{:BOmega})
    end
    return sum(effects)
end
