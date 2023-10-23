"""
    create_transmission_function(from::Int, [to::Int,] condition::SymbolicUtils.BasicSymbolic{Bool})

Create a transmission function based on a boolean condition.

This function generates a function that calculates the transmission effect from
the `from`th variables to the optional `to`th variable. If no `to` is given, the
effect is calculated for all destination nodes. The transmission paths must
satisfy the Boolean statement given in `condition` which can be created using
[`make_condition`](@ref).

## Arguments

- `from::Int`: Index of the variable from which the transmission starts.
- `[to::Int]`: (Optional) Index of the variable to which the transmission ends. 
  If omitted, the returned function calcuates the effect for all destination
  nodes. 
- `condition::SymbolicUtils.BasicSymbolic{Bool}`: Boolean condition specifying 
  the transmission mechanism. The condition is represented as a boolean expression 
  using variables starting with `x` and can be created using [`make_condition`](@ref).

## Returns

- Returns a transmission function that can be applied to impulse response functions 
  (`irfs`) and orthogonalized impulse response functions (`irfs_ortho`) to compute 
  the transmission effects. The first argument in the returned function is
  `irfs` and the second arguement is `irfs_ortho`. Both should be matrices
  corresponding the the structural and orthogonalised IRFs respectively. Both
  should be square and of dimension `k*(h+1)` with `h` being the maximum horizon
  and `k` being the number of variables. The structural version can, for
  example, be obtained via `irfs = inv(I - B) * Qbb`, with `B` and `Qbb`
  obtained from [`to_structural_transmission_model`](@ref).

## Examples

```julia
s = "x2 & !x3"
cond = make_condition(s)
tf = create_transmission_function(1, cond)
irfs = randn(3, 3)
irfs_ortho = randn(3, 3)
tf(irfs, irfs_ortho)
````

"""
function create_transmission_function(from::Int, to::Int, condition::SymbolicUtils.BasicSymbolic{Bool})
    terms, multiplier, _, variable_nums = helper_Q(condition)
    println(variable_nums)
    any(contains_nots.(terms)) && error("Something went wrong in the simplification process. The terms still include nots.\n $terms")
    transmission_function(irfs, irfs_ortho) = begin
        terms = Vector{eltype(irfs)}(undef, length(variable_nums))
        for (i, term) in enumerate(variable_nums)
            if isnothing(term[1]) && length(term) == 1
                terms[i] = irfs[to, from]
                continue
            end

            terms[i] = zero(eltype(irfs))
            terms[i] += irfs[term[1], from]
            for i = 1:(lastindex(term)-1)
                terms[i] *= irfs_ortho[term[i+1], term[i]] / irfs_ortho[term[i], term[i]]
            end
            terms[i] *= irfs_ortho[to, term[end]] / irfs_ortho[term[end], term[end]]
        end
        terms' * multiplier, terms
    end
    
    return transmission_function
end
function create_transmission_function(from::Int, condition::SymbolicUtils.BasicSymbolic{Bool})
    terms, multiplier, _, variable_nums = helper_Q(condition)
    any(contains_nots.(terms)) && error("Something went wrong in the simplification process. The terms still include nots.\n $terms")
    transmission_function(irfs, irfs_ortho) = begin
        terms = zeros(eltype(irfs), size(irfs, 1), length(variable_nums))
        for (i, term) in enumerate(variable_nums)
            if isnothing(term[1]) && length(term) == 1
                terms[:, i] .= irfs[:, from]
                continue
            end

            # terms[:, i] .= zero(eltype(irfs))
            terms[:, i] .+= irfs[term[1], from]
            for j = 1:(lastindex(term)-1)
                terms[:, i] .*= irfs_ortho[term[j+1], term[j]] ./ irfs_ortho[term[j], term[j]]
            end
            terms[:, i] .*= irfs_ortho[:, term[end]] ./ irfs_ortho[term[end], term[end]]
        end
        terms * multiplier
    end
    
    return transmission_function
end