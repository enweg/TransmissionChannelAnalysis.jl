"""
    make_condition(string::String)

Given a `string`, make a `SymbolicUtils.jl` Boolean expression. 

Transmission mechanisms are described using Boolean statements involving the 
variables in the system. Each variable starts with `x` followed by a number.
For example, given a three variable VAR(1), y_{1,t} -> x_1, y_{2, t} -> x_2, 
y_{3, t} -> x_3, y_{1, t+1} -> x_4, y_{2, t+1} -> x_5, ... Boolean statements 
then involve expressions in the `x` variables and define which paths can be
taken. Each path involved in the transmission mechanism must satisfy the Boolean
statement. 

## Arguments

- `string::String`: A Boolean statement given as a string. Variables must start
  with `x` for them to be valid variables. 

## Returns 

- Returns a `SymbolicUtils.jl` Boolean expression that can be used in
  [`create_transmission_function`](@ref). 

## Examples

```julia
s = "x2 & !x3"
cond = make_condition(s)
tf = create_transmission_function(1, cond)
irfs = randn(3, 3)
irfs_ortho = randn(3, 3)
tf(irfs, irfs_ortho)
```

## Notes

- Ensure that variables in the input string are correctly formatted as `x`
  followed by a number. 
- The resulting Boolean expression can be used in [`create_transmission_function`](@ref)
  to create a function that calculates the transmission effect. 

"""
function make_condition(string::String)
    variables = [Symbol(m.match) for m in eachmatch(r"x\d+", string)]
    for v in variables
        eval(Meta.parse("$v = SymbolicUtils.Sym{Bool}(:$v)"))
        # eval(:($v = SymbolicUtils.Sym{Bool}((:)($v))))
    end
    return eval(Meta.parse(string))
end