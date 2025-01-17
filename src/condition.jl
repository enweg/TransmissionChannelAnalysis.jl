"""
    make_condition(s::String)
    make_condition(s_y::String, order::AbstractVector{<:Int})

Make a transmission condition, i.e. Q(b), out of a string. 

Transmission channels are described using Boolean statements involving the 
variables in the dynamic model. `make_condition` allows for specifying these
Boolean conditions as a string which is then converted to an internal 
representation allowing the computation of transmission channels. 

Two ways of specifying the Boolean conditions exist: 

1. `make_condition(s::String)` takes the Boolean condition in the systems form 
    of $WEGNER, i.e. variables must start with `x` followed by a number. 
    For example, given a three variable VAR(1), `y_{1,t} -> x_1`, `y_{2, t} -> x_2`, 
    `y_{3, t} -> x_3`, `y_{1, t+1} -> x_4`, `y_{2, t+1} -> x_5`, ... Boolean statements 
    then involve expressions in the `x` variables and define which paths can be
    taken. Each path involved in the transmission mechanism must satisfy the Boolean
    statement. 
2. `make_condition(s_y::String, order::AbstractVector{<:Int})` does the same 
    as the first method, however the Boolean condition can be specified using 
    the variables of the dynamic systems, i.e. `y`. Variables must then be
    specified using `y_{i,t}` where `i` is the variable number and `t` is 
    the period. At all times `t >= 0` with `0` denoting the contemporaneous 
    horizon.

## Arguments

- `s::String`: A Boolean statement given as a string. Variables must start
  with `x` for them to be valid variables. 
- `s_y::String`: A Boolean statement given as a string. Variabls must have the 
  form `y_{i,t}` where `i` is the variable number and `t >= 0` is the time period. 
  `t=0` corresponds to the contemporaneous horizon.
- `order::AbstractVector{<:Int}`: The variable ordering defined by the transmission
  matrix. 

## Returns 

- Returns a transmission condition. See also `Q`.

## Examples

```julia
s = "x2 & !x3"
cond = make_condition(s)
```

```julia
s_y = "y_{1,0} & !y_{1,1}"
order = [3,1,2]
cond = make_condition(s_y, order)
```

## Notes

- Boolean conditions can consist of AND (&), NOT (!), OR (|), and parentheses. 
- The resulting transmission condition can be used in `transmission` to
  calculate the transmission effect.

"""
function make_condition(s::String)
  vars = collect([m.match for m in eachmatch(r"(x\d+)", s)])
  for v in unique(vars)
      sv = Symbol(string(v))
      eval(:($sv = Q(string($v))))
  end
  return eval(Meta.parse(s))
end
function make_condition(s_y::String, order::AbstractVector{<:Int})
  s_x = map_y_to_x(s_y, order)
  return make_condition(s_x)
end

