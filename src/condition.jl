"""
    make_condition(s::String)

Given a `s::String`, form a transmission condition of the form Q(b) where b is a
Boolean statement. 

Transmission mechanisms are described using Boolean statements involving the 
variables in the system. Each variable starts with `x` followed by a number.
For example, given a three variable VAR(1), y_{1,t} -> x_1, y_{2, t} -> x_2, 
y_{3, t} -> x_3, y_{1, t+1} -> x_4, y_{2, t+1} -> x_5, ... Boolean statements 
then involve expressions in the `x` variables and define which paths can be
taken. Each path involved in the transmission mechanism must satisfy the Boolean
statement. 

## Arguments

- `s::String`: A Boolean statement given as a string. Variables must start
  with `x` for them to be valid variables. 

## Returns 

- Returns a transmission condition. See also [`Q`](@ref).

## Examples

```julia
s = "x2 & !x3"
cond = make_condition(s)
```

## Notes

- Ensure that variables in the input string are correctly formatted as `x`
  followed by a number. 
- The resulting transmission condition can be used in [`transmission`](@ref) to
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