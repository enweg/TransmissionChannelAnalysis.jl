"""
    through(idx, horizons, order) --> Q

All paths must go through variables in `idx` in periods `horizons`. Note, this 
uses the dynamic system notation `y` such that `idx` refers to the index of 
the variable in the original dynamic system, i.e. in the SVARMA.

## Arguments

**For the single variable version**:

- `idx::Int`: Index of variable through which the paths should go. This is the  
  original index in the dynamic system, e.g. the SVAR, before applying the 
  transmission matrix. 
- `horizons::AbstractVector{<:Int}`: Horizons for which the paths must go through 
  the variable. 
- `order::AbstractVector{<:Int}`: Variable ordering determined by the transmission 
  matrix

**For the multiple variable version**:

- `idx::AbstractVector{<:Int}`: Indices of variables through which the paths 
  should go. These are the original indices in the dynamic system, e.g. the SVAR, 
  before applying the transmission matrix. 
- `horizons::Union{AbstractVector{<:Int},Vector{AbstractVector{<:Int}}}`: Horizons 
  for which the paths must go through the variable. Must either be a vector 
  for each variable in `idx` or a single vector. If it is a single vector, then 
  the horizons will be applied to each variable in `idx`.
- `order::AbstractVector{<:Int}`: Variable ordering determined by the transmission 
  matrix

## Returns
- Returns a transmission condition `Q`. 

## Notes
- The transmission effect can be calculated using `transmission`. 

## Examples

Suppose we are interested in the contemporaneous channel in Section 5.1 of 
$WEGNER, i.e. we are interested in the effect going through the federal funds
rate contemporaneously. In our estimated model, the federal funds rate is the
first variable. We would thus define the contemporaneous channel as 

```julia
contemporaneous_channel = through(1, [0], 1:4)  # we have four variables with ffr being first
```

More generally, we could define the following, which would be the effect through
the federal funds rate in the first two periods. 

```julia
q = through(1, [0, 1], 1:4)
```

The following is also allowed which is the effect through the federal funds 
rate contemporaneously and one period later, and through the output gap
contemporaneously and one period later, where the federal funds rate is ordered 
first and the output gap second. 

```julia
q = through([1, 2], [[0, 1], [0, 1]], 1:4)
```

If we want to re-order the federal funds rate after the output gap, we simply
change the ordering to `[2, 1, 3, 4]` where `2` corresponds to the output gap
in the original ordering. 

```julia
q = through([1, 2], [[0, 1], [0, 1]], [2, 1, 3, 4])
```
"""
function through(
    idx::Int,
    horizons::AbstractVector{<:Int},
    order::AbstractVector{<:Int}
)

    s = join(["y_{$idx, $h}" for h in horizons], " & ")
    return make_condition(s, order)
end
function through(
    idx::AbstractVector{<:Int},
    horizons::Vector{<:AbstractVector{<:Int}},
    order::AbstractVector{<:Int}
)

    qs = [through(i, h, order) for (i, h) in zip(idx, horizons)]
    q = qs[1]
    for i = 2:lastindex(qs)
        q = q & qs[i]
    end
    return q
end
function through(
    idx::AbstractVector{<:Int},
    horizons::AbstractVector{<:Int},
    order::AbstractVector{<:Int}
)

    qs = [through(i, horizons, order) for i in idx]
    q = qs[1]
    for i = 2:lastindex(qs)
        q = q & qs[i]
    end
    return q
end

"""
    not_through(idx, horizons, order) --> Q

All paths cannot go through variables in `idx` in periods `horizons`. Note, this 
uses the dynamic system notation `y` such that `idx` refers to the index of 
the variable in the original dynamic system, i.e. in the SVARMA.


## Arguments

**For the single variable version**:

- `idx::Int`: Index of variable through which the paths cannot go. This is the  
  original index in the dynamic system, e.g. the SVAR, before applying the 
  transmission matrix. 
- `horizons::AbstractVector{<:Int}`: Horizons for which the paths cannot go through 
  the variable. 
- `order::AbstractVector{<:Int}`: Variable ordering determined by the transmission 
  matrix

**For the multiple variable version**:

- `idx::AbstractVector{<:Int}`: Indices of variables through which the paths 
  cannot go. These are the original indices in the dynamic system, e.g. the SVAR, 
  before applying the transmission matrix. 
- `horizons::Union{AbstractVector{<:Int},Vector{AbstractVector{<:Int}}}`: Horizons 
  for which the paths cannot go through the variable. Must either be a vector 
  for each variable in `idx` or a single vector. If it is a single vector, then 
  the horizons will be applied to each variable in `idx`.
- `order::AbstractVector{<:Int}`: Variable ordering determined by the transmission 
  matrix

## Returns
- Returns a transmission condition `Q`. 

## Notes
- The transmission effect can be calculated using `transmission`. 

## Examples

The non-contemporaneous channel of monetary policy is defined in Section 5.1 of 
$WEGNER as the effect not going through a contemporaneous adjustment of the 
federal funds rate, where the transmission matrix orders the federal funds rate 
first. Thus, if the original SVAR has the federal funds rate ordered first, the 
non-contemporaneous effect can be obtained in the following way. 

```julia
q = not_through(1, [0], 1:4)
```

Similarly, $WEGNER define the anticipation channel of government defense spending
as the effect not going through government defense spending. With government 
defense spending ordered second in the VAR, the following can be used to obtain 
the anticipation channel. 

```{julia}
q = not_through(2, 0:20, 1:4)
```

"""
function not_through(
    idx::Int,
    horizons::AbstractVector{<:Int},
    order::AbstractVector{<:Int}
)

    s = join(["!y_{$idx, $h}" for h in horizons], " & ")
    return make_condition(s, order)
end
function not_through(
    idx::AbstractVector{<:Int},
    horizons::Vector{<:AbstractVector{<:Int}},
    order::AbstractVector{<:Int}
)

    qs = [not_through(i, h, order) for (i, h) in zip(idx, horizons)]
    q = qs[1]
    for i = 2:lastindex(qs)
        q = q & qs[i]
    end
    return q
end
function not_through(
    idx::AbstractVector{<:Int},
    horizons::AbstractVector{<:Int},
    order::AbstractVector{<:Int}
)
    qs = [not_through(i, horizons, order) for i in idx]
    q = qs[1]
    for i = 2:lastindex(qs)
        q = q & qs[i]
    end
    return q
end
