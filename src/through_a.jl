function _Q_through_a(
    from::Int, 
    to::Int, 
    through::Int, 
    sirfs::Array{T, 3}, 
    horizon::Int, 
    c::Array{Int}
) where {T<:Real}

    c = sort!(c)
    meffect = sirfs[through, from, c[1]+1]
    for i in 2:lastindex(c)
        meffect *= sirfs[through, through, c[i]-c[i-1]+1]
    end
    meffect *= sirfs[to, through, horizon-c[end]+1]
    return meffect
end
"""
    through_a(from::Int, to::Int, through::Int, sirfs::Array{<:Real, 3})
    through_a(from::Int, to::Int, through::Int, sirfs::StructuralImpulseResponseFunction{<:FixedEstimated, M, N})
    through_a(from::Int, to::Int, through::Int, sirfs::StructuralImpulseResponseFunction{<:BayesianEstimated, M, N})

Mediation through a single variable This corresponds to the statement through_t
|| through_{t+1} || ... || through_{t+h} i.e. the effect that must go through
the `through` variable at some horizon. Uses the 'algebra' of path specific
effects. 

## Arguments

- `from::Int`: Variable number of the origin variable
- `to::Int`: Variable number of the destination variable
- `through::Int`: Variable number of the mediating variable
- `sirfs`: Structural IRFs of dimension n×n×H where n is the
  number of variables, and H is the number of horizons, including horizon zero.
  The structural IRFs are of the form `sirfs[to, from, horizon]`. 
"""
function through_a(from::Int, to::Int, through::Int, sirfs::Array{T, 3}) where {T<:Real}
    horizon = size(sirfs, 3)
    meffect = zeros(horizon)
    for h in 0:horizon-1
        combs = combinations(0:h)
        for c in combs
            k = length(c)
            meffect[h+1] += (-1)^(k+1)*_Q_through_a(from, to, through, sirfs, h, c)
        end
    end
    return meffect
end
function through_a(from::Int, to::Int, through::Int, sirfs::StructuralImpulseResponseFunction{E, M, N}) where {E<:FixedEstimated, M, N}
    return through_a(from, to, through, sirfs.irfs.value)
end
function through_a(from::Int, to::Int, through::Int, sirfs::StructuralImpulseResponseFunction{E, M, N}) where {E<:BayesianEstimated, M, N}
    f = x -> through_a(from, to, through, x)
    return map(f, sirfs.irfs)
end
