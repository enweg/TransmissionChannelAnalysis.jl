# TODO: test the functions in this file
"""
Mediation through a single variable
This corresponds to the statement
through_t || through_{t+1} || ... || through_{t+h}
i.e. the effect that must go through the `through` 
variable at some horizon
"""
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
function through_a(from::Int, to::Int, through::Int, sirfs::Array{T, 3}) where {T<:Real}
    horizon = size(sirfs.irfs.value, 3)
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
function through_a(from::Int, to::Int, through::Int, sirfs::StructuralImpulseResponseFunction{E}) where {E<:FixedEstimated}
    return through_a(from, to, through, sirfs.irfs.value)
end
function through_a(from::Int, to::Int, through::Int, sirfs::StructuralImpulseResponseFunction{E}) where {E<:BayesianEstimated}
    f = x -> through_a(from, to, through, x)
    return map(f, sirfs.irfs)
end





function Base.map(f::Function, be::BayesianEstimated; meta=nothing)
    nd = ndims(be)
    return BayesianEstimated(Base.mapslices(f, be.value; dims=collect(1:(nd-2))), meta)
end

function Base.map(f::Function, bes::B...; meta=nothing) where {B<:BayesianEstimated}
    nd = ndims(bes[1])
    nchains = size(bes[1], nd)
    ndraws = size(bes[1], nd-1)
    colons_be = fill(:, nd-2)
    ret = f([be.value[colons_be..., 1, 1] for be in bes]...)
    returns = Array{eltype(ret)}(undef, size(ret)..., ndraws, nchains)
    colons_return = fill(:, ndims(returns)-2)
    for chain in 1:nchains
        for draw in 1:ndraws
            returns[colons_return..., draw, chain] = f([be.value[colons_be..., draw, chain] for be in bes]...)
        end
    end
    return BayesianEstimated(returns, meta)
end






