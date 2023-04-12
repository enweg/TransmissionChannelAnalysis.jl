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

function Base.:+(x::BayesianEstimated, y::BayesianEstimated)
    return BayesianEstimated(x.value + y.value, nothing)
end