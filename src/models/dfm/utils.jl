"""
Convert a lag polynomial `[X_1, X_2, X_3, ...]` to a matrix that 
stacks the lag-matrices horizontally so that it returns `[X_1 X_2 ... ]`
"""
function lag_poly_to_matrix(lag_poly::AbstractVector{<:AbstractMatrix})
    return reduce(hcat, lag_poly)
end

