using DataFrames

mutable struct LP <: Model
    data::DataFrame # input data
    treatment::Union{Symbol,Int}  # treatment variable in data
    p::Int  # number of lags to include
    horizons::AbstractVector{<:Int}  # horizon for LP
    include_constant::Bool  # whether to include a constant

    coeffs::AbstractArray{<:Number}  # estimated coefficients per horizon
    Y::AbstractArray{<:Number}  # LHS  per horizon
    X::AbstractMatrix{<:Number}  # RHS (same for each horizon)
    U::AbstractArray{<:Number}  # Residuals per horizon
    Yhat::AbstractArray{<:Number}  # Fitted values per horizon
end

function LP(
    data::DataFrame,
    treatment::Union{Symbol,Int},
    p::Int,
    horizons::Union{Int,AbstractVector{<:Int}};
    include_constant::Bool=true
)

    if !isa(horizons, AbstractVector)
        horizons = [horizons]
    end

    type = eltype(data[:, 1])
    T, k = size(data)
    n_horizons = length(horizons)
    idx_treatment = _find_variable_idx(treatment, data)

    U = Array{type,3}(undef, 0, 0, 0)
    Yhat = Array{type,3}(undef, 0, 0, 0)
    beta = Array{type,3}(undef, 0, 0, 0)

    mat_data = type.(Matrix(data))
    mat_data_lag = make_lag_matrix(mat_data, p)
    X = hcat(mat_data[(p+1):end, 1:idx_treatment], mat_data_lag[(p+1):end, :])
    if include_constant
        X = hcat(ones(T - p), X)
    end

    max_horizon = maximum(horizons)
    Y = make_lead_matrix(mat_data, max_horizon)
    Y = hcat(mat_data, Y)
    Y = Y[(p+1):end, :]
    Y = reshape(Y, T-p, k, :)
    Y = Y[:, :, horizons .+ 1]

    return LP(data, treatment, p, horizons, include_constant, beta, Y, X, U, Yhat)
end


