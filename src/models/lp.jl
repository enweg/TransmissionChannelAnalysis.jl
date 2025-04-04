using DataFrames

function _create_LP_XY(
    data::DataFrame, 
    treatment::Int,
    p::Int, 
    horizons::AbstractVector{<:Int}, 
    include_constant::Bool=true 
)

    type = eltype(data[:, 1])
    T, k = size(data)
    idx_treatment = _find_variable_idx(treatment, data)

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
    Y = reshape(Y, T - p, k, :)
    Y = Y[:, :, horizons.+1]

    return X, Y
end

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

    U = Array{type,3}(undef, 0, 0, 0)
    Yhat = Array{type,3}(undef, 0, 0, 0)
    coeffs = Array{type,3}(undef, 0, 0, 0)

    X, Y = _create_LP_XY(data, treatment, p, horizons, include_constant)

    return LP(data, treatment, p, horizons, include_constant, coeffs, Y, X, U, Yhat)
end

function coeffs(model::LP, exclude_deterministic::Bool=false)
    require_fitted(model)
    exclude_deterministic || return model.coeffs

    return model.coeffs[:, (model.include_constant+1):end, :]
end
fitted(model::LP) = require_fitted(model) && model.Yhat
residuals(model::LP) = require_fitted(model) && model.U
nobs(model::LP) = size(model.data, 1) - model.p .- model.horizons
get_dependent(model::LP) = model.Y
get_independent(model::LP) = model.X
get_input_data(model::LP) = model.data
is_fitted(model::LP) = size(model.coeffs, 1) > 0
is_structural(model::LP) = true  # we just assume that this is always the case

#-------------------------------------------------------------------------------
# INFORMATION CRITERIA
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# ESTIMATION FUNCTIONS
# Recursive
# LP-IV
#-------------------------------------------------------------------------------

# make this default because LP is anyways a recursive design
function fit!(model::LP, method::AbstractIdentificationMethod=Recursive())
    fit!(model, method)
end
function fit!(model::LP, ::Recursive)
    type = eltype(model.X)
    k = size(model.Y, 2)
    num_coeffs = size(model.X, 2)
    coeffs = Array{type,3}(undef, k, num_coeffs, length(model.horizons))
    Yhat = fill(type(NaN), size(model.Y))
    U = fill(type(NaN), size(model.Y))

    for (i, h) in enumerate(model.horizons)
        X = model.X[1:(end-h), :]
        Y = model.Y[1:(end-h), :, i]
        coeffs[:, :, i] .= Y' * X / (X' * X)

        Yhat[1:(end-h), :, i] .= X * coeffs[:, :, i]'
        U[1:(end-h), :, i] .= Y - Yhat[1:(end-h), :, i]
    end

    model.coeffs = coeffs
    model.Yhat = Yhat
    model.U = U

    return model
end

# there is no good way to select the lag-length yet besides running an 
# auxiliary VAR
function fit_and_select!(model::LP, ::Recursive, ic_function::Function=aic)
    trend_exponents = model.include_constant ? [0] : Real[]
    model_var = VAR(get_input_data(model), model.p; trend_exponents=trend_exponents)
    model_var, ic_table = fit_and_select!(model_var, ic_function)

    model.p = model_var.p
end


#-------------------------------------------------------------------------------
# IMPULSE RESPONSE FUNCTIONS
#-------------------------------------------------------------------------------

function _identify_irfs(model::LP, ::Recursive, max_horizon::Int)
    model.horizons == 0:max_horizon || error("LP horizons do not match IRF horizons.")
    is_fitted(model) || fit!(model, Recursive())

    irfs = coeffs(model, true)[:, model.treatment:model.treatment, :]
    return irfs
end

function IRF(model::LP, method::AbstractIdentificationMethod, max_horizon::Int)
    irfs = _identify_irfs(model, method, max_horizon)
    varnames = Symbol.(names(get_input_data(model)))
    return IRF(irfs, varnames, model, method)
end
