using Statistics

using Base: _NAMEDTUPLE_NAME
mutable struct SDFM <: Model
    # Estimated quantities
    factor_svar::Union{Nothing,SVAR}
    Lambda::AbstractMatrix
    F::AbstractMatrix

    # Input
    input_data::DataFrame
    p::Int  # Number of lags in the VAR
    trend_exponents::AbstractVector{<:Number}  # trends for VAR
    r::Int  # Number of static factors
    scale::Bool  # whether data has been scaled
    center::Bool # whether data has been centered

    # Intermediate data
    Y::AbstractMatrix      # inputs as matrix
    Yhat::AbstractMatrix   # Fitted values
    eta_hat::AbstractMatrix   # measurement error
    # reduced form VAR errors are in the factor_var
end

function SDFM(
    data::DataFrame, p::Int, r::Int;
    trend_exponents::AbstractVector{<:Number}=[0],
    scale::Bool=true,
    center::Bool=true
)

    type = eltype(data[:, 1])
    Lambda = Matrix{type}(undef, 0, 0)
    F = Matrix{type}(undef, 0, 0)
    Yhat = Matrix{type}(undef, 0, 0)
    eta_hat = Matrix{type}(undef, 0, 0)

    Y = Matrix(data)
    if center
        Y .-= mean(Y; dims=1)
    end
    if scale
        Y ./= std(Y; dims=1)
    end

    return SDFM(
        nothing, Lambda, F, data, p, trend_exponents, r, scale, center, Y,
        Yhat, eta_hat
    )
end

is_fitted(model::SDFM) = size(model.Yhat, 1) >= 1 &&
                         (!isnothing(model.factor_svar) ||
                          is_fitted(model.factor_svar))
is_structural(model::SDFM) = true
is_stable(model::SDFM) = require_fitted(model) && is_stable(model.factor_svar)
is_scaled(model::SDFM) = model.scale
is_centered(model::SDFM) = model.center

"""
Returns (Lambda, SVAR coeffs)
"""
function coeffs(model::SDFM, exclude_deterministic::Bool=false)
    require_fitted(model)
    return (model.Lambda, coeffs(model.factor_svar, exclude_deterministic))
end

fitted(model::SDFM) = require_fitted(model) && model.Yhat
factors(model::SDFM) = require_fitted(model) && model.F
loadings(model::SDFM) = require_fitted(model) && model.Lambda
residuals(model::SDFM) = require_fitted(model) &&
                         (model.eta_hat, residuals(model.factor_svar))

shocks(model::SDFM) = require_fitted(model) && shocks(model.factor_svar)

nobs(model::SDFM) = size(model.Y, 1)
get_dependent(model::SDFM) = model.Y
get_input_data(model::SDFM) = model.input_data
get_factor_svar(model::SDFM) = model.factor_svar

#-------------------------------------------------------------------------------
# SIMULATION
#-------------------------------------------------------------------------------

function simulate!(
    ::Type{SDFM},
    merrors::AbstractMatrix{<:Number},
    varshocks::AbstractMatrix{<:Number},
    A0::AbstractMatrix{<:Number},
    A_plus::AbstractMatrix{<:Number},
    Lambda::AbstractMatrix{<:Number};
    trend_exponents::AbstractVector{<:Number}=[0],
    initial_F::Union{Nothing,AbstractVector{<:Number}}=nothing
)

    factor_svar = simulate(
        SVAR, varshocks, A0, A_plus;
        trend_exponents=trend_exponents, initial=initial_F
    )

    F = Matrix(get_input_data(factor_svar))
    T = eltype(merrors)
    mul!(merrors, Lambda, F', one(T), one(T))

    data = DataFrame(merrors', "Y" .* string.(1:size(merrors, 1)))
    p = factor_svar.p
    r = size(F, 2)

    model = SDFM(data, p, r; trend_exponents=trend_exponents)
    model.F = F
    return model
end

function simulate(
    ::Type{SDFM},
    T::Int,
    A0::AbstractMatrix{<:Number},
    A_plus::AbstractMatrix{<:Number},
    Lambda::AbstractMatrix{<:Number},
    Sigma_eta::AbstractMatrix{<:Number}=I(size(Lambda, 1));
    trend_exponents::AbstractVector{<:Int}=[0],
    initial_F::Union{Nothing,AbstractVector{<:Int}}=nothing
)

    N = size(Lambda, 1)
    r = size(A0, 1)

    merrors = cholesky(Sigma_eta).L * randn(N, T)
    varshocks = randn(r, T)
    return simulate!(
        SDFM, merrors, varshocks, A0, A_plus, Lambda;
        trend_exponents=trend_exponents, initial_F=initial_F
    )
end

#-------------------------------------------------------------------------------
# ESTIMATION AND IDENTICATION FUNCTIONS
#-------------------------------------------------------------------------------

function fit!(
    model::SDFM, identification_method::AbstractIdentificationMethod;
    named_factors::Union{Nothing,AbstractVector{<:Int}}=1:model.r
)

    data = get_input_data(model)
    scale = is_scaled(model)
    center = is_centered(model)
    model_dfm = DFM(
        data, model.p, model.r;
        trend_exponents=model.trend_exponents, scale=scale, center=center
    )
    fit!(model_dfm; named_factors=named_factors)

    factor_var = get_factor_var(model_dfm)
    factor_svar = identify(factor_var, identification_method)

    model.factor_svar = factor_svar
    model.Lambda = loadings(model_dfm)
    model.F = factors(model_dfm)
    model.Yhat = fitted(model_dfm)
    model.eta_hat, _ = residuals(model_dfm)
    return model
end

# TODO: implement fit_and_select!

#-------------------------------------------------------------------------------
# IMPULSE RESPONSE FUNCTIONS
#-------------------------------------------------------------------------------

function IRF(model::SDFM, max_horizon::Int)
    require_fitted(model)

    factor_svar = get_factor_svar(model)
    factor_irfs = IRF(factor_svar, max_horizon).irfs
    Lambda = loadings(model)
    variable_irfs = mapslices(x -> Lambda * x, factor_irfs; dims=(1, 2))
    varnames = Symbol.(names(get_input_data(model)))
    factornames = Symbol.("F" .* string.(1:model.r))
    return IRF(variable_irfs, varnames, model), IRF(factor_irfs, factornames, model)
end

function IRF(model::DFM, method::AbstractIdentificationMethod, max_horizon::Int)
    require_fitted(model)

    factor_var = get_factor_var(model)
    factor_irfs = IRF(factor_var, method, max_horizon).irfs
    Lambda = loadings(model)
    variable_irfs = mapslices(x -> Lambda * x, factor_irfs; dims=(1, 2))
    varnames = Symbol.(names(get_input_data(model)))
    factornames = Symbol.("F" .* string.(1:model.r))
    return IRF(variable_irfs, varnames, model), IRF(factor_irfs, factornames, model)
end
