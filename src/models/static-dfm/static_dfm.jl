using Statistics

"""
```math
using Core: Argument
\\begin{split}
Y_t = \\Lambda F_t + \\eta_t \\\\
F_t = B(L)F_{t-1} + A_0^{-1}\\varepsilon_t
\\end{split}
```
"""
mutable struct DFM <: Model
    # Estimated quantities
    factor_var::Union{Nothing,VAR}
    Lambda::AbstractMatrix
    F::AbstractMatrix

    # Input
    input_data::DataFrame
    p::Int  # Number of lags in the VAR
    trend_exponents::AbstractVector{<:Number}  # trends for VAR
    r::Int  # Number of static factors

    # Intermediate data
    Y::AbstractMatrix      # inputs as matrix
    Yhat::AbstractMatrix   # Fitted values
    eta_hat::AbstractMatrix   # measurement error
    # reduced form VAR errors are in the factor_var
end

function DFM(
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

    return DFM(nothing, Lambda, F, data, p, trend_exponents, r, Y, Yhat, eta_hat)
end

is_fitted(model::DFM) = size(model.Yhat, 1) >= 1 &&
                        (!isnothing(model.factor_var) ||
                         is_fitted(model.factor_var))
is_structural(model::DFM) = false

"""
Returns (Lambda, VAR coeffs)
"""
function coeffs(model::DFM, exclude_deterministic::Bool=false)
    require_fitted(model)
    return (model.Lambda, coeffs(model.factor_var, exclude_deterministic))
end

fitted(model::DFM) = require_fitted(model) && model.Yhat
factors(model::DFM) = require_fitted(model) && model.F
loadings(model::DFM) = require_fitted(model) && model.Lambda
residuals(model::DFM) = require_fitted(model) &&
                        (model.eta_hat, residuals(model.factor_var))

nobs(model::DFM) = size(model.Y, 1)
get_dependent(model::DFM) = model.Y
get_input_data(model::DFM) = model.input_data
get_factor_var(model::DFM) = model.factor_var

#-------------------------------------------------------------------------------
# SIMULATION
#-------------------------------------------------------------------------------

function simulate!(
    ::Type{DFM},
    merrors::AbstractMatrix{<:Number},  # measurement errors
    varerrors::AbstractMatrix{<:Number},  # VAR errors
    B::AbstractMatrix{<:Number},   # VAR coefficient matrix
    Lambda::AbstractMatrix{<:Number};  # Factor loadings
    trend_exponents::AbstractVector{<:Number}=[0],  # VAR trend exponents
    initial_F::Union{Nothing,AbstractVector{<:Number}}=nothing
)

    factor_var = simulate!(
        VAR, varerrors, B;
        trend_exponents=trend_exponents, initial=initial_F
    )

    F = Matrix(get_input_data(factor_var))
    T = eltype(merrors)
    mul!(merrors, Lambda, F', one(T), one(T))

    data = DataFrame(merrors', "Y" .* string.(1:size(merrors, 1)))
    p = factor_var.p
    r = size(F, 2)
    return DFM(data, p, r; trend_exponents=trend_exponents)
end

function simulate(
    ::Type{DFM},
    T::Int,
    B::AbstractMatrix{<:Number},
    Lambda::AbstractMatrix{<:Number},
    Sigma_eta::AbstractMatrix{<:Number}=I(size(Lambda, 1)),
    Sigma_u::AbstractMatrix{<:Number}=I(size(B, 1));
    trend_exponents::AbstractVector{<:Number}=[0],
    initial_F::Union{Nothing,AbstractVector{<:Number}}=nothing
)

    r = size(B, 1)
    N = size(Lambda, 1)
    merrors = cholesky(Sigma_eta).L * randn(N, T)
    varerrors = cholesky(Sigma_u).L * randn(r, T)
    return simulate!(
        DFM, merrors, varerrors, B, Lambda;
        trend_exponents=trend_exponents, initial_F=initial_F
    )
end

#-------------------------------------------------------------------------------
# ESTIMATION FUNCTIONS
#-------------------------------------------------------------------------------

"""
    fit!(model::DFM; atol::Real=sqrt(eps()))

Fit a Dynamic Factor Model (DFM) with static factors.

The data must have been centered and scaled beforehand. This can
be achieved by creating the `DFM` with keyword arguments `center=true` and
`scale=true`.

# Arguments
- `model::DFM`: A Dynamic Factor Model with static factors

# Keyword Arguments
- `atol::Real`: Absolute tolerance for the scale and center check

# Notes
- The DFM is estimated in static form using `PCA`.
"""
function fit!(model::DFM; atol::Real=sqrt(eps()))
    Y = get_dependent(model)
    all(isapprox.(std(Y; dims=1), 1.0; atol=atol)) ||
        throw(ArgumentError("Data must be scaled first."))
    all(isapprox.(mean(Y; dims=1), 0.0; atol=atol)) || throw(ArgumentError("Data must be centered first."))

    pca = PCA(Y, model.r)
    F = factors(pca)
    model.F = F
    Lambda = loadings(pca)
    model.Lambda = Lambda
    model.Yhat = fitted(pca)
    model.eta_hat = residuals(pca)

    df_F = DataFrame(F, "F" .* string.(1:model.r))
    factor_var = VAR(df_F, model.p; trend_exponents=model.trend_exponents)
    fit!(factor_var)
    model.factor_var = factor_var

    return model
end
