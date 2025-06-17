using DataFrames

"""
    _create_LP_XY(data::DataFrame, treatment::Int, p::Int,
                  horizons::AbstractVector{<:Int},
                  include_constant::Bool=true)

Constructs the regressor matrix `X` and the response array `Y` for a local
projection (LP) setup.

Local projection regressions are of the form:

```math
    y^h = X\beta + \varepsilon^h
```

where `h` indexes the forecast horizon. The same regressor matrix `X` is used
for all horizons, while the dependent variables `y^h` vary by horizon. This
function stacks these horizon-specific dependent variables along the third
dimension of the returned array `Y`, such that its shape is:
(observations, variables, horizons))

# Arguments
- `data::DataFrame`: the input time series dataset
- `treatment::Int`: index of the treatment variable (used for selecting
  contemporaneous controls)
- `p::Int`: the lag length applied uniformly across all horizons
- `horizons::AbstractVector{<:Int}`: the set of forecast horizons to include
- `include_constant::Bool`: whether to include a constant column in `X`

# Returns
- `X`: a matrix of regressors, common across horizons
- `Y`: a 3-dimensional array of outcomes with horizon `h` stacked along the
  third dimension

# Details
- All variables in columns before the `treatment` variable are used as
  contemporaneous controls.
- Regressors in `X` include (optionally) a contant, contemporaneous controls,
  treatment variable, and lagged values of all variables.
- Dependent variables in `Y` include future values (leads) for all variables
  at each horizon in `horizons`.
"""
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

"""
    LP <: Model

Local Projection (LP) model for estimating impulse response functions (IRFs)
in a flexible and semi-parametric manner.

Each LP regression estimates the dynamic response of an outcome variable at
future horizon `h` to a one-period change in a treatment variable at time `t`,
controlling for contemporaneous and lagged covariates.

The regression model is specified as:

```math
    w_{i,t+h} = \\mu_{i,h} + \\theta_{i,h} x_t + \\gamma_{i,h}' r_t +
                \\sum_{l=1}^p \\delta_{i,h,l} w_{t-l} + \\xi_{i,h,t}
```

where ``w_t = (r_t', x_t, q_t')`` and:
- ``x_t`` is the treatment variable
- ``r_t`` contains contemporaneous controls (all variables before `x_t`)
- ``p`` is the number of lags included
- ``\\theta_{i,h}`` is the relative IRF of `x_t` on the `i`-th variable at
  horizon ``h``.

The treatment variable may be endogenous. Structural interpretation of IRFs
can be achieved using valid instruments—see `ExternalInstrument` for one such
method. If the treatment satisfies a conditional ignorability assumption
(a recursive assumption in macro), then the coefficient has a structural
interpretation even without the use of instruments. For this to hold,
``x_t - E(x_t|r_t, w_{t-1}, ..., w_{t-p})`` must be equal to the structural shock.

# Fields
- `data::DataFrame`: the dataset containing the time series
- `treatment::Union{Symbol, Int}`: column indicating the treatment variable
- `p::Int`: number of lags to include
- `horizons::AbstractVector{<:Int}`: forecast horizons for the projections
- `include_constant::Bool`: whether to include an intercept in each regression
- `coeffs::AbstractArray{<:Number}`: coefficient estimates per horizon
- `Y::AbstractArray{<:Number}`: dependent variables for each horizon
- `X::AbstractMatrix{<:Number}`: common regressor matrix
- `U::AbstractArray{<:Number}`: residuals per horizon
- `Yhat::AbstractArray{<:Number}`: fitted values per horizon
"""
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

"""
    LP(data::DataFrame,
       treatment::Union{Symbol, Int},
       p::Int,
       horizons::Union{Int, AbstractVector{<:Int}};
       include_constant::Bool=true)

Constructs a `LP` model object for estimating local projections with specified
`horizons`, `treatment` variable, and lag length `p`.

All variables before the treatment variable in the dataset are used as
contemporaneous controls.

# Arguments
- `data`: the time series dataset
- `treatment`: the treatment variable (column index or name)
- `p`: number of lags to include
- `horizons`: forecast horizons to compute IRFs for

## Keyword Arguments
- `include_constant`: whether to include a constant in the regressors
"""
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

    idx_treatment = _find_variable_idx(treatment, data)
    X, Y = _create_LP_XY(data, idx_treatment, p, horizons, include_constant)

    return LP(data, treatment, p, horizons, include_constant, coeffs, Y, X, U, Yhat)
end

"""
    coeffs(model::LP, exclude_deterministic::Bool=false)

Returns the coefficient estimates of the local projections. The returned
object is three dimensional with the coefficients for each horizon stacked
along the third dimension. If `exclude_deterministic` is `true`, then
all coefficients for deterministic variables, i.e. the coefficient on the
constant, will be removed.

Coefficients for each horizon have the ordering
[constant, contemporaneous controls, treatment, lagged controls]
"""
function coeffs(model::LP, exclude_deterministic::Bool=false)
    require_fitted(model)
    exclude_deterministic || return model.coeffs

    return model.coeffs[:, (model.include_constant+1):end, :]
end

"""
    fitted(model::LP)

Returns the fitted values of the local projection in a three dimensional array.
The fitted values for each horizon are stacked along the third dimension.
The second dimension corresponds to each variable.
"""
fitted(model::LP) = require_fitted(model) && model.Yhat

"""
    residuals(model::LP)

Returns the residuals of the local projection in a three dimensional array, where
the third dimension corresponds to the various horizons, and the second
dimension corresponds to the various variables.
"""
residuals(model::LP) = require_fitted(model) && model.U

nobs(model::LP) = size(model.data, 1) - model.p .- model.horizons

"""
    get_dependent(model::LP)

The dependent variables are returned as a three dimensional array, with the
third dimension corresponding the the various forecast horizons.
"""
get_dependent(model::LP) = model.Y

get_independent(model::LP) = model.X
get_input_data(model::LP) = model.data
is_fitted(model::LP) = size(model.coeffs, 1) > 0

"""
    is_structural(model::LP)

Always returns `true`. We assume here that local projections are only used
to estimate structural IRFs. This is not always true. `is_structural` should
thus be used with care in the case of local projections.
"""
is_structural(model::LP) = true  # we just assume that this is always the case

function Base.show(io::IO, ::MIME"text/plain", x::LP)
    treatment = isa(x.treatment, Int) ? names(x.data)[x.treatment] : x.treatment
    s = """
        Local Projection
        ================
        variables: $(join(names(x.data), ", "))
        treatment variable: $(treatment)
        p: $(x.p)
        constant: $(x.include_constant)
        horizons: $(x.horizons)
        """
    println(io, s)
end

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

function fit!(model::LP, method::ExternalInstrument)
    idx_treatment = _find_variable_idx(model.treatment, model.data)
    idx_treatment_instrument = _find_variable_idx(method.treatment, model.data)
    idx_treatment == idx_treatment_instrument || error("LP and ExternalInstrument treatment variable differ.")
    size(method.instruments, 1) == size(get_input_data(model), 1) || error("External instruments do not span the same period as model data.")

    Z = method.instruments[(model.p+1):end, :]
    m = model.include_constant
    if method.normalising_horizon > 0
        # lead the treatment column in X to adjust for which horizon
        # the unit effect normalisation applies to
        nlead = method.normalising_horizon
        model.X[:, idx_treatment+m] .= make_lead_matrix(model.X[:, idx_treatment+m], nlead)
        # remove NaNs at end of data
        model.X = model.X[1:(end-nlead), :]
        model.Y = model.Y[1:(end-nlead), :, :]
        Z = Z[1:(end-nlead), :]
    end

    # Adding all other exogenous variables to Z
    # these are all variables that are not the treatment variable
    Z = hcat(Z, model.X[:, setdiff(1:size(model.X, 2), idx_treatment + m)])

    type = eltype(model.X)
    k = size(model.Y, 2)
    num_coeffs = size(model.X, 2)
    coeffs = Array{type,3}(undef, k, num_coeffs, length(model.horizons))
    Yhat = fill(type(NaN), size(model.Y))
    U = fill(type(NaN), size(model.Y))

    for (i, h) in enumerate(model.horizons)
        X = model.X[1:(end-h), :]
        Y = model.Y[1:(end-h), :, i]
        Zi = Z[1:(end-h), :]
        coeffs[:, :, i] .= _2sls(X, Y, Zi)'

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
    model_best = LP(
        get_input_data(model),
        model.treatment,
        model_var.p,
        model.horizons;
        include_constant=model.include_constant
    )
    return model_best, ic_table
end


#-------------------------------------------------------------------------------
# IMPULSE RESPONSE FUNCTIONS
#-------------------------------------------------------------------------------

function _identify_irfs(model::LP, ::Recursive, max_horizon::Int)
    model.horizons == 0:max_horizon || error("LP horizons do not match IRF horizons.")
    is_fitted(model) || fit!(model, Recursive())

    data = get_input_data(model)
    k = size(data, 2)
    idx_treatment = _find_variable_idx(model.treatment, data)
    irfs = fill(NaN, k, k, max_horizon + 1)
    irfs[:, idx_treatment:idx_treatment, :] = coeffs(model, true)[:, idx_treatment:idx_treatment, :]
    return irfs
end

function _identify_irfs(model::LP, method::ExternalInstrument, max_horizon::Int)
    model.horizons == 0:max_horizon || error("LP horizons do not match IRF horizons.")
    is_fitted(model) || fit!(model, method)

    data = get_input_data(model)
    idx_treatment = _find_variable_idx(model.treatment, data)
    k = size(data, 2)
    irfs = fill(NaN, k, k, max_horizon + 1)
    irfs[:, idx_treatment:idx_treatment, :] = coeffs(model, true)[:, idx_treatment:idx_treatment, :]
    return irfs
end

function IRF(model::LP, method::AbstractIdentificationMethod, max_horizon::Int)
    irfs = _identify_irfs(model, method, max_horizon)
    varnames = Symbol.(names(get_input_data(model)))
    return IRF(irfs, varnames, model, method)
end
