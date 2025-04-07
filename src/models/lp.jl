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

    idx_treatment = _find_variable_idx(treatment, data)
    X, Y = _create_LP_XY(data, idx_treatment, p, horizons, include_constant)

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

function _2sls(
    X::AbstractMatrix{<:Number},
    Y::AbstractMatrix{<:Number},
    Z::AbstractMatrix{<:Number}
)

    X_hat = Z * ((Z' * Z) \ (Z' * X))
    return (X_hat' * X_hat) \ (X_hat' * Y)
end
# Enforcing that instrument are ordered before treatment 
# This makes it easy to get Zt because we can just extract it from the 
# contemporaneous block in Xt
# Current implementation allows for that some contemporaneous controls are 
# not used as instruments for the treatment
function fit!(model::LP, method::ExternalInstrument)
    idxs_instruments = _find_variable_idx.(method.instruments, [model.data])
    idx_treatment = _find_variable_idx(model.treatment, model.data)
    idx_treatment_instrument = _find_variable_idx(method.treatment, model.data)
    all(idxs_instruments .< idx_treatment) || error("Instruments must come before treatment variable in data.")
    idx_treatment == idx_treatment_instrument || error("LP and ExternalInstrument treatment variable differ.")

    m = model.include_constant
    # X = [constant contemporaneous lags]
    # Z can be extracted from contemporaneous
    Z = model.X[:, idxs_instruments.+m]
    # Now exclude them from X
    X = model.X[:, filter(x -> !(x in idxs_instruments .+ m), 1:size(model.X, 2))]
    # the treatment index also changes
    # since treatment is always the last contemporanous variable 
    # all we need to know is how many instruments we remove from that block
    model.treatment = idx_treatment - length(idxs_instruments)

    if method.normalising_horizon > 0
        # lead the treatment column in X to adjust for which horizon 
        # the unit effect normalisation applies to
        nlead = method.normalising_horizon
        X[:, model.treatment+m] .= make_lead_matrix(X[:, model.treatment+m], nlead)
        # remove NaNs at end of data
        X = X[1:(end-nlead), :]
        model.Y = model.Y[1:(end-nlead), :, :]
        Z = Z[1:(end-nlead), :]
    end
    model.X = X

    # Adding all other exogenous variables to Z
    # these are all variables that are not the treatment variable
    Z = hcat(Z, X[:, filter(!=(model.treatment + m), 1:size(X, 2))])

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

    irfs = coeffs(model, true)[:, model.treatment:model.treatment, :]
    return irfs
end

function _identify_irfs(model::LP, method::ExternalInstrument, max_horizon::Int)
    model.horizons == 0:max_horizon || error("LP horizons do not match IRF horizons.")
    is_fitted(model) || fit!(model, method)

    irfs = coeffs(model, true)[:, model.treatment:model.treatment, :]
    return irfs
end

function IRF(model::LP, method::AbstractIdentificationMethod, max_horizon::Int)
    irfs = _identify_irfs(model, method, max_horizon)
    varnames = Symbol.(names(get_input_data(model)))
    return IRF(irfs, varnames, model, method)
end

