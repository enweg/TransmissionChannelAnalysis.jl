mutable struct SVAR <: Model
    # defining the DGP
    A_plus::AbstractMatrix{<:Number}              # Coefficient matrices [C, A1, ..., Ap]
    A0::AbstractMatrix{<:Number}                  # Contemporaneous relationships
    p::Int                                        # Lag order
    trend_exponents::AbstractVector{<:Number}     # Determinstic trend as exponents for time trend

    # reduced form model
    var::VAR
end
function SVAR(data::DataFrame, p::Int; trend_exponents::AbstractVector{<:Number}=[0])
    K = size(data, 2)
    type = eltype(data[:, 1])

    # Making missing matrix of zero size
    A_plus = Matrix{type}(undef, 0, 0)
    A0 = Matrix{type}(undef, 0, 0)

    # constructing reduced form model
    var = VAR(data, p; trend_exponents=trend_exponents)
    return SVAR(A_plus, A0, p, trend_exponents, var)
end

function coeffs(model::SVAR, exclude_deterministic::Bool=false)
    require_fitted(model)
    exclude_deterministic || return (model.A0, model.A_plus)

    m = length(model.trend_exponents)
    return (model.A0, model.A_plus[:, (m+1):end])
end
fitted(model::SVAR) = require_fitted(model) && fitted(model.var)
residuals(model::SVAR) = require_fitted(model) && residuals(model.var)
shocks(model::SVAR) = require_fitted(model) && residuals(model.var) * model.A0'
nobs(model::SVAR) = nobs(model.var)
get_dependent(model::SVAR) = get_dependent(model.var)
get_independent(model::SVAR) = get_independent(model.var)
get_input_data(model::SVAR) = get_input_data(model.var)
is_fitted(model::SVAR) = size(model.A0, 1) >= 1
is_structural(model::SVAR) = true

make_companion_matrix(model::SVAR) = require_fitted(model) && make_companion_matrix(model.var)
spectral_radius(model::SVAR) = require_fitted(model) && spectral_radius(model.var)
is_stable(model::SVAR) = require_fitted(model) && is_stable(model.var)

aic(model::SVAR) = require_fitted(model) && aic(model.var)
hqc(model::SVAR) = require_fitted(model) && hqc(model.var)
sic(model::SVAR) = require_fitted(model) && sic(model.var)
bic(model::SVAR) = require_fitted(model) && bic(model.var)

#-------------------------------------------------------------------------------
# ESTIMATION AND IDENTICATION FUNCTIONS
#-------------------------------------------------------------------------------

function _identify(model::VAR, ::Recursive)
    is_fitted(model) || error("Reduced form model must first be estimated.")
    Sigma_u = cov(model)
    L = cholesky(Sigma_u).L
    A0 = inv(L)
    A_plus = A0 * coeffs(model)
    return A0, A_plus
end

function identify!(model::SVAR, method::AbstractIdentificationMethod)
    A0, A_plus = _identify(model.var, method)
    model.A0 = A0
    model.A_plus = A_plus
    return model 
end
function identify(model::VAR, method::AbstractIdentificationMethod)
    A0, A_plus = _identify(model, method)
    return SVAR(A_plus, A0, model.p, model.trend_exponents, model)
end

function fit!(model::SVAR, identification_method::AbstractIdentificationMethod)
    fit!(model.var)
    return identify!(model, identification_method)
end

function fit_and_select!(
    model::SVAR, 
    identification_method::AbstractIdentificationMethod,
    ic_function::Function=aic 
)

    var_best, ic_table = fit_and_select!(model.var, ic_function)
    model.var = var_best
    return identify!(model, identification_method), ic_table
end

#-------------------------------------------------------------------------------
# SIMULATION
#-------------------------------------------------------------------------------

function simulate(                                            # k variables, T periods, p lags
    ::Type{SVAR},                                             # 
    shocks::AbstractMatrix{<:Number},                         # k × T
    A0::AbstractMatrix{<:Number},                             # k × k (invertible) 
    A_plus::AbstractMatrix{<:Number};                         # k × kp + m
    trend_exponents::AbstractVector{<:Number}=[0],            # m × 1 
    initial::Union{Nothing,AbstractVector{<:Number}}=nothing  # kp × 1
)

    m = length(trend_exponents)
    kp = size(A_plus, 2) - m
    p = floor(Int, kp / size(A0, 1))

    Phi0 = inv(A0)
    errors = Phi0 * shocks
    B = Phi0 * A_plus

    var = simulate!(VAR, errors, B; trend_exponents=trend_exponents, initial=initial)

    return SVAR(get_input_data(var), p; trend_exponents=trend_exponents)
end

function simulate(
    ::Type{SVAR},
    T::Int,
    A0::AbstractMatrix{<:Number},
    A_plus::AbstractMatrix{<:Number};
    trend_exponents::AbstractVector{<:Number}=[0],
    initial::Union{Nothing,AbstractVector{<:Number}}=nothing
)
    m = length(trend_exponents)
    kp = size(A_plus, 2) - m
    p = floor(Int, kp / size(A0, 1))

    Phi0 = inv(A0)
    B = Phi0 * A_plus
    Sigma_u = Phi0 * Phi0'
    var = simulate(VAR, T, B, Sigma_u; trend_exponents=trend_exponents, initial=initial)

    return SVAR(get_input_data(var), p; trend_exponents=trend_exponents)
end

#-------------------------------------------------------------------------------
# IMPULSE RESPONSE FUNCTIONS
#-------------------------------------------------------------------------------

function IRF(model::SVAR, max_horizon::Int)
    require_fitted(model)
    Phi0 = inv(model.A0)
    irfs = _var_irf(coeffs(model.var), model.p, max_horizon)
    for h=0:max_horizon
        view(irfs, :, :, h+1) .= view(irfs, :, :, h+1) * Phi0 
    end
    varnames = names(get_input_data(model))
    return IRF(irfs, varnames, model)
end

function _identify_irfs(model::VAR, method::Recursive, max_horizon::Int)
    require_fitted(model)
    irfs = _var_irf(coeffs(model), model.p, max_horizon)
    A0, _ = _identify(model, method)
    Phi0 = inv(A0)
    for h = 0:max_horizon
        view(irfs, :, :, h + 1) .= view(irfs, :, :, h + 1) * Phi0
    end
    return irfs
end

function _identify_irfs(model::VAR, method::InternalInstrument, max_horizon::Int)
    require_fitted(model)
    data = get_input_data(model)
    idx_instrument = _find_variable_idx(method.instrument, data)
    idx_normalising = _find_variable_idx(method.normalising_variable, data)

    # internal instrument IRFs are cumputed as relative IRFs of the 
    # Cholesky shock of the instrument variable on the outcome variable 
    # over the Cholesky shock of the instrument variable on the normalising 
    # variable and the pre-defined horizon
    cholesky_irfs = _identify_irfs(model, Recursive(), max_horizon)
    normalising_factor = cholesky_irfs[idx_normalising, idx_instrument, method.normalising_horizon]
    cholesky_irfs ./= normalising_factor

    # we can only interpret the IRFs to the Cholesky shock of the instrument
    K = size(data, 2)
    T = eltype(cholesky_irfs)
    cholesky_irfs[:, filter(!=(idx_instrument), 1:K), :] .= T(NaN)

    return cholesky_irfs
end

function IRF(model::VAR, method::AbstractIdentificationMethod, max_horizon::Int)
    irfs = _identify_irfs(model, method, max_horizon)
    varnames = names(get_input_data(model))
    return IRF(irfs, varnames, model, method)
end


