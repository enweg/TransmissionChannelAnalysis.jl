"""
    SVAR <: Model

Structural Vector Autoregressive (SVAR) model in matrix form.

A SVAR of lag order `p` is specified as:

```math
    A_0y_t = C e_t + A_1 y_{t-1} + ... + A_p y_{t-p} + \\varepsilon_t
```
where:
- ``e_t`` is a vector of deterministic components (constant, trends, etc)
- ``C, A_i`` are conformable matrices
- ``\\varepsilon_t`` are structural shocks

This can be rewritten compactly as:

```math
    y_t'A0' = z_t' A_+' + u_t'
```

where:
- ``z_t = [e_t; y_{t-1}; ...; y_{t-p}]`` includes deterministic components and
  lagged values of all variables
- ``A_+ = [C, A_1, ..., A_p]`` is the coefficient matrix stacking trend and
  autoregressive terms

Assuming ``A0`` is invertible, the reduced-form VAR can be obtained as 

```math
y_t' = z_t' A_+'(A_0')^{-1} + u_t'(A_0')^{-1}
```

which can be represented using a `VAR` object.

# Fields
- `A_plus::Matrix{<:Number}`: Coefficient matrix `[C, A_1, ..., A_p]`
- `A_0::AbstractMatrix{<:Number}`: Contemporaneous relationships.
- `p::Int`: Lag order of the (S)VAR
- `trend_exponents::Vector{<:Number}`: Time trend exponents (e.g., `[0,1]`
  implies constant and linear trend)
- `var::VAR`: Reduced form representation of SVAR 
"""
mutable struct SVAR <: Model
    # defining the DGP
    A_plus::AbstractMatrix{<:Number}              # Coefficient matrices [C, A1, ..., Ap]
    A0::AbstractMatrix{<:Number}                  # Contemporaneous relationships
    p::Int                                        # Lag order
    trend_exponents::AbstractVector{<:Number}     # Determinstic trend as exponents for time trend

    # reduced form model
    var::VAR
end

"""
    SVAR(data::DataFrame, 
         p::Int;
         trend_exponents::Vector{<:Number} = [0])

Constructs a `SVAR` model object with data, lag length `p`, and
specified time trend exponents. Coefficients and residuals are uninitialised 
but can be estimed using `fit!` and an appropriate `AbstractIdentificationMethod`.

# Arguments
- `data::DataFrame`: Dataset used for the estimation
- `p::Int`: Lag length

## Keyword Arguments
- `trend_exponents::Vector{<:Number}`: Exponents of time trends (default: `[0]`)
"""
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

"""
    shocks(model::SVAR) --> AbstractMatrix{<:Number}

Returns the identified shocks.
"""
shocks(model::SVAR) = require_fitted(model) && residuals(model.var) * model.A0'
nobs(model::SVAR) = nobs(model.var)
get_dependent(model::SVAR) = get_dependent(model.var)
get_independent(model::SVAR) = get_independent(model.var)
get_input_data(model::SVAR) = get_input_data(model.var)
is_fitted(model::SVAR) = size(model.A0, 1) >= 1
is_structural(model::SVAR) = true

"""
    make_companion_matrix(model::SVAR) --> Matrix{<:Number}

Returns the companion matrix for the `VAR` representation of the `SVAR`. 
"""
make_companion_matrix(model::SVAR) = require_fitted(model) && make_companion_matrix(model.var)
spectral_radius(model::SVAR) = require_fitted(model) && spectral_radius(model.var)
is_stable(model::SVAR) = require_fitted(model) && is_stable(model.var)

aic(model::SVAR) = require_fitted(model) && aic(model.var)
hqc(model::SVAR) = require_fitted(model) && hqc(model.var)
sic(model::SVAR) = require_fitted(model) && sic(model.var)
bic(model::SVAR) = require_fitted(model) && bic(model.var)

function Base.show(io::IO, ::MIME"text/plain", x::SVAR)
    varnames = names(get_input_data(x))
    constant = 0 in x.trend_exponents
    trends = ["t^$i" for i in filter(!=(0), x.trend_exponents)]
    s = """
        Structural Vector Autoregression
        ================================
        variables: $(join(varnames, ", "))
        p: $(x.p)
        constant: $(constant)
        trends: $(join(trends, ", "))
        fitted: $(is_fitted(x))
        """
    println(io, s)
end
#-------------------------------------------------------------------------------
# ESTIMATION AND IDENTICATION FUNCTIONS
#-------------------------------------------------------------------------------

"""
    _identify(model::SVAR, 
              method::AbstractIdentificationMethod
             ) --> (Matrix{<:Number}, Matrix{<:Number})

    _identify(B::AbstractMatrix{<:Number}, 
              Sigma_u::AbstractMatrix{<:Number}, 
              method::AbstractIdentificationMethod
             ) --> (Matrix{<:Number}, Matrix{<:Number})

Internal method for identifying structural matrices in SVAR models.

Must return `A0` and `A_plus` matrices. Should not be implemented by
models that do not use structural matrices (e.g., local projections).

For more details, see the `SVAR` documentation.

## Arguments
- `model::VAR`: An estimated VAR model
- `method::AbstractIdentificationMethod`: An identification method to identify 
  the SVAR from the VAR. 
- `B::AbstractMatrix{<:Number}`: Coefficient matrix [C B_1 ... B_p] of the VAR. 
  See the `VAR` documention for more information. 
- `Sigma_u::AbstractMatrix{<:Number}`: Covariance matrix of VAR residuals. See
  the `VAR` documentation for more information. 
"""
function _identify(::VAR, ::I) where {I<:AbstractIdentificationMethod}
    error("_identify is not implemented for SVAR and method $I.")
end

function _identify(
    B::AbstractMatrix{<:Number}, 
    Sigma_u::AbstractMatrix{<:Number}, 
    ::Recursive
)

    L = cholesky(Sigma_u).L
    A0 = inv(L)
    A_plus = A0 * B
    return A0, A_plus
end
function _identify(model::VAR, ::Recursive)
    is_fitted(model) || error("Reduced form model must first be estimated.")
    Sigma_u = cov(model)
    B = coeffs(model)
    return _identify(B, Sigma_u, Recursive())
end

function _identify(
    ::AbstractMatrix{<:Number}, 
    ::AbstractMatrix{<:Number}, 
    ::InternalInstrument
)
    error("Internal instruments can only be used to identify IRFs but not the full SVAR.")
end
function _identify(model::VAR, ::InternalInstrument)
    require_fitted(model)
    error("Internal instruments can only be used to identify IRFs but not the full SVAR.")
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

"""
    fit!(model::SVAR, identification_method::AbstractIdentificationMethod) --> SVAR

Estimate an SVAR using `identification_method`. 

## Argument
- `model::SVAR`: A SVAR model to be estimated
- `identification_method::AbstractIdentificationMethod`: Identification method 
  used to identify SVAR from reduced-form VAR
"""
function fit!(model::SVAR, identification_method::AbstractIdentificationMethod)
    fit!(model.var)
    return identify!(model, identification_method)
end

"""
    fit_and_select!(model::SVAR, 
                    identification_method::AbstractIdentificationMethod, 
                    ic_function::Function=aic) --> (SVAR, DataFrame)

Select and estimate a `SVAR` model by first selecting an estimating a `VAR` 
model and then identifying the `SVAR` from the `VAR` using `identification_method`. 
"""
function fit_and_select!(
    model::SVAR, 
    identification_method::AbstractIdentificationMethod,
    ic_function::Function=aic 
)

    var_best, ic_table = fit_and_select!(model.var, ic_function)
    model.var = var_best
    model.p = var_best.p
    return identify!(model, identification_method), ic_table
end

#-------------------------------------------------------------------------------
# SIMULATION
#-------------------------------------------------------------------------------

"""
    simulate(::Type{SVAR}, shocks::Matrix{<:Number}, A0::Matrix{<:Number},
             A_plus::Matrix{<:Number};
             trend_exponents::Vector{<:Number}=[0],
             initial::Union{Nothing,Vector{<:Number}}=nothing) -> SVAR

    simulate(::Type{SVAR}, T::Int, A0::Matrix{<:Number}, A_plus::Matrix{<:Number};
             trend_exponents::Vector{<:Number}=[0],
             initial::Union{Nothing,Vector{<:Number}}=nothing) -> SVAR

Simulates a structural VAR (SVAR) model using the structural form
characterised by matrices `A0` and `A_plus`.

The SVAR model is given by:

```math
    A_0 y_t = C e_t + A_1 y_{t-1} + ... + A_p y_{t-p} + \\varepsilon_t
```

# Method 1: `simulate(::Type{SVAR}, shocks, A0, A_plus)`
Simulates an SVAR process given structural shocks. 

# Method 2: `simulate(::Type{SVAR}, T, A0, A_plus)`
Generates random structural shocks from a multivariate Guassian with identify 
covariance and simulates the SVAR process using these shocks.

# Arguments
- `A0::Matrix{<:Number}`: Contemporaneous impact matrix (invertible)
- `A_plus::Matrix{<:Number}`: Coefficient matrix `[C A_1 ... A_p]` 
   (size `k × (k*p + m)`) where `m` is the number of exogenous components. 
- `shocks::Matrix{<:Number}`: Structural shocks of shape `k × T`
- `T::Int`: Number of time periods to simulate
- `trend_exponents::Vector{<:Number}`: Exponents for deterministic trends
  (e.g., `[0,1]` for constant and linear trend)
- `initial::Union{Nothing, Vector{<:Number}}`: Optional initial conditions
  for lagged variables (length `k * p`). If `nothing` lags will be initialised 
  at zero.

# Returns
- `SVAR`: A simulated `SVAR` model containing the generated dataset.
  The data can be accessed using `get_input_data`, or the model can be
  estimated using `fit!`.
"""
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

"""
    _svar_irf(A0::Matrix{<:Number}, A_plus::Matrix{<:Number},
              p::Int, max_horizon::Int) -> Array{<:Number, 3}

Internal function to compute structural impulse response functions (IRFs) up to
horizon `max_horizon` from a structural VAR(p) model.

!!! warning
    The matrix `A_plus` must exclude any deterministic components (e.g., constant,
    trend). Use `coeffs(model, true)` to extract it properly from a `SVAR` model.

The resulting IRFs have the shape `(k, k, H + 1)`, where:
- `k`: number of variables
- `H`: maximum impulse response horizon
- First dimension: responding variable
- Second dimension: shock
- Third dimension: horizon

# Arguments
- `A0::Matrix{<:Number}`: Contemporaneous impact matrix (must be invertible)
- `A_plus::Matrix{<:Number}`: Stacked lag coefficient matrix 
  `[A_1 A_2 ... A_p]`, excluding deterministic terms
- `p::Int`: Lag order
- `max_horizon::Int`: Maximum horizon for IRFs

# Returns
- `Array{<:Number, 3}`: Structural IRF array of size `(k, k, max_horizon + 1)`
"""
function _svar_irf(
    A0::AbstractMatrix{<:Number}, 
    A_plus::AbstractMatrix{<:Number},   # excludes any exogenous terms
    p::Int, 
    max_horizon::Int
)

    Phi0 = inv(A0)
    B_plus = Phi0 * A_plus
    irfs_var = _var_irf(B_plus, p, max_horizon)
    irfs = mapslices(x -> x * Phi0, irfs_var; dims = [1, 2])
    return irfs
end

function IRF(model::SVAR, max_horizon::Int)
    require_fitted(model)
    A0, A_plus = coeffs(model, true)
    irfs = _svar_irf(A0, A_plus, model.p, max_horizon)
    varnames = Symbol.(names(get_input_data(model)))
    return IRF(irfs, varnames, model)
end

"""
    _identify_irfs(model::VAR, method::I, max_horizon::Int)
        where {I<:AbstractIdentificationMethod}

Internal method to directly identify IRFs from a `VAR` model using the 
identification method `method`.

Must return a 3-dimensional array with dimensions (variables, shocks,
horizons) where the maximum horizon is given by `max_horizon`. 
"""
function _identify_irfs(::VAR, ::I, max_horizon::Int) where {I<:AbstractIdentificationMethod}
    error("_identify_irfs is not implemented for model VAR and method $I.")
end

"""
    _identify_irfs(B::AbstractMatrix{<:Number}, 
                   Sigma_u::AbstractMatrix{<:Number}, 
                   p::Int, 
                   method::AbstractIdentificationMethod, 
                   max_horizon::Int) --> AbstractArray{<:Number, 3}

Internal method to directly identify IRFs from `VAR` matrices `B` and `Sigma_u`
using identification method `method`. 

## Arguments
- `B::AbstractMatrix{<:Number}`: Stacked matrix of lag matrices [B_1, ..., B_p]
  excluding coefficients on deterministic components. 
- `Sigma_u::AbstractMatrix{<:Number}`: Covariance matrix of `VAR` residuals. 
- `p::Int`: Lag-length of `VAR`. 
- `method::AbstractIdentificationMethod`: Identification method used to identify 
  `SVAR` IRFs from `VAR`. 
- `max_horizon::Int`: Maximum horizon for impulse responses. 
"""
function _identify_irfs(
    B::AbstractMatrix{<:Number},        # excludes deterministic components
    Sigma_u::AbstractMatrix{<:Number}, 
    p::Int,
    ::Recursive, 
    max_horizon::Int
)

    A0, A_plus = _identify(B, Sigma_u, Recursive())
    return _svar_irf(A0, A_plus, p, max_horizon)
end
function _identify_irfs(model::VAR, method::Recursive, max_horizon::Int)
    require_fitted(model)

    B = coeffs(model, true)
    Sigma_u = cov(model)
    return _identify_irfs(B, Sigma_u, model.p, method, max_horizon)
end

function _identify_irfs(
    B::AbstractMatrix{<:Number},         # excludes deterministic components
    Sigma_u::AbstractMatrix{<:Number}, 
    p::Int, 
    method::InternalInstrument, 
    max_horizon::Int
)

    isa(method.instrument, Symbol) && error("Pass the model instead of the coefficient matrices if you want to specify the instrument by name.")
    isa(method.normalising_variable, Symbol) && error("Pass the model instead of the coefficient matrices if you want to specify the normalising variable by name.")
    idx_instrument = method.instrument
    idx_normalising = method.normalising_variable

    # internal instrument IRFs are cumputed as relative IRFs of the 
    # Cholesky shock of the instrument variable on the outcome variable 
    # over the Cholesky shock of the instrument variable on the normalising 
    # variable and the pre-defined horizon
    cholesky_irfs = _identify_irfs(B, Sigma_u, p, Recursive(), max_horizon)
    normalising_factor = cholesky_irfs[idx_normalising, idx_instrument, method.normalising_horizon + 1]
    cholesky_irfs ./= normalising_factor

    K = size(cholesky_irfs, 1)
    T = eltype(cholesky_irfs)
    cholesky_irfs[:, filter(!=(idx_instrument), 1:K), :] .= T(NaN)
    return cholesky_irfs
end
function _identify_irfs(model::VAR, method::InternalInstrument, max_horizon::Int)
    require_fitted(model)
    data = get_input_data(model)
    idx_instrument = _find_variable_idx(method.instrument, data)
    idx_normalising = _find_variable_idx(method.normalising_variable, data)

    method_tmp = InternalInstrument(idx_instrument, idx_normalising, method.normalising_horizon)

    B = coeffs(model, true)
    Sigma_u = cov(model)

    return _identify_irfs(B, Sigma_u, model.p, method_tmp, max_horizon)
end

function IRF(model::VAR, method::AbstractIdentificationMethod, max_horizon::Int)
    irfs = _identify_irfs(model, method, max_horizon)
    varnames = Symbol.(names(get_input_data(model)))
    return IRF(irfs, varnames, model, method)
end


