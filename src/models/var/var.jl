using LinearAlgebra
using DataFrames
import Statistics: cov

"""
    VAR <: Model

Vector Autoregressive (VAR) model in matrix form.

A VAR of lag order `p` is specified as:

```math
    y_t = C e_t + B_1 y_{t-1} + ... + B_p y_{t-p} + u_t
```
where:
- ``e_t`` is a vector of deterministic components (constant, trends, etc)
- ``C, B_i`` are conformable matrices
- ``u_t`` is vector white noise

This can be rewritten compactly as:

```math
    y_t' = z_t' B_+' + u_t'
```

where:
- ``z_t = [e_t; y_{t-1}; ...; y_{t-p}]`` includes deterministic components and
  lagged values of all variables
- ``B_+ = [C, B_1, ..., B_p]`` is the coefficient matrix stacking trend and
  autoregressive terms

Stacking observations from ``t = p+1`` to ``T``, the system can be written in
matrix form:

```math
    Y = X B_+' + U
```

where:
- ``Y``: an ``(T - p) \\times k`` matrix of outcomes
- ``X``: an ``(T - p) \\times (k p + d)`` matrix of regressors (lags and trends)
- ``B_+``: a ``k \\times (kp + d)`` matrix of coefficients where ``d`` is the
  number of deterministic components.
- ``U``: a ``(T - p) \\times k`` matrix of residuals

This structure is represented by the `VAR` object.

# Fields
- `B::Matrix{<:Number}`: Coefficient matrix `[C, B_1, ..., B_p]`
- `Sigma_u::Matrix{<:Number}`: Covariance matrix of the error term
- `p::Int`: Lag order of the VAR
- `trend_exponents::Vector{<:Number}`: Time trend exponents (e.g., `[0,1]`
  implies constant and linear trend)
- `input_data::DataFrame`: Original dataset used to estimate the VAR
- `Y::Matrix{<:Number}`: Left-hand side matrix of outcomes `y_t`, shape `(T - p, k)`
- `X::Matrix{<:Number}`: Right-hand side matrix of regressors `z_t`, shape `(T - p, k * p + d)`
- `U::Matrix{<:Number}`: Residuals `u_t`, shape `(T - p, k)`
- `Yhat::Matrix{<:Number}`: Fitted values `X * B_+'`, shape `(T - p, k)`
"""
mutable struct VAR <: Model

    # defining the DGP
    B::AbstractMatrix{<:Number}                   # Coefficient matrices [C, B1, ..., Bp]
    Sigma_u::AbstractMatrix{<:Number}             # Covariance of error term
    p::Int                                        # Lag order
    trend_exponents::AbstractVector{<:Number}     # Determined trend as exponents for time trend

    # input data
    input_data::DataFrame                         # Data used to fit model

    # intermediate data
    Y::AbstractMatrix{<:Number}                   # LHS matrix
    X::AbstractMatrix{<:Number}                   # RHS matrix
    U::AbstractMatrix{<:Number}                   # Residual Matrix
    Yhat::AbstractMatrix{<:Number}                # Fitted values
end

"""
    VAR(B::Matrix{<:Number},
        Sigma_u::Matrix{<:Number},
        p::Int,
        trend_exponents::Vector{<:Number},
        data::DataFrame)

Constructs a `VAR` model from estimated coefficient matrix `B`, error
covariance matrix `Sigma_u`, lag length `p`, time trend exponents, and
the dataset `data`.

# Arguments
- `B::Matrix{<:Number}`: Coefficient matrix
- `Sigma_u::Matrix{<:Number}`: Covariance matrix of residuals
- `p::Int`: Lag order
- `trend_exponents::Vector{<:Number}`: Trend exponents for time trend terms
- `data::DataFrame`: Dataset used to construct X and Y

"""
function VAR(
    B::AbstractMatrix{<:Number},
    Sigma_u::AbstractMatrix{<:Number},
    p::Int,
    trend_exponents::AbstractVector{<:Number},
    data::DataFrame
)
    type = eltype(data[:, 1])
    U = Matrix{type}(undef, 0, 0)
    Yhat = Matrix{type}(undef, 0, 0)
    Y = Matrix(data[(p+1):end, :])
    # first get lagged values
    X = make_lag_matrix(Matrix(data), p)
    X = X[(p+1):end, :]
    # then add time trends
    for te in reverse(trend_exponents)
        X = hcat(((p+1):size(data, 1)) .^ te, X)
    end

    return VAR(B, Sigma_u, p, trend_exponents, data, Y, X, U, Yhat)
end

"""
    VAR(data::DataFrame, p::Int;
        trend_exponents::Vector{<:Number} = [0])

Constructs a `VAR` model object with data, lag length `p`, and
specified time trend exponents. Coefficients and residuals are uninitialised
but can be estimed using `fit!`.

# Arguments
- `data::DataFrame`: Dataset to construct lag matrix and outcomes
- `p::Int`: Lag length

## Keyword Arguments
- `trend_exponents::Vector{<:Number}`: Exponents of time trends (default: `[0]`)
"""
function VAR(data::DataFrame, p::Int; trend_exponents::AbstractVector{<:Number}=[0])
    type = eltype(data[:, 1])

    # Making missing matrix of zero size
    B = Matrix{type}(undef, 0, 0)
    Sigma_u = Matrix{type}(undef, 0, 0)

    return VAR(B, Sigma_u, p, trend_exponents, data)
end

function coeffs(model::VAR, exclude_deterministic::Bool=false)
    require_fitted(model)
    exclude_deterministic || return model.B

    m = length(model.trend_exponents)
    return model.B[:, (m+1):end]
end
cov(model::VAR) = require_fitted(model) && model.Sigma_u
fitted(model::VAR) = require_fitted(model) && model.Yhat
residuals(model::VAR) = require_fitted(model) && model.U
# effective observations since first p observations are lost
nobs(model::VAR) = size(model.Y, 1)
get_dependent(model::VAR) = model.Y
get_independent(model::VAR) = model.X
get_input_data(model::VAR) = model.input_data
is_fitted(model::VAR) = size(model.Yhat, 1) >= 1
is_structural(model::VAR) = false

function Base.show(io::IO, ::MIME"text/plain", x::VAR)
    varnames = names(get_input_data(x))
    constant = 0 in x.trend_exponents
    trends = ["t^$i" for i in filter(!=(0), x.trend_exponents)]
    s = """
        Vector Autoregression
        =====================
        variables: $(join(varnames, ", "))
        p: $(x.p)
        constant: $(constant)
        trends: $(join(trends, ", "))
        fitted: $(is_fitted(x))
        """
    println(io, s)
end

#-------------------------------------------------------------------------------
# CHECKING MODEL ASSUMPTIONS
#-------------------------------------------------------------------------------

"""
    make_companion_matrix(model::VAR) -> Matrix{<:Number}

Returns the companion matrix of a `VAR(p)` model.

The companion matrix is a square matrix of size `(k * p, k * p)`, where `k`
is the number of variables in the system and `p` is the lag order.

```math
    Z_t = A Z_{t-1} + \\varepsilon_t
```

where ``Z_t`` is a ``k \\times p``-dimensional stacked vector of lagged variables.

The matrix has the following block form:

```math
A = \\begin{bmatrix}
    B_1 & B_2 & ... & B_p \\\\
    I_k & 0   & ... & 0   \\\\
    0   & I_k & ... & 0   \\\\
    ... & ... & ... & ...
\\end{bmatrix}
```

where ``B_1, ..., B_p`` are the VAR coefficient matrices, and ``I_k`` is the
identity matrix of size `k`.

# Arguments
- `model::VAR`: A fitted VAR model.
"""
function make_companion_matrix(model::VAR)
    require_fitted(model)
    return make_companion_matrix(coeffs(model), model.p, length(model.trend_exponents))
end

"""
    spectral_radius(model::Union{VAR,SVAR}) -> Real

Computes the spectral radius of the companion matrix of a fitted `VAR` model.

The spectral radius of a matrix ``A`` is defined as the largest absolute value
among its eigenvalues:

```math
\\rho(A) = \\max_i |\\lambda_i|
```

where ``\\lambda_i`` are the eigenvalues of ``A``.

In the context of VAR models, the spectral radius of the companion matrix is
a key measure of stability. A `VAR` model is stable if the
spectral radius is strictly less than 1.

# Arguments
- `model::VAR`: A fitted VAR model

"""
function spectral_radius(model::VAR)
    require_fitted(model)
    return spectral_radius(make_companion_matrix(model))
end

"""
    is_stable(model::Union{VAR,SVAR}) -> Bool
    is_stable(C::AbstractMatrix) -> Bool

Checks whether a VAR model or a companion matrix corresponds to a stable
VAR process.

A VAR system is stable if all eigenvalues of its companion matrix lie strictly
inside the unit circle. This is equivalent to requiring that the spectral
radius is less than one:

```math
\\rho(C) = \\max_i |\\lambda_i| < 1
```
# Arguments
- `model::VAR`: A fitted VAR model
- `C::AbstractMatrix`: A companion matrix

# Returns
- `Bool`: `true` if the VAR is stable, `false` otherwise
"""
is_stable(C::AbstractMatrix) = (spectral_radius(C) < 1)
is_stable(model::VAR) = require_fitted(model) && (spectral_radius(model) < 1)

#-------------------------------------------------------------------------------
# INFORMATION CRITERIA
#
# Beware that when comparing ICs over various lag orders, we must make
# sure that they are compared over the same dataset. That is, estimations
# with lower p must adjust their Sigma_u estimation to those data points
# also used by the highest p. That's the reason for the two types of functions
# below.
#
# ICs for VARs are discussed in Section 2.6. of Kilia & Lütkepohl (2017).
#
# Kilian, L., & Lütkepohl, H. (2017).
# Structural Vector Autoregressive Analysis: (1st ed.).
# Cambridge University Press. https://doi.org/10.1017/9781108164818
#-------------------------------------------------------------------------------

function _ic(Sigma_u::AbstractMatrix{<:Number}, num_coeffs::Int, ct::Number)
    return logdet(Sigma_u) + ct * num_coeffs
end
function _ic(model::VAR, ct::Number)
    require_fitted(model)
    Sigma_u = model.Sigma_u
    K = size(model.Y, 2)
    num_coeffs = K * size(model.X, 1)
    return _ic(Sigma_u, num_coeffs, ct)
end

function aic(Sigma_u::AbstractMatrix{<:Number}, num_coeffs::Int, T::Int)
    ct = 2 / T
    return _ic(Sigma_u, num_coeffs, ct)
end
function aic(model::VAR)
    T = nobs(model)
    ct = 2 / T
    return _ic(model, ct)
end

function hqc(Sigma_u::AbstractMatrix{<:Number}, num_coeffs::Int, T::Int)
    ct = 2 * log(log(T)) / T
    return _ic(Sigma_u, num_coeffs, ct)
end
function hqc(model::VAR)
    T = nobs(model)
    ct = 2 * log(log(T)) / T
    return _ic(model, ct)
end

function sic(Sigma_u::AbstractMatrix{<:Number}, num_coeffs::Int, T::Int)
    ct = log(T) / T
    return _ic(Sigma_u, num_coeffs, ct)
end
function sic(model::VAR)
    T = nobs(model)
    ct = log(T) / T
    return _ic(model, ct)
end

bic(Sigma_u::AbstractMatrix{<:Number}, num_coeffs::Int, T::Int) = sic(Sigma_u, num_coeffs, T)
bic(model::VAR) = sic(model)

"""
    aic(model::Union{VAR,SVAR}) -> Real
    aic(Sigma_u::Matrix{<:Number}, num_coeffs::Int, T::Int) -> Real

    bic(model::Union{VAR,SVAR}) -> Real
    bic(Sigma_u::Matrix{<:Number}, num_coeffs::Int, T::Int) -> Real

    sic(model::Union{VAR,SVAR}) -> Real
    sic(Sigma_u::Matrix{<:Number}, num_coeffs::Int, T::Int) -> Real

    hqc(model::Union{VAR,SVAR}) -> Real
    hqc(Sigma_u::Matrix{<:Number}, num_coeffs::Int, T::Int) -> Real

Computes information criteria for model selection in VARs.

Given a fitted `VAR` model or directly the covariance matrix `Sigma_u`, the
number of estimated coefficients, and sample size `T`, each function returns
the respective criterion value:

- **AIC** (Akaike Information Criterion):

```math
\\text{AIC} = \\log\\det(\\Sigma_u) + \\frac{2k}{T}
```

- **SIC/BIC** (Schwarz/Bayesian Information Criterion):

```math
\\text{SIC} = \\log\\det(\\Sigma_u) + \\frac{\\log(T) k}{T}
```

- **HQC** (Hannan-Quinn Criterion):

```math
\\text{HQC} = \\log\\det(\\Sigma_u) + \\frac{2 \\log(\\log T) k}{T}
```

Here, ``k`` is the number of estimated parameters and ``T`` is the effective
sample size.

# Arguments
- `model::VAR`: A fitted VAR model
- `Sigma_u::Matrix{<:Number}`: Residual covariance matrix
- `num_coeffs::Int`: Number of estimated coefficients
- `T::Int`: Number of observations used in estimation
"""
aic, sic, bic, hqc

#-------------------------------------------------------------------------------
# ESTIMATION FUNCTIONS
#-------------------------------------------------------------------------------

"""
    fit!(model::VAR)

Estimate a `VAR` model using OLS.
"""
function fit!(model::VAR)
    X, Y = model.X, model.Y
    model.B = (Y' * X) / (X' * X)
    model.Yhat = model.X * model.B'
    model.U = model.Y - model.Yhat
    model.Sigma_u = model.U' * model.U / nobs(model)
    return model
end

"""
    fit_and_select!(model::VAR, ic_function::Function=aic) --> (VAR, DataFrame)

Select and estimate a `VAR` model.

The best model is determined by the model with the smallest `ic_function` value
among all models with `p=1:model.p`. Thus, the lag-length of the provided model
determines the maximum lag length.

Available choices for `ic_function` are `aic`, `bic`, `sic`, `hqc`, but user
defined functions can be provided as long as they have the signature
`ic_function(Sigma_u::Matrix{<:Number}, num_coeffs::Int, T::Int)`
where `Sigma_u` is the VAR error covariance matrix, `num_coeffs` is the number
of estimated coefficients, and `T` is a number of effective observations.

To be correct, the error covariance matrix of all models is estimated over the
same time period. Calling `aic` or other functions on manually estimated models \
with differing lag-lengths will not compare the models on the same time
period -- the model with higher `p` will have fewer effective number of
observations`. It is thus recommended to do model comparison via this function.

# Arguments
- `model::VAR`: VAR model, where the provided lag-length `p` determines the
  maximum lag-length.
- `ic_function::Function`: Information criterion function. See the details above.
  Default is AIC, since it is generally recommended to go with more rather than
  fewer lags.

# Returns
Returns a tuple `(VAR, DataFrame)` where the first element is the best model and
the second element is a table with information regarding the `ic_function` value
for each estimated model. Note that manually calling `aic` or similar functions
on the returned model might not provide that same value, since the covariance
will now be estimated over the full period rather than the common period.
"""
function fit_and_select!(model::VAR, ic_function::Function=aic)
    # model.p provides the maximum lag length in this case
    # ic_function(Sigma_u, num_coeffs, T)
    # returns model_best, ICs table

    p_max = model.p
    ps = 0:p_max
    ics = zeros(length(ps))

    # Given model gives the largest lag length
    # We will go from largest to smallest.
    fit!(model)
    Sigma_u = model.Sigma_u
    n_coeffs = length(coeffs(model))
    T = nobs(model)
    ic_best = ic_function(Sigma_u, n_coeffs, T)
    ics[p_max+1] = ic_best
    model_best = model

    for p = (p_max-1):-1:0
        model_tmp = VAR(get_input_data(model), p; trend_exponents=model.trend_exponents)
        fit!(model_tmp)
        n_coeffs_tmp = length(coeffs(model_tmp))
        U_tmp = residuals(model_tmp)
        U_tmp = U_tmp[(p_max-p+1):end, :]  # model with pmax had fewer observations
        T == size(U_tmp, 1) || error("U adjustment is wrong.")  # TODO: remove?
        Sigma_u_tmp = U_tmp' * U_tmp / T

        ic_tmp = ic_function(Sigma_u_tmp, n_coeffs_tmp, T)
        ics[p+1] = ic_tmp
        if ic_tmp < ic_best
            model_best = model_tmp
            ic_best = ic_tmp
        end
    end

    return model_best, DataFrame(p=ps, IC=ics)
end

#-------------------------------------------------------------------------------
# SIMULATION
#-------------------------------------------------------------------------------

"""
    function _simulate_var!(                                   # k variables, T periods
        errors::AbstractMatrix{M},                         # k × T
        B::AbstractMatrix{M};                              # k × kp+m
        trend_exponents::AbstractVector{<:Real}=[0],       # m × 1
        initial::Union{Nothing,AbstractVector{M}}=nothing  # kp × 1
    ) where {M<:Real}

Simulate a VAR(p) overwriting `errors` with the simulated data.
"""

"""
    _simulate_var!(errors::Matrix{<:Number},
              B::Matrix{<:Number};
              trend_exponents::Vector{<:Number} = [0],
              initial::Union{Nothing, Vector{<:Number}} = nothing) -> Matrix{<:Number}

Simulate a VAR(p) overwriting `errors` with the simulated data.

The simulated process follows the structure:

```math
    y_t = C e_t + B_1 y_{t-1} + ... + B_p y_{t-p} + u_t
```

where the full coefficient matrix is `B = [C B_1 ... B_p]`. The simulation
accounts for deterministic trends via `trend_exponents`. For example,
`trend_exponents=[0,1]` implies that a constant and a linear trend are included.

# Arguments
- `errors::Matrix{<:Number}`: A `(k × T)` matrix that is overwritten with the
  simulated data. Initially contains the error terms `u_t`. `k` is the number
  of endogneous variables.
- `B::Matrix{<:Number}`: The full coefficient matrix of size `(k × (k * p + m))`,
  where `m` is the number of deterministic trend terms.
- `trend_exponents::Vector{<:Number}`: Exponents of time to model deterministic
  components (e.g., `[0, 1]` gives constant and linear trend). Default is `[0]`,
  i.e. a constant.
- `initial::Union{Nothing, Vector{<:Number}}`: Initial values for lagged variables,
  a vector of length `k * p`. If `nothing`, lags are initialised at zero.

# Returns
- The matrix `errors`, now containing the simulated VAR data.

"""
function _simulate_var!(                                             # k variables, T periods
    errors::AbstractMatrix{<:Number},                            # k × T
    B::AbstractMatrix{<:Number};                                 # k × kp+m
    trend_exponents::AbstractVector{<:Number}=[0],               # m × 1
    initial::Union{Nothing,AbstractVector{<:Number}}=nothing     # kp × 1
)

    M = eltype(errors)
    k, T = size(errors)
    m = length(trend_exponents)
    kp = size(B, 2) - m
    # needs this check because otherwise we can have wrong dimensions without
    # getting an error
    kp % k == 0 || error("Dimensions of B are wrong.")
    if isnothing(initial)
        initial = zeros(M, kp)
    end
    Zt = ones(M, m + kp)
    view(Zt, (m+1):(m+kp)) .= view(initial, :)
    for t = 1:T
        _make_trend!(view(Zt, 1:m), t, trend_exponents)
        mul!(view(errors, :, t), B, Zt, one(M), one(M))
        if kp > 0
            # no need to rotate in if we do not have lagged terms
            _rotate_in!(view(Zt, (m+1):(m+kp)), view(errors, :, t))
        end
    end

    return errors
end

# overwrites errors
function simulate!(                                             # k variables, T periods, p lags
    ::Type{VAR},                                                  #
    errors::AbstractMatrix{<:Number},                           # k  × T
    B::AbstractMatrix{<:Number};                                # k  × kp + m
    trend_exponents::AbstractVector{<:Number}=[0],              # m  × 1
    initial::Union{Nothing,AbstractVector{<:Number}}=nothing   # kp × 1
)

    errors = _simulate_var!(errors, B; trend_exponents=trend_exponents, initial=initial)
    data = DataFrame(errors', "Y" .* string.(1:size(errors, 1)))

    m = length(trend_exponents)
    k = size(B, 1)
    kp = size(B, 2) - m
    p = floor(Int, kp / k)

    return VAR(data, p; trend_exponents=trend_exponents)
end
function simulate(
    ::Type{VAR},
    T::Int,
    B::AbstractMatrix{<:Number},
    Sigma_u::AbstractMatrix{<:Number}=I(size(B, 1));
    trend_exponents::AbstractVector{<:Number}=[0],
    initial::Union{Nothing,AbstractVector{<:Number}}=nothing
)
    k = size(B, 1)
    errors = cholesky(Sigma_u).L * randn(k, T)
    return simulate!(VAR, errors, B; trend_exponents=trend_exponents, initial=initial)
end

"""
    simulate(::Type{VAR}, T::Int, B::Matrix{<:Number},
             Sigma_u::Matrix{<:Number}=I(size(B,1));
             trend_exponents::Vector{<:Number}=[0],
             initial::Union{Nothing,Vector{<:Number}}=nothing) -> VAR

    simulate!(::Type{VAR}, errors::Matrix{<:Number}, B::Matrix{<:Number};
              trend_exponents::Vector{<:Number}=[0],
              initial::Union{Nothing,Vector{<:Number}}=nothing) -> VAR

Simulates a `VAR(p)` model using the specified coefficient matrix `B` and
optionally a covariance matrix `Sigma_u` and deterministic trends.

The simulation uses the reduced-form VAR representation:

```math
    y_t = C e_t + B_1 y_{t-1} + ... + B_p y_{t-p} + u_t
```

# Method 1: `simulate!`
Simulates a `VAR` process by overwriting the provided error matrix `errors`
with simulated values. Returns a `VAR` object constructed from the simulated
data.

# Method 2: `simulate`
Generates error terms internally from a Gaussian distribution with covariance
`Sigma_u` (default is the identity matrix). Returns a `VAR` object containing the
simulated series.

# Arguments
- `B::Matrix{<:Number}`: Coefficient matrix `[C B_1 ... B_p]`, size `k × (k * p + m)`,
  where `m` is the number of deterministic components.
- `trend_exponents::Vector{<:Number}`: Exponents used to simulate trends
  (e.g. `[0,1]` implies constant and linear trend)
- `initial::Union{Nothing, Vector{<:Number}}`: Optional initial values for
  lags (length `k * p`). Default is zero.
- `T::Int`: Number of time periods to simulate
- `Sigma_u::Matrix{<:Number}`: Covariance matrix of the error term
  (default is identity)
- `errors::Matrix{<:Number}`: Pre-allocated matrix of shape `k × T`, initially
  filled with innovations, overwritten with simulated data

# Returns
- `VAR`: A new `VAR` object containing the simulated dataset. The data can
  be obtained using `get_input_data`. Alternatively, the model can be directly
  estimated using `fit!`.
"""
simulate, simulate!

#-------------------------------------------------------------------------------
# IMPULSE RESPONSE FUNCTIONS
#-------------------------------------------------------------------------------

"""
    _var_irf(B::Matrix{<:Number}, p::Int, max_horizon::Int) -> Array{<:Number, 3}

Internal function to compute impulse response functions (IRFs) up to
horizon `max_horizon` from a reduced-form VAR(p) model.

!!! warning
    The matrix `B` must **exclude** coefficients for deterministic components
    (constant, trend, etc.). Use `coeffs(model, true)` when extracting
    coefficients from a `model::VAR` object to ensure the correct structure.

The resulting IRFs have the shape `(k, k, H + 1)`, where:
- `k`: number of variables in the system
- `H`: maximum impulse response horizon
- The first dimension is the responding variable
- The second dimension is the shock
- The third dimension indexes the horizon

# Arguments
- `B::Matrix{<:Number}`: VAR coefficient matrix of shape `(k, k * p)`,
  excluding deterministic components
- `p::Int`: Lag order of the VAR
- `max_horizon::Int`: Maximum forecast horizon for the IRFs

# Returns
- `Array{<:Number, 3}`: IRF array of size `(k, k, max_horizon + 1)`
"""
function _var_irf(B::AbstractMatrix{<:Number}, p::Int, max_horizon::Int)
    k = size(B, 1)
    T = eltype(B)
    irfs = zeros(T, k, k, max_horizon + 1)
    copyto!(view(irfs, :, :, 1), diagm(ones(k)))

    for h = 1:max_horizon
        for j = 1:min(h, p)
            Bj = view(B, :, ((j-1)*k+1):(j*k))
            mul!(view(irfs, :, :, h + 1), Bj, view(irfs, :, :, h + 1 - j), 1.0, 1.0)
        end
    end

    return irfs
end

function IRF(model::VAR, max_horizon::Int)
    irfs = _var_irf(coeffs(model, true), model.p, max_horizon)
    varnames = Symbol.(names(get_input_data(model)))
    return IRF(irfs, varnames, model)
end
