using LinearAlgebra
using DataFrames


"""
# Fields
- `B::AbstractMatrix{<:Number}`: Coefficient matrix of the form 
   ``[C, B_1, \\ldots, B_p]``, where ``C`` holds the coefficients for the 
   deterministic components and ``B_i`` is the ``i``th lag matrix. 
- `Sigma_u::AbstractMatrix{<:Number}`: Covariance matrix of the error term. 
- `p::Int`: Lag order of the VAR. 
- `trend_exponents::AbstractVector{<:Number}`: Determines the deterministic trends 
   using exponents of time. E.g. `trend_exponents=[0, 1]` implies a constant and
   a linear trend. 
- `input_data::DataFrame`: Data used to estimate the VAR. 
- `Ylag::AbstractMatrix{<:Number}`

```math
\\begin{split}
y_t &= C e_t + B_1 y_{t-1} + ... + B_p y_{t-p} + u_t \\\\
y_t &= [C, B_1, ..., B_p][e_t; y_{t-1}; ...; y_{t-p}] + u_t
\\end{split}
```

```math
\\begin{split}
B_+ &= [C, B_1, ..., B_p] \\ 
z_t &= [e_t; y_{t-1}; ...; y_{t-p}]
\\end{split}
```

```math
y_t' = z_t'B_+' + u_t'
```

```math
\\begin{split}
Y &= [y_{p+1}'; y_{p+2}'; ...; y_{T-1}'; y_{T}'] \\\\
U &= [u_{p+1}'; u_{p+2}'; ...; u_{T-1}'; u_{T}'] \\\\
X &= [z_{p+1}'; z_{p+2}'; ...; z_{T-1}'; z_{T}']
\\end{split}
```

```math
Y = XB_+' + U
```
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
        X = hcat(((p+1):size(data, 1)).^te, X) 
    end

    return VAR(B, Sigma_u, p, trend_exponents, data, Y, X, U, Yhat)
end
function VAR(data::DataFrame, p::Int; trend_exponents::AbstractVector{<:Number}=[0])
    type = eltype(data[:, 1])

    # Making missing matrix of zero size
    B = Matrix{type}(undef, 0, 0)
    Sigma_u = Matrix{type}(undef, 0, 0)

    return VAR(B, Sigma_u, p, trend_exponents, data)
end

coeffs(model::VAR) = model.B
cov(model::VAR) = model.Sigma_u
fitted(model::VAR) = model.Yhat
residuals(model::VAR) = model.U
# effective observations since first p observations are lost
nobs(model::VAR) = size(model.Y, 1)
get_dependent(model::VAR) = model.Y
get_independent(model::VAR) = model.X
get_input_data(model::VAR) = model.input_data
is_fitted(model::VAR) = size(model.Yhat, 1) >= 1

#-------------------------------------------------------------------------------
# CHECKING MODEL ASSUMPTIONS
#-------------------------------------------------------------------------------

make_companion_matrix(model::VAR) = make_companion_matrix(coeffs(model), model.p, length(model.trend_exponents))
spectral_radius(model::VAR) = spectral_radius(make_companion_matrix(model))
is_stable(C::AbstractMatrix) = (spectral_radius(C) < 1)
is_stable(model::VAR) = (spectral_radius(model) < 1)

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

#-------------------------------------------------------------------------------
# ESTIMATION FUNCTIONS
#-------------------------------------------------------------------------------

function fit!(model::VAR)
    X, Y = model.X, model.Y
    model.B = (Y' * X) / (X' * X)
    model.Yhat = model.X * model.B'
    model.U = model.Y - model.Yhat
    model.Sigma_u = model.U' * model.U / nobs(model)
    return model
end

function fit_and_select!(model::VAR, ic_function::Function)
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
    ics[p_max + 1] = ic_best
    model_best = model

    for p = (p_max-1):-1:0
        model_tmp = VAR(get_input_data(model), p; trend_exponents=model.trend_exponents)
        fit!(model_tmp)
        n_coeffs_tmp = length(coeffs(model_tmp))
        U_tmp = residuals(model_tmp)
        U_tmp = U_tmp[(p_max-p+1):end, :]  # model with pmax had fewer observations
        T == size(U_tmp, 1) || error("U adjustment is wrong.")  # TODO: remove
        Sigma_u_tmp = U_tmp' * U_tmp / T
        
        ic_tmp = ic_function(Sigma_u_tmp, n_coeffs_tmp, T)
        ics[p+1] = ic_tmp
        if ic_tmp < ic_best
            model_best = model_tmp
            ic_best = ic_tmp
        end
    end

    return model_best, DataFrame(p = ps, IC = ics)
end

#-------------------------------------------------------------------------------
# SIMULATION
#-------------------------------------------------------------------------------


"""
    function _simulate!(                                   # k variables, T periods
        errors::AbstractMatrix{M},                         # k × T
        B::AbstractMatrix{M};                              # k × kp+m
        trend_exponents::AbstractVector{<:Real}=[0],       # m × 1
        initial::Union{Nothing,AbstractVector{M}}=nothing  # kp × 1
    ) where {M<:Real}

Simulate a VAR(p) overwriting `errors`.
"""
function _simulate!(                                             # k variables, T periods
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
        _rotate_in!(view(Zt, (m+1):(m+kp)), view(errors, :, t))
    end

    return errors
end

# overwrites errors
function simulate!(                                             # k variables, T periods, p lags
    ::Type{VAR},                                                  #  
    errors::AbstractMatrix{<:Number},                           # k  × T
    B::AbstractMatrix{<:Number};                                # k  × kp + m 
    trend_exponents::AbstractVector{<:Number}=[0],              # m  × 1
    initial::Union{Nothing, AbstractVector{<:Number}}=nothing   # kp × 1
)

    errors = _simulate!(errors, B; trend_exponents=trend_exponents, initial=initial)
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
    initial::Union{Nothing, AbstractVector{<:Number}}=nothing
)
    k = size(B, 1) 
    errors = cholesky(Sigma_u).L * randn(k, T) 
    return simulate!(VAR, errors, B; trend_exponents=trend_exponents, initial=initial)
end
