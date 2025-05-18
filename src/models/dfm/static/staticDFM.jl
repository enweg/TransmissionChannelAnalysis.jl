
"""
Static form of a DFM

```math
\\begin{split}
X_t &= \\Lambda F_t + e_t \\\\
e_{i,t} &= \\delta_i(L) e_{i, t-1} + \\nu_{i, t} \\\\
F_t &= \\Phi(L)F_{t-1} + G\\eta_t \\\\
G &= [I_q; O] \\\\
\\mathbb{E}(\\nu_t\\eta_{t-k}') &= 0 \\quad \\forall k
\\end{split}
```
- ``r`` static factors 
- ``N`` cross section 
- ``T`` time dimension
"""
mutable struct StaticDFM <: DFM
    # DGP
    Lambda::AbstractMatrix                    # Factor loadings
    Phi::AbstractVector{<:AbstractMatrix}     # Factor lag polynomial (p-long)
    p::Int
    delta::AbstractVector{<:AbstractMatrix}   # Observation error lag-polynomial (q-long)
    q::Int
    num_factors::Int                          # Number of Factors (coudl be fixed or estimated)

    # Additional Estimated Quantities
    F::AbstractMatrix                         # Factors
    E_factor::AbstractMatrix                  # Factor residuals
    E_observation::AbstractMatrix             # Observation Residuals
    X_hat::AbstractMatrix                     # Fitted Observations

    # Inputs
    input_data::DataFrame                     # Data used to fit model
end
function StaticDFM(data::DataFrame, num_factors::Int, p::Int, q::Int)
    type = eltype(data[:, 1])

    Lambda = Matrix{type}(undef, 0, 0)
    Phi = [Matrix{type}(undef, 0, 0) for _ in 1:p]
    delta = [Matrix{type}(undef, 0, 0) for _ in 1:q]

    F = Matrix{type}(undef, 0, 0)
    E_factor = Matrix{type}(undef, 0, 0)
    E_observation = Matrix{type}(undef, 0, 0)
    X_hat = Matrix{type}(undef, 0, 0)

    return StaticDFM(
        Lambda, Phi, p, delta, q, num_factors, F, E_factor,
        E_observation, X_hat, data
    )
end

# TODO: test 
is_fitted(model::StaticDFM) = (size(model.X_hat, 1) > 0)
coeffs(model::StaticDFM) = require_fitted(model) && (model.Lambda, model.Phi, model.delta)
get_factor_loadings(model::StaticDFM) = require_fitted(model) && model.Lambda
fitted(model::StaticDFM) = require_fitted(model) && model.X_hat
get_factors(model::StaticDFM) = require_fitted(model) && model.F
residuals(model::StaticDFM) = require_fitted(model) && (model.E_observation, model.E_factor)
nobs(model::StaticDFM) = size(model.input_data, 1)
get_input_data(model::StaticDFM) = model.input_data
is_structural(model::StaticDFM) = false
get_variable_names(model::StaticDFM) = names(get_input_data(model))

"""
Simulate the dynamic form of the DFM
- `nu` becomes the obervations
- `G_eta` becomes the factors (corresponds to ``G\\eta_t``)
"""
function _simulate!(
    ::Type{StaticDFM},
    nu::AbstractMatrix,                      # innovations observation error equation
    G_eta::AbstractMatrix,                   # innovations static factor (G\\eta_t)
    Lambda::AbstractMatrix,                  # factor loading matrix
    Phi::AbstractVector{<:AbstractMatrix},   # lag polynomial for static factors
    delta::AbstractVector{<:AbstractMatrix}  # lag polynomial observation error equation
)

    # FIX: this is an ugly quick fix for the situation in which Phi or delta 
    # are empty. Find a better solution before release.
    if isempty(delta)
        delta = [zeros(size(nu, 1), size(nu, 1))]
    end
    if isempty(Phi)
        Phi = [zeros(size(Lambda, 2), size(Lambda, 2))]
    end


    Phi_mat = lag_poly_to_matrix(Phi)
    delta_mat = lag_poly_to_matrix(delta)

    z_f = zeros(size(Phi_mat, 2))
    z_e = zeros(size(delta_mat, 2))

    size(nu, 2) == size(G_eta, 2) || error("Inconsistent definition of time periods.")
    T = size(nu, 2)
    M = eltype(nu)

    for t = 1:T
        # factors 
        mul!(view(G_eta, :, t), Phi_mat, z_f, one(M), one(M))
        _rotate_in!(z_f, view(G_eta, :, t))
        # observation errors
        mul!(view(nu, :, t), delta_mat, z_e, one(M), one(M))
        _rotate_in!(z_e, view(nu, :, t))
        # observations
        mul!(view(nu, :, t), Lambda, view(G_eta, :, t), one(M), one(M))
    end

    return G_eta, nu  # factors, observations
end

# TODO: document
function simulate!(
    ::Type{StaticDFM},
    nu::AbstractMatrix,                               # innovations observation error equation
    G_eta::AbstractMatrix,                            # innovations static factor (G\\eta_t)
    Lambda::AbstractMatrix,                           # factor loading matrix
    Phi::AbstractVector{<:AbstractMatrix},            # lag polynomial for static factors
    delta::AbstractVector{<:AbstractMatrix}=Matrix[]  # lag polynomial observation error equation
)
    _simulate!(StaticDFM, nu, G_eta, Lambda, Phi, delta)

    F = G_eta
    data = DataFrame(nu', :auto)
    p = length(Phi)
    q = length(delta)

    return StaticDFM(data, size(F, 2), p, q)
end
function simulate(
    ::Type{StaticDFM},
    T::Int,                                           # number of periods
    Lambda::AbstractMatrix,                           # factor loading matrix
    Phi::AbstractVector{<:AbstractMatrix},            # lag polynomial for static factors
    delta::AbstractVector{<:AbstractMatrix}=Matrix[]  # lag polynomial observation error equation
)

    n_observed = size(Lambda, 1)
    n_factors = size(Lambda, 2)

    nu = randn(n_observed, T)
    G_eta = randn(n_factors, T)

    return simulate!(StaticDFM, nu, G_eta, Lambda, Phi, delta)
end
