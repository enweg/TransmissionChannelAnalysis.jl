
"""
```math
\\\\begin{split}
X_t &= \\lambda(L)f_t + e_t \\\\
e_{i,t} &= \\delta_i(L) e_{i, t-1} + \\nu_{i, t} \\\\
f_t &= \\Psi(L)f_{t-1} + \\eta_{t} \\\\
\\mathbb{E}(\\nu_t\\eta_{t-k}') &= 0 \\quad \\forall k
\\\\end{split}
```
- ``q`` dynamic factors 
- ``N`` cross section
- ``T`` time dimension
"""
mutable struct DynamicDFM <: DFM
    lambda::AbstractVector{<:AbstractMatrix}
    Psi::AbstractVector{<:AbstractMatrix}
    delta::AbstractVector{<:AbstractMatrix}
end

"""
Simulate the dynamic form of the DFM
- `nu` becomes the obervations
- `eta` becomes the factors
"""
function _simulate!(
    ::Type{DynamicDFM}, 
    nu::AbstractMatrix,                         # innovations in observation error equation
    eta::AbstractMatrix,                        # innovations in factor equation
    lambda::AbstractVector{<:AbstractMatrix},   # lag polynomial observation equation
    Psi::AbstractVector{<:AbstractMatrix},      # lag polynomial factor equation
    delta::AbstractVector{<:AbstractMatrix}     # lag polynomial observation error equation
)

    delta_mat = lag_poly_to_matrix(delta)
    Psi_mat = lag_poly_to_matrix(Psi)
    lambda_mat = lag_poly_to_matrix(lambda)

    z_f = zeros(size(Psi_mat, 2))
    z_e = zeros(size(delta_mat, 2))
    z_x = zeros(size(lambda_mat, 2))

    size(nu, 2) == size(eta, 2) || error("Inconsistent definition of time periods.")
    T = size(nu, 2)
    M = eltype(nu)
    for t = 1:T
        # factors 
        mul!(view(eta, :, t), Psi_mat, z_f, one(M), one(M))
        _rotate_in!(z_f, view(eta, :, t))
        _rotate_in!(z_x, view(eta, :, t))  # because X also depends on contemporaneous factors
        # observation errors
        mul!(view(nu, :, t), delta_mat, z_e, one(M), one(M))
        _rotate_in!(z_e, view(nu, :, t))
        # observations
        mul!(view(nu, :, t), lambda_mat, z_x, one(M), one(M))
    end

    return eta, nu  # factors, observations    
end
 

