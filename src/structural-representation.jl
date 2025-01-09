"""
Use the Cholesky decomposition to find ``L`` and ``D`` as used in $WEGNER.
"""
function make_L_D(Sigma::AbstractMatrix{<:Real})
    Linv = cholesky(Sigma).L
    L = inv(Linv)
    D = diagm(1 ./ diag(L))
    return L, D
end

"""
Construct B of the static representation ``x = Bx + Omega\\varepsilon``.

# Notes
- See `make_static_model` for further information.
"""
function make_B(
    As::Vector{<:AbstractMatrix},
    Sigma::AbstractMatrix{<:Real},
    order::AbstractVector{<:Int},
    max_horizon::Int
)

    # 1. Creating the transmission matrix
    T = permmatrix(order)
    As = [T * A * T' for A in As]

    # 2. Cholesky decomposition
    L, D = make_L_D(T * Sigma * T')
    K = size(Sigma, 1)  # assuming as many shocks as variables
    As = [D * L * A for A in As]  # this gives DQ'A_i^* in the paper

    # 3. creating B
    row_block = hcat(reduce(hcat, reverse(As)), I - D * L)
    B = zeros(K * (max_horizon + 1), K * (max_horizon + 1))
    return slide_in!(B, row_block)

end

"""
Construct Omega of the static representation ``x = Bx + Omega\\varepsilon``.

# Notes
- See `make_static_model` for further information.
"""
function make_Omega(
    Phi0::AbstractMatrix,
    Psis::Vector{<:AbstractMatrix},
    Sigma::AbstractMatrix{<:Real},
    order::AbstractVector{<:Int},
    max_horizon::Int
)

    T = permmatrix(order)
    K = size(Sigma, 1)
    L, D = make_L_D(T * Sigma * T')
    Qt = L * T * Phi0
    Psis = [D * L * T * Psi * Phi0 for Psi in Psis]
    row_block = hcat(reduce(hcat, reverse(Psis)), D * Qt)
    Omega = zeros(K * (max_horizon + 1), K * (max_horizon + 1))
    return slide_in!(Omega, row_block)
end

"""
    make_static_model(
        Phi0::AbstractMatrix, 
        As::Vector{<:AbstractMatrix}, 
        Psis::Vector{<:AbstractMatrix}, 
        Sigma::AbstractMatrix{<:Real}, 
        order::AbstractVector{<:Int}, 
        max_horizon::Int
    ) -> Tuple{AbstractMatrix, AbstractMatrix}

Transform an SVARMA dynamic model into the static 
representation ``x = Bx + Omega\\varepsilon``.

# Arguments
- `Phi0::AbstractMatrix`: The matrix of contemporaneous structural impulse 
  responses.
- `As::Vector{<:AbstractMatrix}`: A vector of reduced-form autoregressive (AR).
  First entry corresponds to the AR matrix for the first lag, etc. 
- `Psis::Vector{<:AbstractMatrix}`: A vector of reduced-form moving average (MA) 
  matrices. First index corresponds to the first lag, etc.
- `Sigma::AbstractMatrix{<:Real}`: The covariance matrix of reduced-form errors.
- `order::AbstractVector{<:Int}`: The vector indicating the order of variables, 
  typically determined by the transmission matrix.
- `max_horizon::Int`: The maximum time horizon to consider for the static model, 
  with `0` representing the contemporaneous period.

# Returns
- `(B, Omega)` with the meaning being the same as in $WEGNER.

# Notes
- Use `make_B` and `make_Omega` to construct the two matrices seperately.

"""
function make_static_model(
    Phi0::AbstractMatrix, 
    As::Vector{<:AbstractMatrix}, 
    Psis::Vector{<:AbstractMatrix}, 
    Sigma::AbstractMatrix{<:Real}, 
    order::AbstractVector{<:Int}, 
    max_horizon::Int
)

    B = make_B(As, Sigma, order, max_horizon)
    Omega = make_Omega(Phi0, Psis, Sigma, order, max_horizon)
    return B, Omega
end
