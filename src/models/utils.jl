using LinearAlgebra
using DataFrames

"""
    make_companion_matrix(Bs::AbstractVector{<:AbstractMatrix{<:Number}}) -> Matrix{<:Number}

Constructs the companion matrix of a Vector Autoregression (VAR) model.

# Description
In a VAR(p) model, the system can be rewritten in its companion form, which is 
useful for analyzing stability properties and computing impulse responses. 
Given the VAR representation:

```math
y_t = B_1 y_{t-1} + \\dots + B_p y_{t-p} + u_t,
```

we define the alternative representation:

```math
y_t = B^+ x_t + u_t,
```

where:
- ``B^+ = [B_1 B_2 \\dots B_p]`` is the coefficient matrix of lagged values,
- ``x_t = (y_{t-1}', \\dots, y_{t-p}')'`` is the stacked state vector.

The function constructs the standard companion matrix ``C``, which is given by:

```math
C =
\\begin{bmatrix}
B^+ \\\\
I_{(p-1)k} & 0
\\end{bmatrix}
```

where ``k`` is the number of endogenous variables and ``p`` is the number of lags.

# Arguments
- `Bs::AbstractVector{AbstractMatrix{<:Number}}`: A vector containing the 
   coefficient matrices ``[B_1, B_2, \\dots, B_p]``, where each 
   ``B_i`` is a ``k \\times k`` matrix.

# Returns
- `C::Matrix{<:Number}`: The companion matrix of size ``(kp \\times kp)``.

# Example
```julia
Bs = [rand(2,2) for _ in 1:2]
C = make_companion_matrix(Bs)
```
"""
function make_companion_matrix(Bs::AbstractVector{<:AbstractMatrix{<:Number}})
    B_plus = reduce(hcat, Bs)
    p = length(Bs)
    return make_companion_matrix(B_plus, p, 0)
end
function make_companion_matrix(B_plus::AbstractMatrix{<:Number}, p::Int, n_exo::Int)
    B_plus = B_plus[:, (n_exo+1):end]
    k = size(B_plus, 1)
    C = diagm(-k => ones(k * p - k))
    C[1:k, :] .= B_plus
    return C
end



"""
    spectral_radius(X::AbstractMatrix) -> Number

Computes the spectral radius of a matrix `X`.

# Description
The spectral radius of a matrix is defined as the largest absolute eigenvalue 
of the matrix. 

Mathematically, for a given square matrix `X`, the spectral radius is computed as:

```math
\\rho(X) = \\max |\\lambda_i|
```
where `\\lambda_i` are the eigenvalues of `X`.

# Arguments
- `X::AbstractMatrix`: A square matrix. 

# Returns
- The spectral radius of `X`.

# Example
```julia
X = [0.5 0.2; 0.1 0.9]
r = spectral_radius(X)
```
"""
spectral_radius(X::AbstractMatrix) = maximum(abs, eigvals(X))


"""
    make_lag_matrix!(X::AbstractMatrix{<:Number}, Y::AbstractMatrix{<:Number}) -> AbstractMatrix{<:Number}
    make_lag_matrix(Y::AbstractMatrix{<:Number}, nlag::Int) -> AbstractMatrix{<:Number}

Constructs a lag matrix from `Y`, storing the results in `X` (mutating) or 
returning a new matrix.

# Description

If `Y` is a matrix of time series data with `T` observations and `k` variables, 
then the lag matrix `X` has along its first `k` columns data `Y` lagged by 
one period, along columns `(k+1):2k` data `Y` lagged by two periods, etc. 

Missing values (due to lags) are filled with `NaN` to maintain dimensional 
and type consistency.

# Arguments
- `X::AbstractMatrix{<:Number}`: A preallocated matrix of dimensions 
  `T × (k \\times nlag)` to store lagged values.
- `Y::AbstractMatrix{<:Number}`: A matrix of time series data with dimensions 
  `T × k`.
- `nlag::Int`: The number of lags to include.

# Returns
- If using `make_lag_matrix!(X, Y)`, the function modifies `X` in place and returns it.
- If using `make_lag_matrix(Y, nlag)`, a new matrix containing the lagged values 
  is created and returned.

# Examples

**Mutating version**

```julia
k = 2;
T = 10;
Y = collect(1:T);
Y = reduce(hcat, fill(Y, k));
nlag = 2;
X = zeros(T, size(Y, 2) * nlag);
make_lag_matrix!(X, Y)
```

**Non-mutating version**

```julia
k = 2;
T = 10;
Y = collect(1:T);
Y = reduce(hcat, fill(Y, k));
nlag = 2;
X = make_lag_matrix(Y, nlag)
```

# Notes
- The function assumes that `Y` has more rows than `nlag`, ensuring enough data points for lagging.
- The output maintains the same number of rows as `Y`, with leading `NaN` values where necessary.
- The mutating version (`make_lag_matrix!`) is useful for performance-sensitive applications, avoiding extra allocations.
"""
function make_lag_matrix!(X::AbstractMatrix{<:Number}, Y::AbstractVecOrMat{<:Number})
    view(X, :, :) .= eltype(X)(NaN)
    k = size(Y, 2)
    nlag = floor(Int, size(X, 2) / k)
    for l = 1:nlag
        @views X[(l+1):end, ((l-1)*k+1):(l*k)] .= Y[1:(end-l), :]
    end
    return X
end
function make_lag_matrix(Y::AbstractVecOrMat{<:Number}, nlag::Int)
    # do it separately because otherwise it breaks with vector input
    r = size(Y, 1)
    c = size(Y, 2)
    T = eltype(Y) <: Int ? Float64 : eltype(Y)
    X = zeros(T, r, c * nlag)
    return make_lag_matrix!(X, Y)
end

# TODO: document
# TODO: write official tests
function make_lead_matrix!(X::AbstractMatrix{<:Number}, Y::AbstractVecOrMat{<:Number})
    view(X, :, :) .= eltype(X)(NaN)
    k = size(Y,2)
    nlead = floor(Int, size(X, 2) / k)
    for l = 1:nlead
        @views X[1:(end-l), ((l-1)*k+1):(l*k)] .= Y[(l+1):end, :]
    end
    return X
end
function make_lead_matrix(Y::AbstractVecOrMat{<:Number}, nlead::Int)
    r = size(Y, 1)
    c = size(Y, 2)
    T = eltype(Y) <: Int ? Float64 : eltype(Y)
    X = zeros(T, r, c * nlead)
    return make_lead_matrix!(X, Y)
end

"""
    _make_trend!(v::AbstractVector, t::Int, trend_exponents::AbstractVector{<:Real})

Given a set of `trend_exponents`, fill the vector `v` with `[t^e for e in trend_exponents]`
without any allocations. 
"""
function _make_trend!(v::AbstractVector, t::Int, trend_exponents::AbstractVector)
    for i in 1:lastindex(trend_exponents)
        v[i:i] .= t^trend_exponents[i]
    end
    return v
end

"""
Rotate in the new values `x_new` at the beginning of the old `x`. This works
by shifting all values in `x` to higher indices such that the subvector `x_new`
fits into the beginning of `x`. 

Notes:
    - This is used for simulation and other places where the lag vector needs to 
    be updated repeatedly
"""
function _rotate_in!(x::AbstractVector{T}, x_new::AbstractVector{T}) where {T}
    k = length(x_new)
    for i = lastindex(x):-1:(k+1)
        x[i] = x[i-k]
    end
    copyto!(view(x, 1:k), x_new)
    return x
end

function _find_variable_idx(variable::Union{Symbol, Int}, data::DataFrame)
    isa(variable, Int) && return variable
    return findfirst(==(variable), names(data))
end
