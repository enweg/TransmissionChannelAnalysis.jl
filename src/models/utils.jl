using LinearAlgebra
using DataFrames

"""
    make_companion_matrix(B_plus::Matrix{<:Number}) -> Matrix{<:Number}
    make_companion_matrix(Bs::Vector{<:Matrix{<:Number}}, p::Int, n_exo::Int) -> Matrix{<:Number}

Constructs the companion matrix of a VAR(p) model.

# Arguments
- `Bs`: a vector of lag matrices `[B_1, B_2, ..., B_p]`, each `k × k`
- `B_plus`: a single matrix formed by horizontally concatenating ``[C B_1 B_2 ... B_p]``
   where ``C`` is the matrix of coefficients for deterministic (exogeneous)
   components.
- `p::Int`: lag-length of the VAR.
- `n_exo::Int`: number of exogenous components, i.e. number of columns in ``C``.

# Returns
- `Matrix{<:Number}`: The companion matrix of size `(k*p × k*p)`

# Notes
The companion matrix `C` has the block form:
```math
C = \\begin{bmatrix}
    B_+ \\\\
    I_{(p-1)k} & 0
\\end{bmatrix}
```
where ``B_+`` stacks the lag matrices, and ``k`` is the number of variables.
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


"""
    make_lag_matrix!(X::AbstractMatrix{<:Number}, Y::AbstractMatrix{<:Number}) -> AbstractMatrix{<:Number}
    make_lag_matrix(Y::AbstractMatrix{<:Number}, nlag::Int) -> AbstractMatrix{<:Number}

Constructs a lag matrix from `Y`, storing the results in `X` (mutating) or
returning a new matrix.

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
"""
make_lag_matrix, make_lag_matrix!

function make_lead_matrix!(X::AbstractMatrix{<:Number}, Y::AbstractVecOrMat{<:Number})
    view(X, :, :) .= eltype(X)(NaN)
    k = size(Y, 2)
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
    make_lead_matrix!(X::AbstractMatrix{<:Number}, Y::AbstractVecOrMat{<:Number}) -> AbstractMatrix{<:Number}
    make_lead_matrix(Y::AbstractVecOrMat{<:Number}, nlead::Int) -> AbstractMatrix{<:Number}

Constructs a lead matrix from `Y`, storing the results in `X` (mutating) or
returning a new matrix.

If `Y` is a matrix of time series data with `T` observations and `k` variables,
then the lead matrix `X` has along its first `k` columns data `Y` led by
one period, along columns `(k+1):2k` data `Y` led by two periods, etc.

Missing values (due to leads) are filled with `NaN` to preserve shape and
type consistency.

# Arguments
- `X::AbstractMatrix{<:Number}`: A preallocated matrix of dimensions
  `T × (k × nlead)` to store lead values
- `Y::AbstractVecOrMat{<:Number}`: A matrix of time series data with dimensions
  `T × k`
- `nlead::Int`: The number of leads to include

# Returns
- If using `make_lead_matrix!(X, Y)`, the function modifies `X` in place and returns it
- If using `make_lead_matrix(Y, nlead)`, a new matrix containing the lead values is created and returned

# Examples
**Mutating version**
```julia
k = 2;
T = 10;
Y = reduce(hcat, fill(collect(1:T), k));
nlead = 2;
X = zeros(T, size(Y, 2) * nlead);
make_lead_matrix!(X, Y)
```

**Non-mutating version**
```julia
k = 2;
T = 10;
Y = reduce(hcat, fill(collect(1:T), k));
nlead = 2;
X = make_lead_matrix(Y, nlead)
```
"""
make_lead_matrix, make_lead_matrix!

"""
    _make_trend!(v::AbstractVector, t::Int, trend_exponents::AbstractVector{<:Real}) --> AbstractVector

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
    _rotate_in!(x::AbstractVector{T}, x_new::AbstractVector{T}) --> AbstractVector{T}

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

"""
    _find_variable_idx(variable::Union{Symbol, Int}, data::DataFrame) -> Int

Returns the column index of a variable in a `DataFrame`.

If the input `variable` is an `Int`, it is assumed to already represent
the column index and is returned directly.

# Arguments
- `variable::Union{Symbol, Int}`: Column name (`Symbol`) or column index (`Int`)
- `data::DataFrame`: The DataFrame in which to search

# Returns
- `Int`: The column index of the specified variable
"""
function _find_variable_idx(variable::Union{Symbol,Int}, data::DataFrame)
    isa(variable, Int) && return variable
    return findfirst(==(variable), Symbol.(names(data)))
end

"""
    _separate_lag_matrices(B::Matrix{<:Number}, p::Int) -> Vector{Matrix{<:Number}}

Internal utility to convert a stacked VAR coefficient matrix into a vector of
individual lag matrices.

Given a matrix `B_plus = [B_1 B_2 ... B_p]` that horizontally stacks the
lag coefficient matrices of a VAR(p) model (excluding deterministic components),
this function returns a vector containing each `B_i` as a separate matrix.

# Arguments
- `B::Matrix{<:Number}`: The stacked coefficient matrix of size `(k × k*p)`,
  where `k` is the number of variables and `p` the lag order
- `p::Int`: The number of lags in the VAR model

# Returns
- `Vector{Matrix{<:Number}}`: A vector `[B_1, B_2, ..., B_p]` where each
  `B_i` is a `k × k` matrix

"""
function _separate_lag_matrices(B::AbstractMatrix, p::Int)
    k = div(size(B, 2), p)
    Bs = [B[:, ((i-1)*k+1):(i*k)] for i = 1:p]
    return Bs
end

"""
    _2sls(X::AbstractMatrix{<:Number},
         Y::AbstractMatrix{<:Number},
         Z::AbstractMatrix{<:Number})

Internal function for two-stage least squares (2SLS) estimation.

Estimates regression coefficients when some regressors in `X` may be
endogenous, using instruments in `Z`. Commonly used in external
instrument identification.

# Arguments
- `X::Matrix{<:Number}`: regressors (may include endogenous variables)
- `Y::Matrix{<:Number}`: outcomes
- `Z::Matrix{<:Number}`: instruments and exogenous regressors

# Returns
- `Matrix{<:Number}`: 2SLS coefficient estimates
"""
function _2sls(
    X::AbstractMatrix{<:Number},
    Y::AbstractMatrix{<:Number},
    Z::AbstractMatrix{<:Number}
)

    X_hat = Z * ((Z' * Z) \ (Z' * X))
    return (X_hat' * X_hat) \ (X_hat' * Y)
end

function _find_nonmissing_period(X::AbstractMatrix)
    good_obs = map(row -> !any(ismissing.(row)), eachrow(X))
    bad_obs = (!).(good_obs)

    idx_first = findfirst(good_obs)
    idx_last = findfirst(bad_obs[(idx_first+1):end])
    idx_last = isnothing(idx_last) ? size(X, 1) : idx_last + idx_first - 1

    all(good_obs[idx_first:idx_last]) || error("Some values are missing inbetween.")

    return idx_first, idx_last
end

function _find_data_overlap(mats::AbstractMatrix...)
    N = size(mats[1], 1)
    idx_first = 1
    idx_last = N
    for mat in mats
        size(mat, 1) == N || throw(ArgumentError("All matrices must have same number of rows."))

        idx_first_tmp, idx_last_tmp = _find_nonmissing_period(mat)
        idx_first = max(idx_first, idx_first_tmp)
        idx_last = min(idx_last, idx_last_tmp)
    end

    return (mat[idx_first:idx_last, :] for mat in mats)
end
