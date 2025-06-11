
function _get_lambda_F(X::AbstractMatrix, eigvecs::AbstractMatrix, m::Int)
    lambda = view(eigvecs, :, 1:m)
    F = X * lambda
    return (lambda, F)
end

# TODO: document
mutable struct PCA
    # X = F * lambda' + E
    # X_t = lambda * F_t + e_t
    X::AbstractMatrix
    eigvals::AbstractVector
    eigvecs::AbstractMatrix
    lambda::AbstractMatrix
    F::AbstractMatrix
end
function PCA(X::AbstractMatrix, m::Int)
    T = size(X, 1)
    Sigma_X = Symmetric(X' * X / T)
    evals, evecs = eigen(Sigma_X; sortby=x -> -x)  # ordering largest to smallest
    lambda, F = _get_lambda_F(X, evecs, m)
    # lambda = view(evecs, :, 1:m)
    # F = X * lambda

    return PCA(X, evals, evecs, lambda, F)
end
PCA(data::DataFrame, m::Int) = PCA(Matrix(data), m)


factors(pca::PCA) = pca.F
loadings(pca::PCA) = pca.lambda
get_input_data(pca::PCA) = pca.X
fitted(pca::PCA) = pca.F * pca.lambda'
residuals(pca::PCA) = pca.X - fitted(pca)

# TODO: document
function PCA!(pca::PCA, m::Int)
    # changes the number of selected PCs
    lambda, F = _get_lambda_F(pca.X, pca.eigvecs, m)
    pca.lambda = lambda
    pca.F = F
end

function _msr(pca::PCA)
    m = size(pca.F, 2)
    return _msr(pca, m)
end
function _msr(pca::PCA, m::Int)
    # m is the selected number of PCs
    lambda, F = _get_lambda_F(pca.X, pca.eigvecs, m)
    Xhat = F * lambda'
    return sum((pca.X - Xhat) .^ 2) / length(Xhat)
end

# TODO: document all the functions below
function BaiNg_g1(n, T)
    # n = number of variables
    # T = number of time periods / observations
    return ((n + T) / (n * T)) * log((n * T) / (n + T))
end
function BaiNg_g2(n, T)
    # n = number of variables
    # T = number of time periods / observations
    return ((n + T) / (n * T)) * log(min(n, T))
end
function BaiNg_g3(n, T)
    # n = number of variables
    # T = number of time periods / observations
    return log(min(n, T)) / min(n, T)
end

function BaiNgPC(pca::PCA, m::Int, m_max::Int, g::Function)
    # T = time periods
    # n = number of variables
    T, n = size(pca.X)
    return _msr(pca, m) + m * _msr(pca, m_max) * g(n, T)
end
function BaiNgPC1(pca::PCA, m::Int, m_max::Int)
    return BaiNgPC(pca, m, m_max, BaiNg_g1)
end
function BaiNgPC2(pca::PCA, m::Int, m_max::Int)
    return BaiNgPC(pca, m, m_max, BaiNg_g2)
end
function BaiNgPC3(pca::PCA, m::Int, m_max::Int)
    return BaiNgPC(pca, m, m_max, BaiNg_g3)
end

function BaiNgIC(pca::PCA, m::Int, g::Function)
    T, n = size(pca.X)
    return log(_msr(pca, m)) + m * g(n, T)
end
function BaiNgIC1(pca::PCA, m::Int, m_max::Int=size(pca.X, 2))
    # last argument to be consistent with PC
    return BaiNgIC(pca, m, BaiNg_g1)
end
function BaiNgIC2(pca::PCA, m::Int, m_max::Int=size(pca.X, 2))
    return BaiNgIC(pca, m, BaiNg_g2)
end
function BaiNgIC3(pca::PCA, m::Int, m_max::Int=size(pca.X, 2))
    return BaiNgIC(pca, m, BaiNg_g3)
end

function select_factors!(pca::PCA, m_max::Int, ic::Function)
    # ic(pca, m, m_max) --> Number
    ic_values = zeros(m_max)
    best_m = m_max
    for m = m_max:-1:1
        ic_values[m] = ic(pca, m, m_max)
        if ic_values[m] < ic_values[best_m]
            best_m = m
        end
    end

    PCA!(pca, best_m)
    return pca, DataFrame(:number_factors => 1:m_max, :IC => ic_values)
end

# TODO: document
# Requires CairoMakie to be loaded. In that case PlotsExt will load
function screeplot end

function apply_named_factors(
    Lambda::AbstractMatrix{<:Number},  # factor loadings
    F::AbstractMatrix{<:Number};  # factors
    naming_idx::AbstractVector{<:Int}=1:size(F, 2)  # naming variables
)
    # X_t = Lambda * F_t
    # X = F * Lambda'

    length(unique(naming_idx)) == size(F, 2) || throw(ArgumentError("Number of naming variables < number of factors."))

    Lambda_1 = Lambda[naming_idx, :]
    Lambda_1_inv = inv(Lambda_1)
    Lambda = Lambda * Lambda_1_inv
    F = F * Lambda_1'
    return Lambda, F
end

function apply_named_factors!(pca::PCA; kwargs...)
    Lambda = loadings(pca)
    F = factors(pca)
    Lambda, F = apply_named_factors(Lambda, F; kwargs...)
    pca.lambda = Lambda
    pca.F = F
    return pca
end
