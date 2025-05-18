@testset "DFM utils" begin
    @testset "lag_poly_to_matrix" begin
        lag_poly = [ones(3, 3), 2 * ones(3, 3), 3 * ones(3, 3)]
        lag_matrix = TransmissionChannelAnalysis.lag_poly_to_matrix(lag_poly)
        @test lag_matrix == [lag_poly[1] lag_poly[2] lag_poly[3]]
    end
end

@testset "DFM PCA" begin
    # 1. Factor loadings are orthonormal
    X = randn(10_000, 4)
    pca = PCA(X, 2)
    @test maximum(abs, pca.lambda' * pca.lambda - I) < 1e-6
    # 2. Factors are uncorrelated
    @test maximum(abs, tril(pca.F' * pca.F, -1)) < 1e-6

    # 3. If number of factors is the same as number of variables, 
    # the same reconstructed X should be approximately the true X.
    X = randn(10_000, 4)
    pca = PCA(X, 4)
    Xhat = pca.F * pca.lambda'
    @test maximum(abs, Xhat - X) < sqrt(eps())

    # 4. Recovery of factors and loadings up to rotation
    # If Fhat are the recovered factors, then the rotation can be 
    # recovered via H = (F' * F) \ (F' * Fhat). We can then test 
    # whether we recover the correct loading. 
    F = randn(10_000, 2)
    lambda = randn(4, 2)
    X = F * lambda'
    pca = PCA(X, 2)
    Fhat = pca.F
    H = (F' * F) \ (F' * Fhat)
    @test maximum(abs, H * pca.lambda' - lambda') < sqrt(eps())
end

@testset "DFM PCA IC" begin
    Random.seed!(6150533)
    T = 10_000
    m = 3
    n = 20

    F = randn(T, m)
    lambda = randn(n, m)
    X = F * lambda'

    pca = PCA(X, m)
    @test TransmissionChannelAnalysis._msr(pca) < sqrt(eps())
    @test TransmissionChannelAnalysis._msr(pca, m) < sqrt(eps())
    TransmissionChannelAnalysis._msr(pca, 2)
    TransmissionChannelAnalysis._msr(pca, 1)

    BaiNgPC1(pca, 1, m)
    BaiNgPC2(pca, 1, m)
    BaiNgPC3(pca, 1, m)

    BaiNgIC1(pca, 1, m)
    BaiNgIC2(pca, 1, m)
    BaiNgIC3(pca, 1, m)

    pca = PCA(X, 1)
    pca, ic_table = select_factors!(pca, m+2, BaiNgPC1)
    @test size(pca.F, 2) == m
    pca, ic_table = select_factors!(pca, m+2, BaiNgPC2)
    @test size(pca.F, 2) == m
    pca, ic_table = select_factors!(pca, m+2, BaiNgPC3)
    @test size(pca.F, 2) == m
    pca, ic_table = select_factors!(pca, m+2, BaiNgIC1)
    @test size(pca.F, 2) == m
    pca, ic_table = select_factors!(pca, m+2, BaiNgIC2)
    @test size(pca.F, 2) == m
    pca, ic_table = select_factors!(pca, m+2, BaiNgIC3)
    @test size(pca.F, 2) == m

    # adding a tiny bit of errors
    Random.seed!(6150533)
    T = 10_000
    m = 3
    n = 20
    F = randn(T, m)
    lambda = randn(n, m)
    X = F * lambda' + 0.1 * randn(T, n)
    pca = PCA(X, 1)
    pca, ic_table = select_factors!(pca, m+2, BaiNgPC1)
    @test size(pca.F, 2) == m
    pca, ic_table = select_factors!(pca, m+2, BaiNgPC2)
    @test size(pca.F, 2) == m
    pca, ic_table = select_factors!(pca, m+2, BaiNgPC3)
    @test size(pca.F, 2) == m
    pca, ic_table = select_factors!(pca, m+2, BaiNgIC1)
    @test size(pca.F, 2) == m
    pca, ic_table = select_factors!(pca, m+2, BaiNgIC2)
    @test size(pca.F, 2) == m
    pca, ic_table = select_factors!(pca, m+2, BaiNgIC3)
    @test size(pca.F, 2) == m

    # Below is just a check if it also runs if there is some noise
    # No official test is being run. 
    # Consider this an implementation test.
    X = F * lambda' + 0.1 * randn(T, n)
    pca = PCA(X, 1)
    pca, ic_table = select_factors!(pca, m+2, BaiNgPC1)
end

@testset "StaticDFM Basic Functions" begin
    # ALL TESTS HERE ARE ONLY IMPLEMENTATION TESTS
    
    N = 10
    r = 4
    T = 1_000
    nu = randn(N, T)
    G_eta = randn(r, T)

    Lambda = 0.2 * randn(N, r)
    Phi = [0.2 * randn(r, r), 0.2 * randn(r, r)]
    delta = [0.5 * diagm(randn(N))]

    model = simulate!(StaticDFM, nu, G_eta, Lambda, Phi, delta)

    @test !is_fitted(model) 
    @test_throws "StaticDFM must first" coeffs(model)
    @test_throws "StaticDFM must first" get_factor_loadings(model)
    @test_throws "StaticDFM must first" fitted(model)
    @test_throws "StaticDFM must first" get_factors(model)
    @test_throws "StaticDFM must first" residuals(model)
    @test nobs(model) == T
    get_input_data(model)
    @test !is_structural(model)
    get_variable_names(model)

    # TODO: test after fitting the model
end

function _simulation_staticDFM_alt(
    nu::AbstractMatrix,
    G_eta::AbstractMatrix,
    Lambda::AbstractMatrix,
    Phi::AbstractVector{<:AbstractMatrix},
    delta::AbstractVector{<:AbstractMatrix}
)

    T = size(nu, 2)
    N = size(nu, 1)
    r = size(G_eta, 1)

    X = zeros(N, T)
    e = zeros(N, T + length(delta))
    F = zeros(r, T + length(Phi))

    for t = 1:T
        # factors
        for i = 1:length(Phi)
            F[:, t+length(Phi)] .+= Phi[i] * F[:, t+length(Phi)-i]
        end
        F[:, t+length(Phi)] .+= G_eta[:, t]
        # observation innovations
        for i = 1:length(delta)
            e[:, t+length(delta)] .+= delta[i] * e[:, t+length(delta)-i]
        end
        e[:, t+length(delta)] .+= nu[:, t]
        # observations
        X[:, t] .= Lambda * F[:, t+length(Phi)] + e[:, t+length(delta)]
    end
    
    return X, F[:, (length(Phi)+1):end]
end

@testset "StaticDFM Simulation" begin
    Random.seed!(6150533)
    N = 10
    r = 4
    T = 1_000
    nu = randn(N, T)
    G_eta = randn(r, T)

    Lambda = 0.2 * randn(N, r)
    Phi = [0.2 * randn(r, r), 0.2 * randn(r, r)]
    delta = [0.5 * diagm(randn(N))]

    F, X = TransmissionChannelAnalysis._simulate!(StaticDFM, copy(nu), copy(G_eta), Lambda, Phi, delta) 
    X_test, F_test = _simulation_staticDFM_alt(nu, G_eta, Lambda, Phi, delta)

    @test maximum(abs, F_test - F) < sqrt(eps())
    @test maximum(abs, X - X_test) < sqrt(eps())

    # no autocorrelation in observation errors
    Random.seed!(6150533)
    N = 10
    r = 4
    T = 1_000
    nu = randn(N, T)
    G_eta = randn(r, T)

    Lambda = 0.2 * randn(N, r)
    Phi = [0.2 * randn(r, r), 0.2 * randn(r, r)]
    delta = Matrix{Float64}[]

    F, X = TransmissionChannelAnalysis._simulate!(StaticDFM, copy(nu), copy(G_eta), Lambda, Phi, delta) 
    X_test, F_test = _simulation_staticDFM_alt(nu, G_eta, Lambda, Phi, delta)

    @test maximum(abs, F_test - F) < sqrt(eps())
    @test maximum(abs, X - X_test) < sqrt(eps())
    
    # VAR(0) for factors
    Random.seed!(6150533)
    N = 10
    r = 4
    T = 1_000
    nu = randn(N, T)
    G_eta = randn(r, T)

    Lambda = 0.2 * randn(N, r)
    Phi = Matrix{Float64}[]
    delta = [0.5 * diagm(randn(N))]

    F, X = TransmissionChannelAnalysis._simulate!(StaticDFM, copy(nu), copy(G_eta), Lambda, Phi, delta) 
    X_test, F_test = _simulation_staticDFM_alt(nu, G_eta, Lambda, Phi, delta)

    @test maximum(abs, F_test - F) < sqrt(eps())
    @test maximum(abs, X - X_test) < sqrt(eps())

    # implementation tests
    N = 10
    r = 4
    T = 1_000
    nu = randn(N, T)
    G_eta = randn(r, T)

    Lambda = 0.2 * randn(N, r)
    Phi = [0.2 * randn(r, r), 0.2 * randn(r, r)]
    delta = [0.5 * diagm(randn(N))]

    model = simulate!(StaticDFM, nu, G_eta, Lambda, Phi, delta)
    model = simulate(StaticDFM, T, Lambda, Phi, delta)
end

function _simulation_dynamicDFM_alt(
    nu::AbstractMatrix,
    eta::AbstractMatrix,
    lambda::AbstractVector{<:AbstractMatrix},
    Psi::AbstractVector{<:AbstractMatrix},
    delta::AbstractVector{<:AbstractMatrix}
)

    T = size(nu, 2)
    N = size(nu, 1)
    q = size(eta, 1)

    X = zeros(N, T)
    e = zeros(N, T + length(delta))
    f = zeros(q, T + length(Psi))

    for t = 1:T
        # factors 
        for i = 1:length(Psi)
            f[:, t+length(Psi)] .+= Psi[i] * f[:, t+length(Psi)-i]
        end
        f[:, t+length(Psi)] .+= eta[:, t]
        # observation errors
        for i = 1:length(delta)
            e[:, t+length(delta)] .+= delta[i] * e[:, t+length(delta)-i]
        end
        e[:, t+length(delta)] .+= nu[:, t]
        # observations
        for i = 1:length(lambda)
            X[:, t] .+= lambda[i] * f[:, t+length(Psi)-i+1]  # + 1 because of contemporaneous term
        end
        X[:, t] .+= e[:, t+length(delta)]
    end

    return X, f[:, (end-T+1):end]
end

@testset "DynamicDFM Simulation" begin
    Random.seed!(6150533)
    N = 10
    q = 2
    T = 1_000
    nu = randn(N, T)
    eta = randn(q, T)
    delta = [0.5 * diagm(randn(N))]
    lambda = [0.2 * randn(N, q)]
    Psi = [0.2 * randn(q, q), 0.02 * randn(q, q)]
    f, X = TransmissionChannelAnalysis._simulate!(DynamicDFM, copy(nu), copy(eta), lambda, Psi, delta) 
    X_test, f_test = _simulation_dynamicDFM_alt(nu, eta, lambda, Psi, delta)
    @test maximum(abs, f_test - f) < sqrt(eps())
    @test maximum(abs, X - X_test) < sqrt(eps())

    Random.seed!(6150533)
    N = 200
    q = 10
    T = 1_000
    nu = randn(N, T)
    eta = randn(q, T)
    delta = [0.2 * diagm(randn(N))]
    lambda = [0.02 * randn(N, q)]
    Psi = [0.02 * randn(q, q), 0.02 * randn(q, q)]
    f, X = TransmissionChannelAnalysis._simulate!(DynamicDFM, copy(nu), copy(eta), lambda, Psi, delta) 
    X_test, f_test = _simulation_dynamicDFM_alt(nu, eta, lambda, Psi, delta)
    @test maximum(abs, f_test - f) < sqrt(eps())
    @test maximum(abs, X - X_test) < sqrt(eps())
end
