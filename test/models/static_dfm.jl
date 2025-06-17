using TransmissionChannelAnalysis
using LinearAlgebra
using Random
using Statistics

@testset "DFM simulation + estimation" begin
    Random.seed!(6150533)
    N = 100  # number of variables
    r = 2     # number of factors
    T = 100_000  # number of time periods
    trend_exponents = Real[]
    p = 1  # number of lags in the factor VAR

    merrors = 1e-4 * randn(N, T)
    varerrors = randn(r, T)
    B = diagm(fill(0.3, r))  # Makes sure that factors are uncorrelated
    Lambda = Matrix(qr(randn(N, r)).Q)  # creating an orthogonal matrix

    model = simulate!(
        DFM, merrors, varerrors, B, Lambda;
        trend_exponents=trend_exponents
    )

    data = get_input_data(model)
    F = model.F
    S = diagm(1 ./ vec(std(Matrix(data); dims=1)))
    Lambda = S * Lambda
    Y = Matrix(data) * S

    model = DFM(data, p, r)
    fit!(model)
    Y_hat = fitted(model)
    Y_hat_true = F * Lambda'

    # Given the variance of the input data is roughly one, a mistake of 1e-2
    # can be considered very small and unlikely to happen out of luck
    @test maximum(abs, Y_hat .- Y_hat_true) < 1e-2

    # --- CHECKING WHETHER TREND EXPONENTS AND HIGHER P ARE ALSO WORKING
    Random.seed!(6150533)
    N = 100  # number of variables
    r = 2     # number of factors
    T = 100_000  # number of time periods
    trend_exponents = Real[]
    p = 2  # number of lags in the factor VAR
    merrors = 1e-4 * randn(N, T)
    varerrors = randn(r, T)
    B = reduce(hcat, diagm(fill(0.3, r)) for _ in 1:p)
    B = hcat(randn(r, length(trend_exponents)), B)
    Lambda = Matrix(qr(randn(N, r)).Q)  # creating an orthogonal matrix
    model = simulate!(
        DFM, merrors, varerrors, B, Lambda;
        trend_exponents=trend_exponents
    )

    data = get_input_data(model)
    F = model.F
    S = diagm(1 ./ vec(std(Matrix(data); dims=1)))
    Lambda = S * Lambda
    Y = Matrix(data) * S

    model = DFM(data, p, r)
    fit!(model)
    Y_hat = fitted(model)
    Y_hat_true = F * Lambda'

    @test maximum(abs, Y_hat .- Y_hat_true) < 1e-2

    # --- IMPLEMENTATION TEST: CHECKING IF SIMULATE ALSO WORKS
    model = simulate(DFM, T, B, Lambda; trend_exponents=trend_exponents)
    # --- IMPLEMENTATION TEST: CHECKING IF OTHER TREND EXPONENTS WORK
    Random.seed!(6150533)
    N = 100  # number of variables
    r = 2     # number of factors
    T = 100_000  # number of time periods
    trend_exponents = [0]
    p = 2  # number of lags in the factor VAR
    B = reduce(hcat, diagm(fill(0.3, r)) for _ in 1:p)
    B = hcat(randn(r, length(trend_exponents)), B)
    Lambda = Matrix(qr(randn(N, r)).Q)  # creating an orthogonal matrix
    model = simulate(DFM, T, B, Lambda; trend_exponents=trend_exponents)
end

@testset "DFM Implementation Tests" begin
    Random.seed!(6150533)
    N = 100  # number of variables
    r = 2     # number of factors
    T = 100_000  # number of time periods
    trend_exponents = [0]
    p = 2  # number of lags in the factor VAR
    B = reduce(hcat, diagm(fill(0.3, r)) for _ in 1:p)
    B = hcat(randn(r, length(trend_exponents)), B)
    Lambda = Matrix(qr(randn(N, r)).Q)  # creating an orthogonal matrix
    model = simulate(DFM, T, B, Lambda; trend_exponents=trend_exponents)

    @test !is_fitted(model)
    @test_throws "DFM must first" factors(model)
    @test_throws "DFM must first" loadings(model)
    @test_throws "DFM must first" coeffs(model)
    @test_throws "DFM must first" fitted(model)
    @test_throws "DFM must first" residuals(model)
    @test_throws "DFM must first" get_factor_var(model)

    fit!(model)
    @test is_fitted(model)
    @test factors(model) == model.F
    @test loadings(model) == model.Lambda
    @test coeffs(model) == (model.Lambda, coeffs(model.factor_var))
    @test fitted(model) == model.Yhat
    @test residuals(model) == (model.eta_hat, residuals(model.factor_var))
    @test get_factor_var(model) === model.factor_var
end

@testset "DFM IRF" begin
    Random.seed!(6150533)
    N = 100  # number of variables
    r = 2     # number of factors
    T = 100_000  # number of time periods
    trend_exponents = Real[]
    p = 1  # number of lags in the factor VAR

    merrors = 1e-4 * randn(N, T)
    varerrors = randn(r, T)
    B = diagm(fill(0.3, r))  # Makes sure that factors are uncorrelated
    Lambda = Matrix(qr(randn(N, r)).Q)  # creating an orthogonal matrix

    model = simulate!(
        DFM, merrors, varerrors, B, Lambda;
        trend_exponents=trend_exponents
    )

    data = get_input_data(model)
    F = model.F
    S = diagm(1 ./ vec(std(Matrix(data); dims=1)))
    Lambda = S * Lambda

    max_horizon = 4
    @test_throws "DFM must first" IRF(model, max_horizon)
    fit!(model)
    # Giving it the true values
    model.factor_var.B = B
    model.Lambda = Lambda
    irfs_vars, irfs_factors = IRF(model, max_horizon)

    # VAR irfs are tested elsewhere
    irfs_factors_true = TransmissionChannelAnalysis._var_irf(B, p, max_horizon)
    irfs_vars_true = mapslices(x -> Lambda * x, irfs_factors_true; dims=[1, 2])

    @test isapprox(irfs_vars.irfs, irfs_vars_true; atol=sqrt(eps()))
    @test isapprox(irfs_factors.irfs, irfs_factors_true; atol=sqrt(eps()))
end
