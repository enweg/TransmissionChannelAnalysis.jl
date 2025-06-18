using TransmissionChannelAnalysis
using LinearAlgebra
using Random
using Statistics
using DataFrames

@testset "SDFM simulation + estimation" begin
    # Note that the actual SVAR estimation and identification is tested
    # elsewhere. Thus, all that remains to be tested is whether the
    # estimation runs without errors, whether the named factors are applied,
    # and, even though already tested for PCA seperately, whether the
    # Lambda * F_t is close to the truth -- this part is identified even if
    # the factors and loadings are not seperately identified.

    Random.seed!(6150533)
    N = 100  # number of variables
    r = 2     # number of factors
    T = 100_000  # number of time periods
    trend_exponents = Real[]
    p = 1  # number of lags in the factor VAR

    # Note that for SDFM, factors are named and thus not orthogonal
    merrors = 1e-4 * randn(N, T)
    varshocks = randn(r, T)
    A0 = LowerTriangular(10 * randn(r, r))
    A0 = diagm(sign.(diag(A0))) * A0
    B_plus = diagm(fill(0.3, r))
    A_plus = A0 * B_plus
    Lambda = Matrix(qr(randn(N, r)).Q)  # creating an orthogonal matrix
    Lambda1 = Lambda[1:r, :]
    Lambda = Lambda * inv(Lambda1)

    model = simulate!(
        SDFM, merrors, varshocks, A0, A_plus, Lambda;
        trend_exponents=trend_exponents
    )

    data = get_input_data(model)
    F = model.F
    S = diagm(1 ./ vec(std(Matrix(data); dims=1)))
    Lambda = S * Lambda
    Y = Matrix(data) * S

    model = SDFM(data, p, r; trend_exponents=trend_exponents)
    method = Recursive()
    fit!(model, method)
    Y_hat = fitted(model)
    Y_hat_true = F * Lambda'
    # Variance of data is about one, so an error of 1e-2 is tiny
    @test maximum(abs, Y_hat - Y_hat_true) < 1e-2
    # Factors are named, so first r rows in Lambda must be identity matrix
    @test maximum(abs, loadings(model)[1:r, :] - I) < sqrt(eps())

    # --- CHECKING WHETHER TREND EXPONENTS AND HIGHER P ARE ALSO WORKING
    Random.seed!(6150533)
    N = 100  # number of variables
    r = 2     # number of factors
    T = 100_000  # number of time periods
    trend_exponents = Real[]
    p = 2  # number of lags in the factor VAR
    merrors = 1e-4 * randn(N, T)
    varshocks = randn(r, T)
    A0 = LowerTriangular(10 * randn(r, r))
    A0 = diagm(sign.(diag(A0))) * A0
    B_plus = reduce(hcat, diagm(fill(0.3, r)) for _ in 1:p)
    B_plus = hcat(randn(r, length(trend_exponents)), B_plus)
    A_plus = A0 * B_plus
    Lambda = Matrix(qr(randn(N, r)).Q)  # creating an orthogonal matrix
    Lambda1 = Lambda[1:r, :]
    Lambda = Lambda * inv(Lambda1)
    model = simulate!(
        SDFM, merrors, varshocks, A0, A_plus, Lambda;
        trend_exponents=trend_exponents
    )

    data = get_input_data(model)
    F = model.F
    S = diagm(1 ./ vec(std(Matrix(data); dims=1)))
    Lambda = S * Lambda
    Y = Matrix(data) * S

    model = SDFM(data, p, r; trend_exponents=trend_exponents)
    method = Recursive()
    fit!(model, method)
    Y_hat = fitted(model)
    Y_hat_true = F * Lambda'
    # Variance of data is about one, so an error of 1e-2 is tiny
    @test maximum(abs, Y_hat - Y_hat_true) < 1e-2
    # Factors are named, so first r rows in Lambda must be identity matrix
    @test maximum(abs, loadings(model)[1:r, :] - I) < sqrt(eps())

    # --- IMPLEMENTATION TEST: CHECKING IF SIMULATE ALSO WORKS
    model = simulate(SDFM, T, A0, A_plus, Lambda; trend_exponents=trend_exponents)
    # --- IMPLEMENTATION TEST: CHECKING IF OTHER TREND EXPONENTS WORK
    trend_exponents = [0]
    B_plus = reduce(hcat, diagm(fill(0.3, r)) for _ in 1:p)
    B_plus = hcat(randn(r, length(trend_exponents)), B_plus)
    A_plus = A0 * B_plus
    model = simulate(SDFM, T, A0, A_plus, Lambda; trend_exponents=trend_exponents)
end

@testset "SDFM Implementation Tests" begin
    Random.seed!(6150533)
    N = 100  # number of variables
    r = 2     # number of factors
    T = 100_000  # number of time periods
    trend_exponents = Real[]
    p = 1  # number of lags in the factor VAR
    merrors = 1e-4 * randn(N, T)
    varshocks = randn(r, T)
    A0 = LowerTriangular(10 * randn(r, r))
    A0 = diagm(sign.(diag(A0))) * A0
    B_plus = diagm(fill(0.3, r))
    A_plus = A0 * B_plus
    Lambda = Matrix(qr(randn(N, r)).Q)  # creating an orthogonal matrix
    Lambda1 = Lambda[1:r, :]
    Lambda = Lambda * inv(Lambda1)
    model = simulate!(
        SDFM, merrors, varshocks, A0, A_plus, Lambda;
        trend_exponents=trend_exponents
    )

    @test !is_fitted(model)
    @test_throws "SDFM must first" factors(model)
    @test_throws "SDFM must first" loadings(model)
    @test_throws "SDFM must first" coeffs(model)
    @test_throws "SDFM must first" fitted(model)
    @test_throws "SDFM must first" residuals(model)
    @test_throws "SDFM must first" get_factor_svar(model)

    method = Recursive()
    fit!(model, method)
    @test is_fitted(model)
    @test factors(model) == model.F
    @test loadings(model) == model.Lambda
    @test coeffs(model) == (model.Lambda, coeffs(model.factor_svar))
    @test fitted(model) == model.Yhat
    @test residuals(model) == (model.eta_hat, residuals(model.factor_svar))
    @test get_factor_svar(model) === model.factor_svar
end

@testset "SDFM IRF" begin
    # Note that the actual identification schemes for the factor SVAR are already
    # tested in the tests for SVARs. Therefore, these are not tested here.
    # What we are testing here is whether the implementation of IRFs is correct.

    Random.seed!(6150533)
    N = 100  # number of variables
    r = 2     # number of factors
    T = 100_000  # number of time periods
    trend_exponents = Real[]
    p = 1  # number of lags in the factor VAR
    merrors = 1e-4 * randn(N, T)
    varshocks = randn(r, T)
    varshocks_copy = copy(varshocks)
    A0 = LowerTriangular(10 * randn(r, r))
    A0 = diagm(sign.(diag(A0))) * A0
    B_plus = diagm(fill(0.3, r))
    A_plus = A0 * B_plus
    Lambda = Matrix(qr(randn(N, r)).Q)  # creating an orthogonal matrix
    Lambda1 = Lambda[1:r, :]
    Lambda = Lambda * inv(Lambda1)
    model = simulate!(
        SDFM, merrors, varshocks, A0, A_plus, Lambda;
        trend_exponents=trend_exponents
    )

    data = get_input_data(model)
    max_horizon = 4
    irfs_factors_true = TransmissionChannelAnalysis._svar_irf(
        A0, A_plus, p, max_horizon
    )
    irfs_vars_true = mapslices(x -> Lambda * x, irfs_factors_true; dims=[1, 2])

    # --- DIRECTLY FROM SDFM
    # Only supports Recursive
    method = Recursive()
    model = SDFM(data, p, r; trend_exponents=trend_exponents)
    @test_throws "SDFM must first" IRF(model, max_horizon)
    fit!(model, method)
    # Giving it the true coefficients
    model.factor_svar.A0 = A0
    model.factor_svar.A_plus = A_plus
    model.Lambda = Lambda
    irfs_vars, irfs_factors = IRF(model, max_horizon)
    @test isapprox(irfs_vars.irfs, irfs_vars_true; atol=sqrt(eps()))
    @test isapprox(irfs_factors.irfs, irfs_factors_true; atol=sqrt(eps()))

    # --- FROM DFM USING RECURSIVE
    model = DFM(data, p, r; trend_exponents=trend_exponents)
    method = Recursive()
    @test_throws "DFM must first" IRF(model, method, max_horizon)
    fit!(model)
    # Giving it the true coefficients
    model.factor_var.B = B_plus
    model.factor_var.Sigma_u = inv(A0) * inv(A0)'
    model.Lambda = Lambda
    irfs_vars, irfs_factors = IRF(model, method, max_horizon)

    @test isapprox(irfs_vars.irfs, irfs_vars_true; atol=sqrt(eps()))
    @test isapprox(irfs_factors.irfs, irfs_factors_true; atol=sqrt(eps()))

    # --- FROM DFM USING EXTERNAL INSTRUMENT
    # Only tests if no errors are thrown. Actual identificatio method is tested
    # in the tests for SVAR
    model = DFM(data, p, r; trend_exponents=trend_exponents)
    fit!(model)
    instruments = DataFrame(instrument=varshocks_copy[1, :])
    # Factor names are by default Fi
    method = ExternalInstrument(:F1, instruments)
    IRF(model, method, max_horizon)
    method = ExternalInstrument(1, instruments)
    IRF(model, method, max_horizon)
end
