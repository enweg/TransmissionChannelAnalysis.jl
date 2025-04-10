using TransmissionChannelAnalysis
using Test
using Random
using LinearAlgebra
Random.seed!(6150533)

@testset "Basic VAR functions" begin
    k = 3
    p = 2
    T = 10_000
    trend_exponents = [0]
    B = 0.2 * randn(k, k*p + length(trend_exponents))

    model = simulate(VAR, T, B; trend_exponents=trend_exponents)
    @test is_fitted(model) == false
    get_dependent(model)
    get_independent(model)
    get_input_data(model)
    @test nobs(model) == T - p

    # following functions should all through an error because model has not been 
    # fitted yet
    @test_throws "VAR must first be estimated" coeffs(model)
    @test_throws "VAR must first" cov(model)
    @test_throws "VAR must first" fitted(model)
    @test_throws "VAR must first" residuals(model)
    @test_throws "VAR must first" aic(model)
    @test_throws "VAR must first" is_stable(model)
    @test_throws "VAR must first" make_companion_matrix(model)

    # fitting the model
    fit!(model)
    @test is_fitted(model) == true
    # should no-longer throw an error
    coeffs(model)
    cov(model)
    fitted(model)
    residuals(model)
    aic(model)
    is_stable(model)
    make_companion_matrix(model)
end

@testset "VAR coefficient recovery" begin
    k = 3
    p = 2
    T = 10_000
    trend_exponents = [0]
    B = 0.2 * randn(k, k*p + length(trend_exponents))

    # testing if estimation and simulation are correct by trying to recover
    # original coefficients if errors are all zero
    errors = zeros(k, T)
    initial = fill(100, k*p)
    model = simulate!(VAR, errors, B; trend_exponents=trend_exponents, initial = initial)
    fit!(model)
    @test maximum(abs, coeffs(model) - B) < 1e-6
end

@testset "VAR covariance recovery" begin
    k = 3
    p = 2
    T = 10_000_000
    trend_exponents = 0:1
    B = 0.2 * randn(k, k*p + length(trend_exponents))
    Sigma_u = [
        10 5 5
        5 10 5 
        5 5 10
    ]
    model = simulate(VAR, T, B, Sigma_u; trend_exponents=trend_exponents)
    fit!(model)
    @test maximum(abs, cov(model) - Sigma_u) < 1e-2
end

@testset "VAR information criteria" begin
    k = 3
    p = 2
    T = 10_000
    trend_exponents = [0]
    B = 0.2 * randn(k, k*p + length(trend_exponents))

    # Checking information criteria
    model = simulate(VAR, T, B; trend_exponents=trend_exponents)
    data = get_input_data(model)
    model = VAR(data, 10; trend_exponents=trend_exponents)
    model_best, ic_table = fit_and_select!(model, aic)
    @test model_best.p == p
    model_best, ic_table = fit_and_select!(model, bic)
    @test model_best.p == p
    model_best, ic_table = fit_and_select!(model, hqc)
    @test model_best.p == p
    model_best, ic_table = fit_and_select!(model, sic)
    @test model_best.p == p

    # Same checks but this time non identity covariance
    Sigma_u = [
        1 0.5 0.5;
        0.5 1 0.5; 
        0.5 0.5 1
    ]
    model = simulate(VAR, T, B, Sigma_u; trend_exponents=trend_exponents)
    data = get_input_data(model)
    model = VAR(data, 10; trend_exponents=trend_exponents)
    model_best, ic_table = fit_and_select!(model, aic)
    @test model_best.p == p
    model_best, ic_table = fit_and_select!(model, bic)
    @test model_best.p == p
    model_best, ic_table = fit_and_select!(model, hqc)
    @test model_best.p == p
    model_best, ic_table = fit_and_select!(model, sic)
    @test model_best.p == p
end

@testset "VAR trend exponents implementation test" begin

    k = 3
    p = 2
    T = 1_000

    # constant and linear trend
    trend_exponents = 0:1
    B = 0.2 * randn(k, k*p + length(trend_exponents))
    model = simulate(VAR, T, B; trend_exponents=trend_exponents)
    fit!(model)

    # polynomial trend
    trend_exponents = 0:2
    B = 0.2 * randn(k, k*p + length(trend_exponents))
    model = simulate(VAR, T, B; trend_exponents=trend_exponents)
    fit!(model)

    # no constant no trend
    trend_exponents = Real[]
    B = 0.2 * randn(k, k*p + length(trend_exponents))
    model = simulate(VAR, T, B; trend_exponents=trend_exponents)
    fit!(model)
    
    # only trend
    trend_exponents = [1]
    B = 0.2 * randn(k, k*p + length(trend_exponents))
    model = simulate(VAR, T, B; trend_exponents=trend_exponents)
    fit!(model)
end

@testset "VAR lag length implementation test" begin
    
    k = 3
    trend_exponents = [0]
    T = 10_000

    # no lags
    p = 0
    B = 2 * randn(k, k*p + length(trend_exponents))
    model = simulate(VAR, T, B; trend_exponents=trend_exponents)
    fit!(model)
    coeffs(model) - B

    # no lags
    p = 1
    B = 0.2 * randn(k, k*p + length(trend_exponents))
    model = simulate(VAR, T, B; trend_exponents=trend_exponents)
    fit!(model)

    # no lags
    p = 10
    B = 0.05 * randn(k, k*p + length(trend_exponents))
    model = simulate(VAR, T, B; trend_exponents=trend_exponents)
    fit!(model)
end

@testset "VAR number of variables implementation test" begin
    p = 2
    trend_exponents = [0]
    T = 10_000

    # AR
    k = 1
    B = 0.2 * randn(k, k*p + length(trend_exponents))
    model = simulate(VAR, T, B; trend_exponents=trend_exponents)
    fit!(model)

    k = 2
    B = 0.2 * randn(k, k*p + length(trend_exponents))
    model = simulate(VAR, T, B; trend_exponents=trend_exponents)
    fit!(model)

    # no lags
    k = 20
    B = 0.05 * randn(k, k*p + length(trend_exponents))
    model = simulate(VAR, T, B; trend_exponents=trend_exponents)
    fit!(model)
end

@testset "VAR IRFs" begin
    
    k = 3
    p = 2
    T = 10_000
    trend_exponents = [0]
    B = 0.2 * randn(k, k*p + length(trend_exponents))

    model = simulate(VAR, T, B; trend_exponents=trend_exponents)
    fit!(model)

    max_horizon = 10
    irf = IRF(model, max_horizon)
    # computing IRFs manually using companion form
    B_hat = coeffs(model, false)
    C = make_companion_matrix(B_hat, p, length(trend_exponents))
    irfs_manual = [(C^i)[1:k, 1:k] for i=0:max_horizon]
    irfs_manual = cat(irfs_manual...; dims=3)
    @test maximum(abs, irf.irfs - irfs_manual) < sqrt(eps())

    # Same thing but not deterministic
    trend_exponents = Real[]
    B = 0.2 * randn(k, k*p + length(trend_exponents))
    model = simulate(VAR, T, B; trend_exponents=trend_exponents)
    fit!(model)
    max_horizon = 10
    irf = IRF(model, max_horizon)
    B_hat = coeffs(model, false)
    C = make_companion_matrix(B_hat, p, length(trend_exponents))
    irfs_manual = [(C^i)[1:k, 1:k] for i=0:max_horizon]
    irfs_manual = cat(irfs_manual...; dims=3)
    @test maximum(abs, irf.irfs - irfs_manual) < sqrt(eps())
    
    # Same thing but more deterministic
    trend_exponents = 0:6
    B = 0.2 * randn(k, k*p + length(trend_exponents))
    model = simulate(VAR, T, B; trend_exponents=trend_exponents)
    fit!(model)
    max_horizon = 10
    irf = IRF(model, max_horizon)
    B_hat = coeffs(model, false)
    C = make_companion_matrix(B_hat, p, length(trend_exponents))
    irfs_manual = [(C^i)[1:k, 1:k] for i=0:max_horizon]
    irfs_manual = cat(irfs_manual...; dims=3)
    @test maximum(abs, irf.irfs - irfs_manual) < sqrt(eps())
    
    # only one variable
    trend_exponents = 0:1
    k = 1
    B = 0.2 * randn(k, k*p + length(trend_exponents))
    model = simulate(VAR, T, B; trend_exponents=trend_exponents)
    fit!(model)
    max_horizon = 10
    irf = IRF(model, max_horizon)
    B_hat = coeffs(model, false)
    C = make_companion_matrix(B_hat, p, length(trend_exponents))
    irfs_manual = [(C^i)[1:k, 1:k] for i=0:max_horizon]
    irfs_manual = cat(irfs_manual...; dims=3)
    @test maximum(abs, irf.irfs - irfs_manual) < sqrt(eps())
end

@testset "VAR transmission implementation" begin
    # This is only a implementation test. All underlying functions 
    # have previously been tested. 

    k = 3
    p = 2
    T = 10_000
    trend_exponents = [0]
    B = 0.2 * randn(k, k*p + length(trend_exponents))

    model = simulate(VAR, T, B; trend_exponents=trend_exponents)
    fit!(model)

    transmission_order = [3, 1, 2]
    q = make_condition("!y_{1,0} & !y_{1,1}", transmission_order)
    transmission_effects = transmission(2, model, Recursive(), q, transmission_order, 3)

    # Obviously there is no valid instrument here, but we check the correctness 
    # of the instrument estimation elsewhere. Here we really just care about 
    # whether the functions run without error. Correctness, again, has been 
    # shown elsewhere. 
    method = InternalInstrument(2)
    transmission_effects = transmission(1, model, method, q, transmission_order, 3)
end
