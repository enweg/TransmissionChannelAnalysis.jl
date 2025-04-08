using TransmissionChannelAnalysis
using LinearAlgebra
using Test
using DataFrames

using Random
Random.seed!(6150533)

@testset "SVAR basic functions" begin
    k = 3
    p = 2
    T = 100_000
    trend_exponents = [0]

    # First create RF model
    Sigma_u = fill(0.5, k, k)
    Sigma_u[diagind(Sigma_u)] .= 1
    B_plus = 0.2 * randn(k, k * p + length(trend_exponents))
    # Obtaining structural system from this
    A0 = inv(cholesky(Sigma_u).L)
    A_plus = A0 * B_plus

    model = simulate(SVAR, T, A0, A_plus; trend_exponents=trend_exponents)
    @test is_fitted(model) == false
    get_dependent(model)
    get_independent(model)
    get_input_data(model)
    @test nobs(model) == T - p

    # should all throw an error because SVAR has not been estimated
    @test_throws "SVAR must first" coeffs(model)
    @test_throws "SVAR must first" fitted(model)
    @test_throws "SVAR must first" residuals(model)
    @test_throws "SVAR must first" aic(model)
    @test_throws "SVAR must first" is_stable(model)

    # fitting  using recursive
    fit!(model, Recursive())
    fitted(model)
    residuals(model)
    aic(model)
    is_stable(model)

    # fitting using internal instrument throws error
    @test_throws "Internal instruments can only" fit!(model, InternalInstrument(2))
end

@testset "SVAR coefficient recovery Recursive identification" begin
    k = 3
    p = 2
    T = 100_000
    trend_exponents = [0]

    # First create RF model
    Sigma_u = fill(0.5, k, k)
    Sigma_u[diagind(Sigma_u)] .= 1
    B_plus = 0.2 * randn(k, k * p + length(trend_exponents))
    # Obtaining structural system from this
    A0 = inv(cholesky(Sigma_u).L)
    A_plus = A0 * B_plus

    # testing recovery of coefficients using Recursive
    initial = fill(100, k * p)
    shocks = zeros(k, T)
    model = simulate(SVAR, shocks, A0, A_plus; trend_exponents=trend_exponents, initial=initial)
    # cannot do identification because covariance is zero. But can fit the RF model
    model_var = model.var
    fit!(model_var)
    @test maximum(abs, coeffs(model_var) - inv(A0) * A_plus) < 1e-6
    # now give it correct covariance
    Phi0 = inv(A0)
    model_var.Sigma_u = Phi0 * Phi0'
    # and use internal function to get identified params back
    # because fit!, identify, identify! simply use this function, we do not need 
    # to test them separately (besides that they run)
    A0_hat, A_plus_hat = TransmissionChannelAnalysis._identify(model_var, Recursive())
    @test maximum(abs, A0_hat - A0) < 1e-6
    @test maximum(abs, A_plus_hat - A_plus) < 1e-6
    # we can also check it better by using the true matrices directly
    m = length(trend_exponents)
    A0_hat, A_plus_hat = TransmissionChannelAnalysis._identify(B_plus, Sigma_u, Recursive())
    @test maximum(abs, A0_hat - A0) < sqrt(eps())
    @test maximum(abs, A_plus_hat - A_plus) < sqrt(eps())
end

@testset "SVAR Information Criteria" begin
    k = 3
    p = 2
    T = 100_000
    trend_exponents = [0]

    # First create RF model
    Sigma_u = fill(0.5, k, k)
    Sigma_u[diagind(Sigma_u)] .= 1
    B_plus = 0.2 * randn(k, k * p + length(trend_exponents))
    # Obtaining structural system from this
    A0 = inv(cholesky(Sigma_u).L)
    A_plus = A0 * B_plus

    model = simulate(SVAR, T, A0, A_plus; trend_exponents=trend_exponents)

    # checking information criteria
    data = get_input_data(model)
    model = SVAR(data, 10; trend_exponents=trend_exponents)
    model_best, ic_table = fit_and_select!(model, Recursive(), aic)
    @test model_best.p == p
    model = SVAR(data, 10; trend_exponents=trend_exponents)
    model_best, ic_table = fit_and_select!(model, Recursive(), bic)
    @test model_best.p == p
    model = SVAR(data, 10; trend_exponents=trend_exponents)
    model_best, ic_table = fit_and_select!(model, Recursive(), hqc)
    @test model_best.p == p
    model = SVAR(data, 10; trend_exponents=trend_exponents)
    model_best, ic_table = fit_and_select!(model, Recursive(), sic)
    @test model_best.p == p
end

@testset "SVAR IRFs Recursive Identification" begin

    k = 3
    p = 2
    T = 100_000
    trend_exponents = [0]

    # First create RF model
    Sigma_u = fill(0.5, k, k)
    Sigma_u[diagind(Sigma_u)] .= 1
    B_plus = 0.2 * randn(k, k * p + length(trend_exponents))
    # Obtaining structural system from this
    A0 = inv(cholesky(Sigma_u).L)
    A_plus = A0 * B_plus

    # Testing IRFs Recursive
    max_horizon = 10
    model = simulate(SVAR, T, A0, A_plus; trend_exponents=trend_exponents)
    @test_throws "SVAR must first" IRF(model, max_horizon)
    fit!(model, Recursive())
    irf = IRF(model, max_horizon)

    A0_hat, _ = coeffs(model)
    Phi0 = inv(A0_hat)
    C = make_companion_matrix(coeffs(model.var), p, length(trend_exponents))
    irfs_companion = [(C^i)[1:k, 1:k] * Phi0 for i = 0:max_horizon]
    irfs_companion = cat(irfs_companion...; dims=3)
    @test maximum(abs, irf.irfs - irfs_companion) < sqrt(eps())

    # testing internal recursive irf functions separately
    max_horizon = 10
    m = length(trend_exponents)
    irfs = TransmissionChannelAnalysis._identify_irfs(B_plus[:, (m+1):end], Sigma_u, p, Recursive(), max_horizon)

    Phi0 = inv(A0)
    C = make_companion_matrix(B_plus, p, m)
    irfs_companion = [(C^i)[1:k, 1:k] * Phi0 for i = 0:max_horizon]
    irfs_companion = cat(irfs_companion...; dims=3)
    @test maximum(abs, irfs_companion - irfs) < sqrt(eps())

    # IRFs can also be directly identified from a VAR using Recursive
    model = simulate(VAR, T, B_plus, Sigma_u; trend_exponents=trend_exponents)
    @test_throws "VAR must first" TransmissionChannelAnalysis._identify_irfs(model, Recursive(), max_horizon)
    fit!(model)
    irfs = TransmissionChannelAnalysis._identify_irfs(model, Recursive(), max_horizon)
    C = make_companion_matrix(coeffs(model), p, length(trend_exponents))
    A0_hat, _ = TransmissionChannelAnalysis._identify(model, Recursive())
    Phi0 = inv(A0_hat)
    irfs_companion = [(C^i)[1:k, 1:k] * Phi0 for i = 0:max_horizon]
    irfs_companion = cat(irfs_companion...; dims=3)
    @test maximum(abs, irfs - irfs_companion) < sqrt(eps())

    # Lastly the implementation test
    irf = IRF(model, Recursive(), max_horizon)
end

@testset "SVAR IRFs InternalInstrument Identification" begin
    k = 3
    p = 2
    T = 1_000_000
    trend_exponents = [0]
    Phi0 = randn(k, k)
    A0 = inv(Phi0)
    # next three steps make sure that diag of A0 is positive
    S = diagm(sign.(diag(A0)))
    Phi0 = S' * Phi0
    A0 = inv(Phi0)

    B_plus = 0.2 * randn(k, k * p + length(trend_exponents))
    A_plus = A0 * B_plus
    Sigma_u = Phi0 * Phi0'
    max_horizon = 10

    # introducing instrument into SVAR
    # first variable is instrument
    # first "shock" is measurement error
    Phi0_tilde = zeros(k + 1, k + 1)
    Phi0_tilde[2:end, 2:end] = Phi0
    Phi0_tilde[1, [1, 2]] = [0.1, 1]

    B_plus_tilde = zeros(k + 1, (k + 1) * p + length(trend_exponents))
    m = length(trend_exponents)
    B_plus_tilde[2:end, 1:m] = B_plus[:, 1:m]
    B = view(B_plus, :, (m+1):size(B_plus, 2))
    B_tilde = view(B_plus_tilde, 2:size(B_plus_tilde, 1), (m+1):size(B_plus_tilde, 2))
    for pp = 1:p
        Bi = view(B, :, (k*(pp-1)+1):(pp*k))
        Bi_tilde = view(B_tilde, :, ((k+1)*(pp-1)+1):(pp*(k+1)))
        view(Bi_tilde, :, 2:size(Bi_tilde, 2)) .= Bi
    end
    B_plus_tilde

    Sigma_u_tilde = Phi0_tilde * Phi0_tilde'

    method = InternalInstrument(2)
    irfs = TransmissionChannelAnalysis._identify_irfs(B_plus_tilde[:, (m+1):end], Sigma_u_tilde, p, method, max_horizon)

    irfs_true = TransmissionChannelAnalysis._svar_irf(A0, A_plus[:, (m+1):end], p, max_horizon)

    @test maximum(abs, irfs[2:end, 1:1, :] - irfs_true[:, 1:1, :] ./ irfs_true[1, 1, 1]) < sqrt(eps())

    # testing it all through simulation
    T = 1_000_000
    shocks = randn(k, T)
    model = simulate(SVAR, shocks, A0, A_plus; trend_exponents=trend_exponents)
    data = get_input_data(model)
    data[!, :instrument] = shocks[1, :] + 0.1 * randn(T)
    select!(data, :instrument, :)
    model_var = VAR(data, p; trend_exponents=trend_exponents)
    fit!(model_var)

    method = InternalInstrument(2)
    irfs_simulated = IRF(model_var, method, max_horizon).irfs

    @test maximum(abs, irfs_simulated[2:end, 1:1, :] - irfs_true[:, 1:1, :] ./ irfs_true[1, 1, 1]) < 1e-2
end

@testset "SVAR transmission implementation" begin
    # This is just an implementation test. All underlying functions 
    # have already been tested elsewhere. 

    k = 3
    p = 2
    T = 100_000
    trend_exponents = [0]

    # First create RF model
    Sigma_u = fill(0.5, k, k)
    Sigma_u[diagind(Sigma_u)] .= 1
    B_plus = 0.2 * randn(k, k * p + length(trend_exponents))
    # Obtaining structural system from this
    A0 = inv(cholesky(Sigma_u).L)
    A_plus = A0 * B_plus

    model = simulate(SVAR, T, A0, A_plus; trend_exponents=trend_exponents)
    fit!(model, Recursive())

    transmission_order = [1, 3, 2]
    q = make_condition("!y_{2,0}", transmission_order)
    transmission_effects = transmission(1, model, q, transmission_order, 3)
end
