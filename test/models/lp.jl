using TransmissionChannelAnalysis
using Test
using DataFrames
using LinearAlgebra
using Random

@testset "LP construction" begin
    k = 4
    T = 10
    p = 2

    data = DataFrame(reshape(1:(T*k), T, k), :auto)

    treatment = 1
    horizons = 0:3
    model = LP(data, treatment, p, horizons; include_constant=true)

    # testing if constant was correctly defined
    @test all(model.X[:, 1] .== 1)
    # testing the rest of the matrices
    for (i, h) in enumerate(horizons)
        @test isequal(model.X[1:(end-h), 2:(treatment+1)], model.Y[1:(end-h), 1:treatment, i] .- h)
        for j = 1:p
            @test isequal(model.X[1:(end-h), ((j-1)*k+3):(j*k+2)], model.Y[1:(end-h), :, i] .- h .- j)
        end
    end
end

@testset "LP basic functions" begin
    k = 4
    T = 100
    p = 2

    data = DataFrame(randn(T, k), :auto)

    treatment = 1
    horizons = 0:3
    model = LP(data, treatment, p, horizons; include_constant=true)
    # Not fitted, so the following should throw errors
    @test_throws "LP must first" coeffs(model)
    @test_throws "LP must first" fitted(model)
    @test_throws "LP must first" residuals(model)
    # different observations by horizon
    @test nobs(model) == T - p .- horizons
    # implementation tests (do they error?)
    get_dependent(model)
    get_independent(model)
    get_input_data(model)

    # fitting model (by default Recursive)
    fit!(model)
    # should no-longer throw errors
    coeffs(model)
    fitted(model)
    residuals(model)
    # explicitly defining Recursive
    fit!(model, Recursive())
end

@testset "LP coefficient + IRF estimat Recursive Identification" begin
    Random.seed!(6150533)
    # Following assumes that SVAR is correctly implemented
    k = 3
    p = 2
    T = 10_000
    trend_exponents = [0]
    A0 = LowerTriangular(randn(k, k))
    S = diagm(sign.(diag(A0)))
    A0 = A0 * S
    B_plus = 0.2 * randn(k, k * p + length(trend_exponents))
    A_plus = A0 * B_plus

    model_svar = simulate(SVAR, T, A0, A_plus; trend_exponents=trend_exponents)
    data = get_input_data(model_svar)

    p_large = 200
    model_svar = SVAR(data, p_large; trend_exponents=trend_exponents)
    fit!(model_svar, Recursive())
    max_horizon = 10
    irfs_svar = IRF(model_svar, max_horizon)

    for treatment = 1:k
        model_lp = LP(data, treatment, p_large, 0:max_horizon; include_constant=true)
        irfs_lp = TransmissionChannelAnalysis._identify_irfs(model_lp, Recursive(), max_horizon)

        @test maximum(abs, irfs_lp - irfs_svar.irfs[:, treatment:treatment, :] ./ irfs_svar.irfs[treatment, treatment, 1]) < 1e-2
        # for contemporaneous horizon they have to be the same
        @test maximum(abs, irfs_lp[:, :, 1] - irfs_svar.irfs[:, treatment:treatment, 1] ./ irfs_svar.irfs[treatment, treatment, 1]) < sqrt(eps())

        # testing only sub-horizons
        model_lp = LP(data, treatment, p_large, [0, 4, 8]; include_constant=true)
        fit!(model_lp, Recursive())
        irfs_lp = coeffs(model_lp, true)[:, treatment:treatment, :]

        @test maximum(abs, irfs_lp[:, :, 1] - irfs_svar.irfs[:, treatment:treatment, 1] ./ irfs_svar.irfs[treatment, treatment, 1]) < sqrt(eps())

        @test maximum(abs, irfs_lp[:, :, 2] - irfs_svar.irfs[:, treatment:treatment, model_lp.horizons[2]+1] ./ irfs_svar.irfs[treatment, treatment, 1]) < 1e-2

        @test maximum(abs, irfs_lp[:, :, 3] - irfs_svar.irfs[:, treatment:treatment, model_lp.horizons[3]+1] ./ irfs_svar.irfs[treatment, treatment, 1]) < 1e-2
    end

    # implementation test
    treatment = 1
    model_lp = LP(data, treatment, p_large, 0:max_horizon; include_constant=true)
    IRF(model_lp, Recursive(), max_horizon)
    model_lp = LP(data, treatment, p_large, [0, 3, 5]; include_constant=true)
    @test_throws "LP horizons do not" IRF(model_lp, Recursive(), max_horizon)
end

@testset "LP Information Critaria" begin
    Random.seed!(6150533)
    # Following assumes that SVAR is correctly implemented
    k = 3
    p = 2
    T = 10_000
    trend_exponents = [0]
    A0 = LowerTriangular(randn(k, k))
    S = diagm(sign.(diag(A0)))
    A0 = A0 * S
    B_plus = 0.2 * randn(k, k * p + length(trend_exponents))
    A_plus = A0 * B_plus

    model_svar = simulate(SVAR, T, A0, A_plus; trend_exponents=trend_exponents)
    data = get_input_data(model_svar)

    p_max = 10
    model_lp = LP(data, 1, p_max, 0:10; include_constant=true)
    model_best, ic_table = fit_and_select!(model_lp, Recursive(), aic)
    @test model_best.p == 2
    model_lp = LP(data, 1, p_max, 0:10; include_constant=true)
    model_best, ic_table = fit_and_select!(model_lp, Recursive(), sic)
    @test model_best.p == 2
    model_lp = LP(data, 1, p_max, 0:10; include_constant=true)
    model_best, ic_table = fit_and_select!(model_lp, Recursive(), bic)
    @test model_best.p == 2
    model_lp = LP(data, 1, p_max, 0:10; include_constant=true)
    model_best, ic_table = fit_and_select!(model_lp, Recursive(), hqc)
    @test model_best.p == 2
end

@testset "LP coefficient + IRF External Instrument" begin
    Random.seed!(6150533)
    # Following assumes that SVAR is correctly implemented
    k = 3
    p = 2
    T = 1_000_000
    trend_exponents = [0]
    A0 = LowerTriangular(randn(k, k))
    S = diagm(sign.(diag(A0)))
    A0 = A0 * S
    B_plus = 0.2 * randn(k, k * p + length(trend_exponents))
    A_plus = A0 * B_plus
    shocks = randn(k, T)

    model_svar = simulate(SVAR, shocks, A0, A_plus; trend_exponents=trend_exponents)
    data = get_input_data(model_svar)
    data[!, :instrument] = shocks[1, :]
    select!(data, :instrument, :Y1, :)

    max_horizon = 3
    m = length(trend_exponents)
    irfs_true = TransmissionChannelAnalysis._svar_irf(A0, A_plus[:, (m+1):end], p, max_horizon)
    irfs_true = irfs_true[:, 1:1, :] ./ irfs_true[1, 1, 1]

    treatment = :Y1
    model = LP(data, :Y1, p, 0:max_horizon; include_constant=true)
    method = ExternalInstrument(treatment, [1])
    irfs_lp = TransmissionChannelAnalysis._identify_irfs(model, method, max_horizon)
    irfs_lp = irfs_lp[2:end, :, :]
    @test maximum(abs, irfs_lp - irfs_true) < 1e-2

    # creating more than one instrument
    data[!, :instrument2] = shocks[1, :] + 0.1 * randn(T)
    data[!, :instrument3] = shocks[1, :] + 0.1 * randn(T)
    data[!, :instrument4] = shocks[1, :] + 0.1 * randn(T)
    select!(data, r"instrument", :Y1, :)

    treatment = :Y1
    model = LP(data, :Y1, p, 0:max_horizon; include_constant=true)
    method = ExternalInstrument(treatment, 1:4)
    irfs_lp = TransmissionChannelAnalysis._identify_irfs(model, method, max_horizon)
    irfs_lp = irfs_lp[5:end, :, :]
    @test maximum(abs, irfs_lp - irfs_true) < 1e-2

    select!(data, r"Y")
    data[!, :instrument] = shocks[1, :]
    select!(data, :instrument, :Y1, :)

    model = LP(data, :Y1, p, 0:max_horizon; include_constant=true)
    method = ExternalInstrument(treatment, [1]; normalising_horizon=1)
    irfs_lp = TransmissionChannelAnalysis._identify_irfs(model, method, max_horizon)
    # The Y1 irf is no-longer the same as SVAR IRF because Y1 is lead by 1
    irfs_lp = irfs_lp[2:end, :, :]

    irfs_true = TransmissionChannelAnalysis._svar_irf(A0, A_plus[:, (m+1):end], p, max_horizon)
    irfs_true = irfs_true[1:end, 1:1, :] ./ irfs_true[1, 1, 2]

    # Estimation is better judged as percentage error here. 
    @test maximum(abs, (irfs_lp - irfs_true) ./ irfs_true) < 1e-2
end

@testset "LP transmission implementation" begin
    # This is just an implementation test. Correctness of intermediate functions 
    # has been tested elsewhere. 

    k = 4
    T = 100
    p = 2

    data = DataFrame(randn(T, k), :auto)

    treatment = 1
    horizons = 0:3
    model = LP(data, treatment, p, horizons; include_constant=true)

    transmission_order = [3, 1, 2, 4]
    q = make_condition("!y_{1,0} & !y_{1,1}", transmission_order)
    transmission_effect = transmission(model, Recursive(), 1, q, transmission_order, maximum(horizons))

    # Oviously the first variable is not a valid instrument, but we are 
    # really just testing whether the function runs. We are not testing for 
    # correctness. That is done elsewhere. 
    method = ExternalInstrument(2, [1])
    transmission_effect = transmission(model, method, 1, q, transmission_order, maximum(horizons))

end
