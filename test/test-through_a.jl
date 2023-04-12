
@testset "through_a" begin
    rng = StableRNG(123)
    n = 3
    p = 2

    A = rand(rng, n, n)
    A = Matrix(UnitLowerTriangular(A))
    B = 0.2*randn(rng, n, n*p)
    b0 = randn(rng, n)
    Σ = diagm(rand(rng, Uniform(0, 10), n))

    svar = SVAR(n, p, FixedEstimated(A), FixedEstimated(B), FixedEstimated(b0), FixedEstimated(Σ))

    Ainv = inv(A)
    B = Ainv*B
    b0 = Ainv*b0 
    Σ = Ainv*Σ*Ainv'

    ndraws = 10
    nchains = 3
    B_bayes = cat([B for _ in 1:ndraws]...; dims=3)
    B_bayes = cat([B_bayes for _ in 1:nchains]...; dims=4)
    b0_bayes = cat([b0 for _ in 1:ndraws]...; dims=3)
    b0_bayes = cat([b0_bayes for _ in 1:nchains]...; dims=4)
    Σ_bayes = cat([Σ for _ in 1:ndraws]...; dims=3)
    Σ_bayes = cat([Σ_bayes for _ in 1:nchains]...; dims=4)


    ts = TSFrame(randn(100, n), Dates.today():Day(1):(Dates.today()+Day(99)))
    var = VAR(n, p, BayesianEstimated(B_bayes, nothing), BayesianEstimated(b0_bayes, nothing), BayesianEstimated(Σ_bayes, nothing), ts)
    sirfs = StructuralImpulseResponseFunction(var, 5, CholeskyVAR())
    sirfs = to_impact_normalisation(sirfs)

    sirfs_array = sirfs.irfs[:, :, :, 1, 1]
    mediating_effect = through_a(1, 3, 2, sirfs_array)
    path_effects = similar(mediating_effect)
    condition = x -> ((2,0) in x) || ((2,1) in x) || ((2,2) in x) || ((2,3) in x) || ((2,4) in x) || ((2,5) in x)
    for t in 0:5
        path_effects[t+1] = mediation(svar, (1,0), (3,t), condition)
    end
    @test all(isapprox.(mediating_effect - path_effects, 0, atol=1e-10))
    @info "Maximum distance: $(maximum(abs.(mediating_effect - path_effects)))"


    mediating_effect_bayes = through_a(1, 3, 2, sirfs)
    @test all(size(mediating_effect_bayes) .== (6, ndraws, nchains))
    @test all(mediating_effect_bayes.value[:, 1, 1] .== mediating_effect)
end

