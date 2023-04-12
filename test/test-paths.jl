@testset "find_paths" begin
    B = [
        1 0 0;
        1 1 1; 
        0 1 1
    ]
    A = [
        1 0 0;
        1 1 0; 
        0 1 1 
    ]

    paths = find_paths(A, B, (1, 0), (3, 2))
    paths_compare = [[(1, 0), (2, 0), (3, 0), (3, 1), (3, 2)], [(1, 0), (2, 0), (3, 1), (3, 2)], [(1, 0), (2, 0), (3, 0), (2, 1), (3, 1), (3, 2)], [(1, 0), (2, 0), (2, 1), (3, 1), (3, 2)], [(1, 0), (2, 1), (3, 1), (3, 2)], [(1, 0), (1, 1), (2, 1), (3, 1), (3, 2)], [(1, 0), (2, 0), (3, 0), (2, 1), (3, 2)], [(1, 0), (2, 0), (2, 1), (3, 2)], [(1, 0), (2, 1), (3, 2)], [(1, 0), (1, 1), (2, 1), (3, 2)], [(1, 0), (2, 0), (3, 0), (3, 1), (2, 2), (3, 2)], [(1, 0), (2, 0), (3, 1), (2, 2), (3, 2)], [(1, 0), (2, 0), (3, 0), (2, 1), (3, 1), (2, 2), (3, 2)], [(1, 0), (2, 0), (2, 1), (3, 1), (2, 2), (3, 2)], [(1, 0), (2, 1), (3, 1), (2, 2), (3, 2)], [(1, 0), (1, 1), (2, 1), (3, 1), (2, 2), (3, 2)], [(1, 0), (2, 0), (3, 0), (2, 1), (2, 2), (3, 2)], [(1, 0), (2, 0), (2, 1), (2, 2), (3, 2)], [(1, 0), (2, 1), (2, 2), (3, 2)], [(1, 0), (1, 1), (2, 1), (2, 2), (3, 2)], [(1, 0), (1, 1), (2, 2), (3, 2)], [(1, 0), (1, 1), (1, 2), (2, 2), (3, 2)]]
    @test length(paths) == length(paths_compare)
    @test all([x in paths_compare for x in paths])
    svar = SVAR(3, 1, FixedEstimated(Float64.(A)), FixedEstimated(Float64.(B)), FixedEstimated(randn(3)), FixedEstimated(randn(3, 3)))
    paths_svar = find_paths(svar, (1, 0), (3, 2))
    @test length(paths) == length(paths_svar)
    @test all([x in paths_svar for x in paths])

    paths = find_paths(A, B, (2, 3), (1, 4))
    # None should exist
    @test length(paths) == 0

    paths = find_paths(A, B, (1, 3), (2, 3))
    paths_compare = [[(1, 3), (2, 3)]]
    @test length(paths) == length(paths_compare)
    @test all([x in paths_compare for x in paths])

    paths = find_paths(A, B, (1, 5), (3, 2))
    @test length(paths) == 0


end

@testset "calculate_path_effects" begin
    B = [
        0.5 0 0; 
        0.5 0.4 0.3; 
        0 0.5 0.4
    ]
    A = [
        1.5 0 0; 
        1.5 1.0 0; 
        0 1.5 1.0;
    ]


    path_effect = calculate_path_effect(A, B, [(1,0), (2,0), (2, 1), (3, 1)])
    path_effect_manual = 1.5*0.4*1.5
    @test isapprox(path_effect, path_effect_manual, atol=1e-10)

    ndraws = 10
    B_bayes = cat([B for _ in 1:ndraws]...; dims=3)
    A_bayes = cat([A for _ in 1:ndraws]...; dims=3)
    nchains = 2
    B_bayes = cat([B_bayes for _ in 1:nchains]...; dims=4)
    A_bayes = cat([A_bayes for _ in 1:nchains]...; dims=4)

    svar = SVAR(3, 1, BayesianEstimated(A_bayes, nothing), BayesianEstimated(B_bayes, nothing), BayesianEstimated(randn(3), nothing), BayesianEstimated(randn(3, 3), nothing))
    path_effect_bayes = calculate_path_effect(svar, [(1,0), (2,0), (2,1), (3,1)])
    @test allequal(path_effect_bayes.value)
    @test isapprox(path_effect_bayes[1, 1], path_effect_manual, atol=1e-10)
    @test all(size(path_effect_bayes) .== (ndraws, nchains))
end

@testset "mediation" begin
    B = [
        0.5 0 0; 
        0.5 0.4 0.3; 
        0 0.5 0.4
    ]
    A = [
        1.5 0 0; 
        1.5 1.0 0; 
        0 1.5 1.0;
    ]
    n = 3
    p = 1

    paths = [
        [(1,0), (2,0), (2,1), (3,1)], 
        [(1,0), (1,1), (2,1), (3,1)] 
    ]

    svar = SVAR(n, p, FixedEstimated(A), FixedEstimated(B), FixedEstimated(randn(n)), FixedEstimated(randn(n, n)))
    condition = x -> true
    total_mediation = mediation(svar, paths, condition)
    total_mediation_manual = calculate_path_effect(svar, paths[1]) + calculate_path_effect(svar, paths[2])
    @test isapprox(total_mediation, total_mediation_manual, atol=1e-10)
    filtered_mediation = mediation(svar, paths, x -> (1,1) in x)
    filtered_mediation_manual = calculate_path_effect(svar, paths[2])
    @test isapprox(filtered_mediation, filtered_mediation_manual, atol=1e-10)


    ndraws = 10
    B_bayes = cat([B for _ in 1:ndraws]...; dims=3)
    A_bayes = cat([A for _ in 1:ndraws]...; dims=3)
    nchains = 2
    B_bayes = cat([B_bayes for _ in 1:nchains]...; dims=4)
    A_bayes = cat([A_bayes for _ in 1:nchains]...; dims=4)
    svar = SVAR(3, 1, BayesianEstimated(A_bayes, nothing), BayesianEstimated(B_bayes, nothing), BayesianEstimated(randn(3), nothing), BayesianEstimated(randn(3, 3), nothing))
    condition = x -> true
    total_mediation = mediation(svar, paths, condition)
    @test all(size(total_mediation) .== (ndraws, nchains))
    @test allequal(total_mediation.value)
    @test isapprox(total_mediation[1, 1], total_mediation_manual, atol=1e-10)
    filtered_mediation = mediation(svar, paths, x -> (1,1) in x)
    @test all(size(filtered_mediation) .== (ndraws, nchains))
    @test allequal(filtered_mediation.value)
    @test isapprox(filtered_mediation[1, 1], filtered_mediation_manual, atol=1e-10)
end