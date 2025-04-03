using DataFrames
using Test

@testset "LP construction" begin
    k = 4
    T = 10
    p = 2

    data = DataFrame(reshape(repeat(Float64.(1:T), k), :, k), :auto)

    treatment = 1
    horizons = 0:3
    model = LP(data, treatment, p, horizons; include_constant=true)

    # testing if constant was correctly defined
    @test all(model.X[:, 1] .== 1)
    # testing the rest of the matrices
    for (i, h) in enumerate(horizons)
        @test isequal(model.X[1:(end-h), 2], model.Y[1:(end-h), 2, i] .- h)
        for j=1:p
            @test isequal(model.X[1:(end-h), ((j-1)*k+3):(j*k+2)], model.Y[1:(end-h), :, i] .- h .- j)
        end
    end
end





