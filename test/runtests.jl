using TransmissionChannelAnalysis
using Test
using Serialization
using JSON
using LinearAlgebra

# common functions used for tests in "transmission-function.jl"
include("./transmission-function-tests.jl")

@testset "Utils" begin
    include("./utils.jl")
end

@testset "TransmissionChannelAnalysis.jl REMOVE_CONTRADICTIONS = true" begin
    # Write your tests here.
    include("./simplifying.jl")
    include("./transmission-function.jl")
end

@testset "TransmissionChannelAnalysis.jl REMOVE_CONTRADICTIONS = false" begin
    # Write your tests here.
    TransmissionChannelAnalysis.REMOVE_CONTRADICTIONS[] = false
    include("./simplifying.jl")
    include("./transmission-function.jl")
end
