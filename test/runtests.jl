using TransmissionChannelAnalysis
using Test
using Serialization
using JSON
using LinearAlgebra

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