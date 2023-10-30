using TransmissionMechanisms
using Test
using Serialization
using JSON
using LinearAlgebra

@testset "TransmissionMechanisms.jl REMOVE_CONTRADICTIONS = true" begin
    # Write your tests here.
    include("./simplifying.jl")
    include("./transmission-function.jl")
end

@testset "TransmissionMechanisms.jl REMOVE_CONTRADICTIONS = false" begin
    # Write your tests here.
    TransmissionMechanisms.REMOVE_CONTRADICTIONS[] = false
    include("./simplifying.jl")
    include("./transmission-function.jl")
end