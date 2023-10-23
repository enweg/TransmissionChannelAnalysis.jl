using TransmissionMechanisms
using Test
using Serialization
using LinearAlgebra

@testset "TransmissionMechanisms.jl" begin
    # Write your tests here.
    include("./simplifying.jl")
    include("./transmission-function.jl")
end
