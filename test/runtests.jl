using TransmissionMechanisms
using MacroEconometrics
using Test
using StableRNGs
using LinearAlgebra
using Distributions
using TSFrames, Dates

@testset "TransmissionMechanisms.jl" begin
    include("test-paths.jl")
    include("test-through_a.jl")
end