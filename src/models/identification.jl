using LinearAlgebra

abstract type AbstractIdentificationMethod end

function _identify(::M, ::I) where {M<:Model, I<:AbstractIdentificationMethod} 
    error("_identify is not implemented for model $M and method $I.")
end

function _identify_irfs(::M, ::I) where {M<:Model, I<:AbstractIdentificationMethod}
    error("_identify_irfs is not implemented for model $M and method $I.")
end

#-------------------------------------------------------------------------------
# RECURSIVE IDENTIFICATION
# (Same as conditional ignorability assumption in micro)
#-------------------------------------------------------------------------------

struct Recursive <: AbstractIdentificationMethod end

#-------------------------------------------------------------------------------
# INTERNAL INSTRUMENT IDENTIFICATION
# (Can only handle a single instrument)
#-------------------------------------------------------------------------------

struct InternalInstrument <: AbstractIdentificationMethod
    instrument::Union{Symbol,Int}
    normalising_variable::Union{Symbol,Int}
    normalising_horizon::Int
end
function InternalInstrument(
    normalising_variable::Union{Symbol,Int};
    instrument::Union{Symbol,Int}=1,
    normalising_horizon::Int=0
)

    return InternalInstrument(instrument, normalising_variable, normalising_horizon)
end

#-------------------------------------------------------------------------------
# EXTERNAL INTSTRUMENT IDENTIFICATION
#-------------------------------------------------------------------------------

# TODO: allow for other normalisation horizons than contemporaneous 
# this is better for news shocks
# - technically speaking this is already possible. The user just 
# needs to manually lead the treatment variable before providing 
# it to LP
struct ExternalInstrument <: AbstractIdentificationMethod
    instruments::Union{AbstractVector{<:Symbol}, AbstractVector{<:Int}}
end
