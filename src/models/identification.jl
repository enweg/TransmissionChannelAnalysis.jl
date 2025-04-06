using LinearAlgebra

abstract type AbstractIdentificationMethod end

function _identify(::M, ::I) where {M<:Model,I<:AbstractIdentificationMethod}
    error("_identify is not implemented for model $M and method $I.")
end

function _identify_irfs(::M, ::I) where {M<:Model,I<:AbstractIdentificationMethod}
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

function Base.show(io::IO, ::MIME"text/plain", x::InternalInstrument)
    s = """
        Internal Instrument Identificaton Method 
        ========================================
        Instrument: $(x.instrument)
        Unit effect normalisation horizon: $(x.normalising_horizon)
        Unit effect normalising variable: $(x.normalising_variable)
        """
    println(io, s)
end

#-------------------------------------------------------------------------------
# EXTERNAL INTSTRUMENT IDENTIFICATION
#-------------------------------------------------------------------------------

struct ExternalInstrument <: AbstractIdentificationMethod
    instruments::Union{AbstractVector{<:Symbol},AbstractVector{<:Int}}
    normalising_horizon::Int
end
function ExternalInstrument(
    instruments::Union{Symbol,Int,AbstractVector{<:Symbol},AbstractVector{<:Int}};
    normalising_horizon::Int=0
)

    if !isa(instruments, AbstractVector)
        instruments = [instruments]
    end

    return ExternalInstrument(instruments, normalising_horizon)
end

function Base.show(io::IO, ::MIME"text/plain", x::ExternalInstrument)
    s = """
        External Instrument Identificaton Method 
        ========================================
        Instruments: $(x.instruments)
        Unit effect normalisation horizon: $(x.normalising_horizon)
        """
    println(io, s)
end
