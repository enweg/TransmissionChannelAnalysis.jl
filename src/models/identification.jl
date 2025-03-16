abstract type AbstractIdentificationMethod end

struct Recursive <: AbstractIdentificationMethod end

struct InternalInstrument <: AbstractIdentificationMethod 
  x::AbstractVector{<:Number}  # holds the instrument
end

struct ExternalInstrument <: AbstractIdentificationMethod
  x::AbstractVecOrMat{<:Number}  # holds the instrument(s)
end
