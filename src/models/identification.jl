using LinearAlgebra

abstract type AbstractIdentificationMethod end

#-------------------------------------------------------------------------------
# RECURSIVE IDENTIFICATION
# (Same as conditional ignorability assumption in micro)
#-------------------------------------------------------------------------------

struct Recursive <: AbstractIdentificationMethod end

function _identify(model::VAR, ::Recursive)
    is_fitted(model) || error("Reduced form model must first be estimated.")
    Sigma_u = cov(model)
    L = cholesky(Sigma_u).L
    A0 = inv(L)
    A_plus = A0 * coeffs(model)
    return A0, A_plus
end

struct InternalInstrument <: AbstractIdentificationMethod
    x::AbstractVector{<:Number}  # holds the instrument
    normalising_variable::Union{Symbol, Int}  
end

identify(model::VAR, method::InternalInstrument) = error("Internal instruments can only be used to identify IRFs and not the coefficients of the model.")
