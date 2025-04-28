using LinearAlgebra

"""
    AbstractIdentificationMethod

Abstract supertype for all identification methods.

New identification methods must subtype `AbstractIdentificationMethod`.
If the identification method can be used to estimate the structural model, i.e. 
identify the structural model's coefficients, then functions `fit!` and 
`fit_and_select!` should be implemented. If the identification method can be 
used to identify structural IRFs from a reduced form model, e.g. structural 
IRFs from a VAR, then the method `IRF` should be implemented. For further 
details regarding the latter point, see the documentation for `IRF`. 

For structural VARs (SVARs), developers may alternatively extend the
internal methods `_identify` and `_identify_irfs`.
"""
abstract type AbstractIdentificationMethod end

"""
    fit!(model::Model, method::AbstractIdentificationMethod)

Fits the model `model` using the identification method `method`. All data 
needed for the estimation should be contained in the `model` struct.
"""
function fit!(::M, ::I) where {M<:Model,I<:AbstractIdentificationMethod}
    error("$(M) does not implement identification method $(I)")
end

"""
    fit_and_select!(model::Model, method::AbstractIdentificationMethod, selection_func::Function)

Fits the model `model` using the identification method `method` and
selects the best specification based on the `selection_func`.

The selection function should return a scalar value, with lower values
indicating better model fit. Examples include information criteria such
as AIC.
"""
function fit_and_select!(::M, ::I, selection_func::Function) where {M<:Model,I<:AbstractIdentificationMethod}
    error("$(M) does not implement identification method $(I)")
end


#-------------------------------------------------------------------------------
# RECURSIVE IDENTIFICATION
# (Same as conditional ignorability assumption in micro)
#-------------------------------------------------------------------------------

"""
    Recursive <: AbstractIdentificationMethod

General-purpose identification method based on recursive, or
conditional ignorability, assumptions.

This method assumes that, conditional on variables ordered before, a
given variable (often interpreted as a treatment or shock) is as good
as random. This enables causal interpretation of the identified shock
as structural.

In macroeconomic applications, this is typically operationalized through
a Cholesky decomposition in structural VARs (SVARs), but the method is
not limited to this setting. It applies to any model in which recursive
ordering (conditioning on other variables) can be used to justify identification.

Conditioning variables are those that are ordered before the treatment variable
in the provided dataset.
"""
struct Recursive <: AbstractIdentificationMethod end

#-------------------------------------------------------------------------------
# INTERNAL INSTRUMENT IDENTIFICATION
# (Can only handle a single instrument)
#-------------------------------------------------------------------------------

"""
    InternalInstrument <: AbstractIdentificationMethod

Identification method using internal instruments within a VAR framework.

This method estimates relative structural impulse response functions (IRFs)
using an internal instrument—i.e., a variable within the system that is
assumed to proxy for a structural shock of interest. The IRFs are normalised
such that the response of a designated normalising variable equals one at
a specific horizon.

This approach is based on:

Plagborg-Møller, M., & Wolf, C. K. (2021). Local Projections and VARs Estimate
    the Same Impulse Responses. *Econometrica*, 89(2), 955–980.
    https://doi.org/10.3982/ecta17813

# Fields
- `instrument::Union{Symbol, Int}`: the variable to be used as the instrument
- `normalising_variable::Union{Symbol, Int}`: the variable used to normalise
  the IRF
- `normalising_horizon::Int`: the horizon at which the IRF of the normalising
  variable is set to one
"""
struct InternalInstrument <: AbstractIdentificationMethod
    instrument::Union{Symbol,Int}
    normalising_variable::Union{Symbol,Int}
    normalising_horizon::Int
end
"""
    InternalInstrument(normalising_variable::Union{Symbol, Int};
                       instrument::Union{Symbol, Int}=1,
                       normalising_horizon::Int=0)

Constructs an `InternalInstrument` identification method with the given
normalising variable, optional instrument (default is the first variable),
and normalisation horizon (default is 0).
"""
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
"""
    ExternalInstrument <: AbstractIdentificationMethod

Identification method using external instruments for structural shock
identification in time series models.

External instruments serve as proxies for structural shocks and are not
included in the estimation model (e.g., a VAR), but are used in
instrumental variable estimators such as two-stage least squares (2SLS).

This method identifies *relative* structural impulse response functions
(IRFs), meaning that the scale of the identified shock is normalised.
Specifically, the IRF is scaled such that the response of the treatment
variable equals one at a given normalising horizon.

This approach also generalises to the Local Projections Instrumental
Variables (LP-IV) framework and supports the use of multiple instruments.

The method is based on:

Stock, J. H., & Watson, M. W. (2018). Identification and Estimation of
Dynamic Causal Effects in Macroeconomics Using External Instruments.
*The Economic Journal*, 128(610), 917–948.
https://doi.org/10.1111/ecoj.12593

# Fields
- `treatment::Union{Symbol, Int}`: the variable used to normalise the IRF
- `instruments::AbstractMatrix`: matrix of external instruments
- `normalising_horizon::Int`: the horizon at which the IRF of the treatment
  variable is set to one

"""
struct ExternalInstrument <: AbstractIdentificationMethod
    treatment::Union{Symbol, Int}
    instruments::AbstractMatrix
    normalising_horizon::Int
end

"""
    ExternalInstrument(treatment::Union{Symbol, Int},
                       instruments::Union{AbstractMatrix, DataFrame};
                       normalising_horizon::Int=0)

Constructs an `ExternalInstrument` identification method using a specified
`treatment` variable, one or more external instruments, and an optional
normalisation horizon (default is 0).

## Arguments 
- `treatment::Union{Symbol, Int}`: treatment variable in the model data specified 
  either as the column index in the data or by variable name. 
- `instruments::Union{AbstractMatrix, DataFrame}`: observations of external 
  instruments. Instruments must be observed over the same period as the model data. 

## Keyword Arguments 
- `normalising_horizon::Int`: horizon of the treatment variable with respect to 
  which all IRFs are normalised. 
"""
function ExternalInstrument(
    treatment::Union{Symbol, Int},
    instruments::Union{AbstractMatrix, DataFrame};
    normalising_horizon::Int=0
)

    if isa(instruments, DataFrame)
        instruments = Matrix(instruments)
    end

    return ExternalInstrument(treatment, instruments, normalising_horizon)
end

function Base.show(io::IO, ::MIME"text/plain", x::ExternalInstrument)
    s = """
        External Instrument Identificaton Method 
        ========================================
        Treatment: $(x.treatment)
        Number of Instruments: $(size(x.instruments, 2))
        Unit effect normalisation horizon: $(x.normalising_horizon)
        """
    println(io, s)
end
