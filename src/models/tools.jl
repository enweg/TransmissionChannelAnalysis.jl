"""
    IRF(irfs::Array{<:Number,3}, varnames::Vector{Symbol},
        model::Model, ident_method::Union{Nothing, AbstractIdentificationMethod}=nothing)

Container for storing impulse response functions (IRFs) computed from a model.

This structure captures both reduced-form and structural IRFs. For IRFs
identified from a reduced-form model using an identification method
(e.g., recursive or external instruments), the `ident_method` should be
explicitly provided.

For models that are already structural (e.g., identified SVARs) or for
reduced-form IRFs (without identification), the `ident_method` should be left
as `nothing`.

# Fields
- `irfs::Array{<:Number,3}`: The IRFs with shape `(variables, shocks, horizons)`
- `varnames::Vector{Symbol}`: Variable names
- `model::Model`: The model object used to generate the IRFs
- `ident_method::Union{Nothing, AbstractIdentificationMethod}`:
  The identification method used, if applicable (only for RF models)

# Constructors
- `IRF(irfs, varnames, model)`
  Creates an IRF object for a structural model or reduced-form IRFs

- `IRF(irfs, varnames, model, ident_method)`
  Use this when IRFs are structurally identified from a reduced-form model
"""
struct IRF
    irfs::AbstractArray{<:Number}                              # contains the actual irfs 
    varnames::AbstractVector{Symbol}                           # variable names 
    model::Model                                               # model used to compute IRFs
    ident_method::Union{Nothing,AbstractIdentificationMethod}  # identification method if model is RF
end
function IRF(
    irfs::AbstractArray{<:Number},
    varnames::AbstractVector{Symbol},
    model::Model
)

    return IRF(irfs, varnames, model, nothing)
end

function Base.show(io::IO, mime::MIME"text/plain", x::IRF)
    return Base.show(io, mime, x.irfs)
end

Base.getindex(x::IRF, inds...) = x.irfs[inds...]
Base.:+(x::IRF, y) = x.irfs + y
Base.:+(x::IRF, y::Number) = x.irfs .+ y
Base.:+(x, y::IRF) = x + y.irfs
Base.:+(x::Number, y::IRF) = x .+ y.irfs
Base.:-(x::IRF, y::Number) = x.irfs .- y 
Base.:-(x::IRF, y) = x.irfs - y 
Base.:-(x, y::IRF) = x - y.irfs
Base.:-(x::Number, y::IRF) = x .- y.irfs
Base.:*(x::IRF, y) = x.irfs * y
Base.:*(x::IRF, y::Number) = x.irfs .* y
Base.:*(x, y::IRF) = x * y.irfs
Base.:*(x::Number, y::IRF) = x .* y.irfs
Base.:/(x::IRF, y) = x.irfs / y
Base.:/(x::IRF, y::Number) = x.irfs ./ y
Base.:/(x, y::IRF) = x / y.irfs
Base.:/(x::Number, y::IRF) = x ./ y.irfs
Base.size(x::IRF) = size(x.irfs)
Base.length(x::IRF) = length(x.irfs)

#-------------------------------------------------------------------------------
# API DESIGNS
#-------------------------------------------------------------------------------

"""
    IRF(model::Model, max_horizon::Int)

Returns the impulse response functions (IRFs) from `model` up to horizon 
`max_horizon`.
"""
function IRF(::M, ::Int) where {M<:Model}
    error("$(M) does not implement IRF.")
end

"""
    IRF(model::Model, method::AbstractIdentificationMethod, max_horizon::Int)

Returns the impulse response functions (IRFs) identified from `model`
using the identification method `method`, up to the specified
`max_horizon`.

Must return an object of type `IRF`. 
"""
function IRF(::M, ::I, max_horizon::Int) where {M<:Model,I<:AbstractIdentificationMethod}
    error("IRFs cannot be identified from $(M) using identification method $(I)")
end
