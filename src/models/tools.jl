
struct IRF
    irfs::AbstractArray{<:Number}                              # contains the actual data 
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

#-------------------------------------------------------------------------------
# API DESIGNS
#-------------------------------------------------------------------------------

function IRF(::M, ::Int) where {M<:Model}
    error("$(M) does not implement IRF.")
end

function IRF(::M, ::I, ::Int) where {M<:Model,I<:AbstractIdentificationMethod}
    error("$(M) does not implement IRF with identification method $(I).")
end
