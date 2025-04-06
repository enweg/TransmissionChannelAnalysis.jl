
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

function IRF(::M, ::Int) where {M<:Model}
    error("$(M) does not implement IRF.")
end

function IRF(::M, ::I, ::Int) where {M<:Model,I<:AbstractIdentificationMethod}
    error("$(M) does not implement IRF with identification method $(I).")
end
