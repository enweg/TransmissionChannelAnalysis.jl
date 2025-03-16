struct SVAR{I<:AbstractIdentificationMethod} <: Model

    # defining the DGP
    As::AbstractVector{AbstractMatrix{<:Number}}  # Coefficient matrices [A0, A1, ..., Ap]
    Phi0::AbstractMatrix{<:Number}                # Contemporaneous IRFs
    p::Int                                        # Lag order
    identification_method::I                      # Method used to identify SVAR 

    # additional data
    input_data::AbstractMatrix{<:Number}          # Data used to fit model
    var_names::AbstractVector{String}             # Variable names
end


