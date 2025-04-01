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

function _identify_irfs(model::VAR, method::Recursive, max_horizon::Int)
    require_fitted(model)
    irfs = _var_irf(coeffs(model), model.p, max_horizon)
    A0, _ = _identify(model, method)
    Phi0 = inv(A0)
    for h = 0:max_horizon
        view(irfs, :, :, h + 1) .= view(irfs, :, :, h + 1) * Phi0
    end
    return irfs
end

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

_identify(::VAR, ::InternalInstrument) = error("Internal instruments can only be used to identify IRFs and not the coefficients of the model.")

function _identify_irfs(model::VAR, method::InternalInstrument, max_horizon::Int)
    require_fitted(model)
    data = get_input_data(model)
    idx_instrument = _find_variable_idx(method.instrument, data)
    idx_normalising = _find_variable_idx(method.normalising_variable, data)

    # internal instrument IRFs are cumputed as relative IRFs of the 
    # Cholesky shock of the instrument variable on the outcome variable 
    # over the Cholesky shock of the instrument variable on the normalising 
    # variable and the pre-defined horizon
    cholesky_irfs = _identify_irfs(model, Recursive(), max_horizon)
    normalising_factor = cholesky_irfs[idx_normalising, idx_instrument, method.normalising_horizon]
    cholesky_irfs ./= normalising_factor

    # we can only interpret the IRFs to the Cholesky shock of the instrument
    K = size(data, 2)
    T = eltype(cholesky_irfs)
    cholesky_irfs[:, filter(!=(idx_instrument), 1:K), :] .= T(NaN)

    return cholesky_irfs
end

