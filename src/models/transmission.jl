"""
    _inverse_order(order::AbstractVector{<:Int}) --> AbstractVector{<:Int}

If `order` defines an ordering, then `_inverse_order(order)` is the inverse 
order. 
"""
function _inverse_order(order::AbstractVector{Int})
    inv = similar(order)
    for (i, val) in enumerate(order)
        inv[val] = i
    end
    return inv
end

"""
    _irf_vec_to_array(irfs::AbstractVector, 
                      k::Int, 
                      order::AbstractVector{<:Int}) --> Array{<:Number, 3}

Obtain IRFs in the standard three-dimensional array representation from IRFs 
in transmission representation, i.e. represented as a matrix in which the 
third dimension of the standard representation is stacked above each other. 
If `order` is not `1:k`, then the endogenous variables will be re-ordered 
using the inverse of `order`. 

## Arguments
- `irfs::AbstractVector`: Vector of IRFs of a single shock with dimension 
  (k*h, 1). In transmission representation, i.e. horizons are stacked above 
  each other. 
- `k::Int`: Number of endogenous variables. 
- `order::AbstractVector{<:Int}`: Ordering defined by the transmission matrix. 
"""
function _irf_vec_to_array(
    irfs::AbstractVector,
    k::Int,
    order::AbstractVector{<:Int}=1:k
)

    horizons = div(length(irfs), k) - 1
    irfs_array = zeros(eltype(irfs), k, 1, horizons + 1)
    inv_order = _inverse_order(order)

    for h = 0:horizons
        irfs_array[:, :, h+1] .= (irfs[(h*k+1):((h+1)*k)])[inv_order]
    end
    return irfs_array
end

"""
    transmission(model::Model, 
                 from::Int, 
                 q::Q, 
                 order::AbstractVector{<:Int}, 
                 max_horizon::Int) --> Array{<:Number, 3}

    transmission(model::Model, 
                 method::AbstractIdentificationMethod, 
                 from::Int, 
                 q::Q,
                 order::AbstractVector{<:Int}, 
                 max_horizon::Int) --> Array{<:Number, 3}

Compute the transmission effect of a transmission channel defined by the 
condition `q`. If `model` is a reduced-form model, `method` will be used to 
identify the required structural shock. 

## Arguments
- `model::Model`: A model, such as an `SVAR`, `VAR`, or `LP`. 
- `method::AbstractIdentificationMethod`: An identification method to identify 
    the `from`-th structural shock. 
- `from::Int`: Shock number. 
- `q::Q`: A transmission condition. See also `Q` and `make_condition`. 
- `order::AbnstractVector{<:Int}`: order of variables defined by the transmission 
  matrix. 
- `max_horizon::Int`: Maximum horizon for the transmission effect.  

## Returns
- Returns a three dimensional array with the first dimension correspondin 
  to the endogenous variables (in original order), the second to the shock, 
  and the third to the horizon (from 0 to `max_horizon`).
"""
function transmission(
    model::SVAR,
    from::Int,
    q::Q,
    order::AbstractVector{<:Int},
    max_horizon::Int
)

    require_fitted(model)

    k = size(get_input_data(model), 2)
    k == length(order) || error("length(order) != k")

    A0, _ = coeffs(model, true)
    Phi0 = inv(A0)
    As = coeffs(model.var, true)
    As = _separate_lag_matrices(As, model.p)
    Sigma = cov(model.var)
    Psis = Vector{Matrix}(undef, 0)

    B, Omega = make_systems_form(Phi0, As, Psis, Sigma, order, max_horizon)
    transmission_effects = transmission(from, B, Omega, q; method=:BOmega)
    return _irf_vec_to_array(transmission_effects, k, order)
end

function transmission(
    model::VAR,
    method::Recursive,
    from::Int,
    q::Q,
    order::AbstractVector{<:Int},
    max_horizon::Int
)

    require_fitted(model)
    k = size(get_input_data(model), 2)
    k == length(order) || error("length(order) != k")

    model_svar = identify(model, method)
    return transmission(model_svar, from, q, order, max_horizon)
end

function transmission(
    model::VAR,
    method::InternalInstrument,
    from::Int,
    q::Q,
    order::AbstractVector{<:Int},
    max_horizon::Int
)


    require_fitted(model)
    k = size(get_input_data(model), 2)
    k == length(order) || error("length(order) != k")
    from == 1 || error("Internal Instruments only identify a single shock. `from` must thus be equal to `1`.")

    # TODO: should we remove the instrument from the irfs? 
    irfs = _identify_irfs(model, method, max_horizon)
    # re-ordering variables
    irfs = irfs[order, :, :]
    # TODO: this is not the most efficient way. We could also just re-order 
    # the already estimated coefficients and covariance matrix. However, 
    # this is the most general and most robust to making mistakes.
    data = get_input_data(model)
    model_tmp = VAR(data[:, order], model.p; trend_exponents=model.trend_exponents)
    fit!(model_tmp)
    irfs_ortho = _identify_irfs(model_tmp, Recursive(), max_horizon)

    irfs = to_transmission_irfs(irfs)
    irfs_ortho = to_transmission_irfs(irfs_ortho)

    transmission_effects = transmission(from, irfs, irfs_ortho, q; method=:irfs)
    return _irf_vec_to_array(transmission_effects, k, order)
end

function transmission(
    model::LP,
    method::Union{Recursive,ExternalInstrument},
    from::Int,
    q::Q,
    order::AbstractVector{<:Int},
    max_horizon::Int
)

    data = get_input_data(model)
    k = size(data, 2)

    irfs = _identify_irfs(model, method, max_horizon)
    # re-ordering according to transmission ordering
    irfs = irfs[order, :, :]
    # Getting orthogonalised IRFs via a recursive identification method 
    # using data re-ordered according to the transmission ordering
    irfs_ortho = zeros(eltype(irfs), k, k, max_horizon + 1)
    for treatment = 1:k
        model_tmp = LP(data[:, order], treatment, model.p, 0:max_horizon; include_constant=model.include_constant)
        irfs_ortho[:, treatment:treatment, :] .= _identify_irfs(model_tmp, Recursive(), max_horizon)[:, treatment:treatment, :]
    end

    irfs = to_transmission_irfs(irfs)
    irfs_ortho = to_transmission_irfs(irfs_ortho)

    transmission_effects = transmission(from, irfs, irfs_ortho, q; method=:irfs)
    return _irf_vec_to_array(transmission_effects, k, order)
end
