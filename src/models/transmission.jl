
function _inverse_order(order::AbstractVector{Int})
    inv = similar(order)
    for (i, val) in enumerate(order)
        inv[val] = i
    end
    return inv
end

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

function transmission(
    from::Int,
    model::SVAR,
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

