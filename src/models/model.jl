abstract type Model end

"""
    is_fitted(::Model)

Return `true` if the model as been estimated. Returns `false` otherwise.
"""
function is_fitted end

"""
    require_fitted(model::Model)

Returns `true` if the `model` is estimated. Throws an error otherwise.
"""
function require_fitted(model::M) where {M<:Model}
    is_fitted(model) && return true
    error("$(M) must first be estimated.")
end

"""
    coeffs(::Model, args...)

Returns the coefficient estimates of the model.
"""
function coeffs end

"""
    fitted(::Model)

Returns the fitted values of the of the estimated model. 
"""
function fitted end

"""
    residuals(::Model)

Returns the residuals of the estimated model.
"""
function residuals end

"""
    nobs(::Model)

Returns the effective number of observations used for the estimation of the
model.
"""
function nobs end

"""
    get_dependent(::Model)

Returns the dependent variables. 
"""
function get_dependent end

"""
    get_independent(::Model)

Returns the independent variables.
"""
function get_independent end

"""
    get_input_data(::Model)

Return the data used for the estimation.
"""
function get_input_data end

"""
    is_structural(::Model)

Returns `true` if the model is a structural model and `false` otherwise.
"""
function is_structural end

"""
    fit!(::Model, args...; kwargs...)

Estimates the model. 
"""
function fit! end

"""
    fit_and_select!(::Model, args..., selection_function::Function)

Estimates the model and selects among various models using `selection_function`. 
The best model has the smallest `selection_function` value, where 
`selection_function` must return a scalar.
"""
function fit_and_select! end

"""
    get_variable_names(model::Model) --> Vector{Symbol}
Obtain the variable names of the model's input data.
"""
function get_variable_names(model::Model)
    return Symbol.(names(get_input_data(model)))
end

"""
    through(
        model::Model, 
        variables::Union{AbstractVector{Symbol}, Symbol}, 
        horizons::Union{AbstractVector{<:Int},Vector{<:AbstractVector{<:Int}}}, 
        order::AbstractVector{Symbol}
    ) --> Q

Define a transmission channel through `variables` for `horizons`. Contrary to 
the other `through` methods, `variables` and `order` are defined using the 
variable names of the model -- other `through` methods use indices. 
"""
function through(
    model::Model, 
    variables::Union{AbstractVector{Symbol}, Symbol}, 
    horizons::Union{AbstractVector{<:Int},Vector{<:AbstractVector{<:Int}}}, 
    order::AbstractVector{Symbol}
)

    data = get_input_data(model)
    if isa(variables, AbstractVector)
        idx = _find_variable_idx.(variables, [data])
    else
        idx = _find_variable_idx(variables, data)
    end
    idx_order = _find_variable_idx.(order, [data])
    return through(idx, horizons, idx_order)
end

"""
    not_through(
        model::Model, 
        variables::Union{AbstractVector{Symbol}, Symbol}, 
        horizons::Union{AbstractVector{<:Int},Vector{<:AbstractVector{<:Int}}}, 
        order::AbstractVector{Symbol}
    ) --> Q

Define a transmission channel not through `variables` for `horizons`. Contrary to 
the other `not_through` methods, `variables` and `order` are defined using the 
variable names of the model -- other `not_through` methods use indices. 
"""
function not_through(
    model::Model, 
    variables::Union{AbstractVector{Symbol}, Symbol}, 
    horizons::Union{AbstractVector{<:Int},Vector{<:AbstractVector{<:Int}}}, 
    order::AbstractVector{Symbol}
)
    data = get_input_data(model)
    if isa(variables, AbstractVector)
        idx = _find_variable_idx.(variables, [data])
    else
        idx = _find_variable_idx(variables, data)
    end
    idx_order = _find_variable_idx.(order, [data])
    return not_through(idx, horizons, idx_order)
end

