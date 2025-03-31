abstract type Model end

function is_fitted end

function require_fitted(model::M) where {M<:Model}
    is_fitted(model) && return true
    error("$(M) must first be estimated.")
end
