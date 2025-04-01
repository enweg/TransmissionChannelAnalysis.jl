module TransmissionChannelAnalysis

WEGNER = "Wegner et al (2024)"
REMOVE_CONTRADICTIONS = Ref(true)

include("./utils.jl")
include("./simplifying.jl")
include("./condition.jl")
include("./transmission-function.jl")
include("./transmission-function-BOmega.jl")
include("./transmission-function-irfs.jl")
include("./structural-representation.jl")
include("./transmission-helper-functions.jl")

# export Q  # Not exported due to easy mistakes by users.
export show_y, @show_y
export map_x_to_y, map_y_to_x
export make_condition
export to_transmission_irfs, transmission
export make_B, make_Omega, make_systems_form
export through, not_through
export set_remove_contradictions

# Models
include("./models/utils.jl")
include("./models/model.jl")
include("./models/var.jl")
include("./models/identification.jl")
include("./models/tools.jl")
include("./models/svar.jl")

export Model, VAR
export coeffs, cov, fitted, residuals, nobs, get_dependent, get_independent, get_input_data
export is_fitted
export make_companion_matrix, spectral_radius, is_stable
export aic, hqc, sic, bic
export fit!, fit_and_select!
export simulate!, simulate

end
