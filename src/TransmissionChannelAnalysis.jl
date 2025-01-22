module TransmissionChannelAnalysis

WEGNER = "Wegner et al (2024)"

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
export remove_contradictions
export make_condition
export to_transmission_irfs, transmission
export make_structural_B, make_structural_Omega, to_structural_transmission_model
export through, not_through

end
