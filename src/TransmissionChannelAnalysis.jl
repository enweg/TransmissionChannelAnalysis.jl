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

end
