module TransmissionMechanisms

WEGNER = "Wegner et al (2024)"

include("./simplifying.jl")
include("./condition.jl")
include("./transmission-function.jl")
include("./transmission-function-BQbb.jl")
include("./transmission-function-irfs.jl")
include("./structural-representation.jl")

# export Q  # Not exported due to easy mistakes by users.
export remove_contradictions
export make_condition
export to_transmission_irfs, transmission
export make_structural_B, make_structural_Qbb, to_structural_transmission_model

end
