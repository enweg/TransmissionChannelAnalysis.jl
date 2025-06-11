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
export make_companion_matrix, spectral_radius
include("./models/model.jl")
export Model, get_variable_names
include("./models/identification.jl")
export AbstractIdentificationMethod, Recursive
export InternalInstrument, ExternalInstrument
include("./models/tools.jl")
export IRF
include("./models/var/var.jl")
export VAR, coeffs, cov, fitted, residuals, nobs, get_dependent, get_independent
export get_input_data, is_structural, is_fitted, is_stable
export aic, hqc, sic, bic
export fit!, fit_and_select!
export simulate!, simulate
include("./models/svar/svar.jl")
export SVAR
export identify, identify!
include("./models/lp/lp.jl")
export LP
# include("./models/dfm/utils.jl")
# include("./models/dfm/abstract.jl")
# include("./models/dfm/pca.jl")
# include("./models/dfm/dynamic/dynamicDFM.jl")
# include("./models/dfm/static/staticDFM.jl")
# export PCA, DynamicDFM, StaticDFM, BaiNgPC1, BaiNgPC2, BaiNgPC3, BaiNgIC1,
#        BaiNgIC2, BaiNgIC3, select_factors!, get_factor_loadings, get_factors
include("models/static-dfm/pca.jl")
export PCA, BaiNgPC1, BaiNgPC2, BaiNgPC3, BaiNgIC1, BaiNgIC2, BaiNgIC3,
    select_factors!, factors, loadings
include("models/static-dfm/static_dfm.jl")
export DFM, get_factor_var, factors, loadings, is_scaled, is_centered
include("models/static-dfm/static_sdfm.jl")
export SDFM


include("./models/transmission.jl")

# Plotting
include("./plots.jl")
export plot_decomposition, plot_decomposition!, add_decomposition_legend!,
    plot_decomposition_comparison!, plot_decomposition_comparison,
    add_decompare_legend!


end
