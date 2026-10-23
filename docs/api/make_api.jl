# Extracts the API documentation of TransmissionChannelAnalysis.jl into the JSON
# file that the website (tca-website) renders with ApiDocs.jl.
#
# From the package root:
#   julia --project=docs/api docs/api/make_api.jl ../tca-website/documentation/_api/api-julia.json
#
# Source links point to the commit recorded in the JSON, so run this from a clean,
# pushed checkout (ideally a release tag).
using ApiDocs
using TransmissionChannelAnalysis
using CairoMakie                      # trigger the package extension

out = get(ARGS, 1, joinpath(@__DIR__, "api-julia.json"))

api = extract_julia(TransmissionChannelAnalysis;
    library = LibraryInfo(TransmissionChannelAnalysis;
        source_url = "https://github.com/enweg/TransmissionChannelAnalysis.jl/blob/{commit}/{path}#L{line}"))
resolve_refs!(api)
has_errors = report(stdout, lint(api))
write_api(out, api)
exit(has_errors ? 1 : 0)
