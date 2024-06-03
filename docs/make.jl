using TransmissionChannelAnalysis
using Documenter

DocMeta.setdocmeta!(TransmissionChannelAnalysis, :DocTestSetup, :(using TransmissionChannelAnalysis); recursive=true)

makedocs(;
    modules=[TransmissionChannelAnalysis],
    authors="Enrico Wegner",
    repo="https://github.com/enweg/TransmissionChannelAnalysis.jl/blob/{commit}{path}#{line}",
    sitename="TransmissionChannelAnalysis.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://enweg.github.io/TransmissionChannelAnalysis.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/enweg/TransmissionChannelAnalysis.jl",
    devbranch="main",
)
