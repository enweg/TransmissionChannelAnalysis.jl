using TransmissionMechanisms
using Documenter

DocMeta.setdocmeta!(TransmissionMechanisms, :DocTestSetup, :(using TransmissionMechanisms); recursive=true)

makedocs(;
    modules=[TransmissionMechanisms],
    authors="Enrico Wegner",
    repo="https://github.com/enweg/TransmissionMechanisms.jl/blob/{commit}{path}#{line}",
    sitename="TransmissionMechanisms.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://enweg.github.io/TransmissionMechanisms.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/enweg/TransmissionMechanisms.jl",
    devbranch="main",
)
