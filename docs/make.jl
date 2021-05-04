using ActiveContours
using Documenter

DocMeta.setdocmeta!(ActiveContours, :DocTestSetup, :(using ActiveContours); recursive=true)

makedocs(;
    modules=[ActiveContours],
    authors="Dale <djblack@uci.edu> and contributors",
    repo="https://github.com/Dale-Black/ActiveContours.jl/blob/{commit}{path}#{line}",
    sitename="ActiveContours.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Dale-Black.github.io/ActiveContours.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Dale-Black/ActiveContours.jl",
)
