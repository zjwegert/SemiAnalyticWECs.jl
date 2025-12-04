using Documenter
using SemiAnalyticWECs

makedocs(;
    modules=[GridapPardiso],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/zjwegert/SemiAnalyticWECs.jl/blob/{commit}{path}#L{line}",
    sitename="SemiAnalyticWECs.jl",
    authors="Zachary J Wegert <wegert@qut.edu.au>",
    assets=String[],
)

deploydocs(;
    repo="github.com/zjwegert/SemiAnalyticWECs.jl",
)