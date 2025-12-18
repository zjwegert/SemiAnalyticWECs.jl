using Documenter
using SemiAnalyticWECs

makedocs(;
    sitename = "SemiAnalyticWECs.jl",
    format = Documenter.HTML(
      prettyurls = true,
      collapselevel = 1,
    ),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/zjwegert/SemiAnalyticWECs.jl/blob/{commit}{path}#L{line}",
)

deploydocs(;
    repo="github.com/zjwegert/SemiAnalyticWECs.jl",
)
