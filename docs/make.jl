push!(LOAD_PATH, "../src/")

using BellDiagonalQudits
using Documenter

makedocs(
    sitename="BellDiagonalQudits.jl",
    modules=[BellDiagonalQudits],
    pages=[
        "Home" => "index.md",
        "Manual" => "manual.md",
        "Library" => "library.md"
    ])

deploydocs(;
    repo="github.com/kungfugo/BellDiagonalQudits.jl"
)
