using Documenter
using Daisho

makedocs(
    sitename = "Daisho.jl",
    modules = [Daisho],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://mmbell.github.io/Daisho.jl",
    ),
    pages = [
        "Home" => "index.md",
        "Guide" => [
            "guide/radar_io.md",
            "guide/quality_control.md",
            "guide/gridding.md",
            "guide/srtm.md",
        ],
        "Theory" => [
            "theory/beam_geometry.md",
            "theory/gridding_algo.md",
        ],
        "API Reference" => "api.md",
    ],
    checkdocs = :none,
)
