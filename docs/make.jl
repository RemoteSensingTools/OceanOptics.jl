using Documenter
using OceanOptics

# Interactive charts: `@example` blocks call `Main.DocCharts.linechart`,
# rendered in the browser by docs/src/assets/oo_charts.js (Chart.js).
include(joinpath(@__DIR__, "DocCharts.jl"))

DocMeta.setdocmeta!(OceanOptics, :DocTestSetup, :(using OceanOptics); recursive = true)

makedocs(
    sitename = "OceanOptics.jl",
    authors  = "RemoteSensingTools",
    modules  = [OceanOptics],
    format   = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical  = "https://RemoteSensingTools.github.io/OceanOptics.jl/",
        assets     = ["assets/oo_charts.js"],
        # Chart pages embed their data (JSON + table) and have no image
        # fallback, so allow larger pages and example outputs.
        size_threshold         = 400 * 2^10,
        size_threshold_warn    = 250 * 2^10,
        example_size_threshold = 96 * 2^10,
    ),
    pages    = [
        "Home"                 => "index.md",
        "End-to-end example"   => "example.md",
        "Gallery"              => "gallery.md",
        "Comparison with HydroLight" => "hydrolight_comparison.md",
        "Consumer API"         => "consumer_api.md",
        "vSmartMOM integration" => "vsmartmom_integration.md",
        "Reference"            => [
            "Types"            => "reference/types.md",
            "Materials"        => "reference/materials.md",
            "Phase functions"  => "reference/phases.md",
            "Inelastic"        => "reference/inelastic.md",
            "Data loading"     => "reference/io.md",
        ],
    ],
    warnonly   = [:missing_docs, :cross_references],
    checkdocs  = :exports,
)

deploydocs(
    repo      = "github.com/RemoteSensingTools/OceanOptics.jl.git",
    devbranch = "main",
    push_preview = true,
)
