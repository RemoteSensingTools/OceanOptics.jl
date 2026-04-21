using Documenter
using OceanOptics

# GR backend needs an off-screen target in headless CI
ENV["GKSwstype"] = "100"

DocMeta.setdocmeta!(OceanOptics, :DocTestSetup, :(using OceanOptics); recursive = true)

makedocs(
    sitename = "OceanOptics.jl",
    authors  = "RemoteSensingTools",
    modules  = [OceanOptics],
    format   = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical  = "https://RemoteSensingTools.github.io/OceanOptics.jl/",
        assets     = String[],
    ),
    pages    = [
        "Home"                 => "index.md",
        "End-to-end example"   => "example.md",
        "Gallery"              => "gallery.md",
        "Consumer API"         => "consumer_api.md",
        "vSmartMOM integration audit" => "vsmartmom_integration.md",
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
