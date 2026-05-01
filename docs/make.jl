using Pkg

Pkg.develop(PackageSpec(path = dirname(@__DIR__)))

using ACMG
using Documenter

DocMeta.setdocmeta!(ACMG, :DocTestSetup, :(using ACMG); recursive = true)

makedocs(;
    modules = [ACMG],
    sitename = "ACMG.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        size_threshold = 300 * 1024,
        size_threshold_warn = 250 * 1024,
    ),
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "Concepts" => "concepts.md",
        "Modular Data" => "modular_data.md",
        "F/R Symbols" => "fr_symbols.md",
        "Gauge Fixing" => "gauge_fixing.md",
        "Braid Representations" => "braid_representations.md",
        "Finite-Field FRData" => "finite_field_fr_data.md",
        "Finite Fields" => "finite_fields.md",
        "Block-U Search" => "search_blocku.md",
        "Zariski Diagnostics" => "zariski_diagnostics.md",
        "API Stability" => "api_stability.md",
        "API Reference" => "api.md",
    ],
    warnonly = true,
)

deploydocs(;
    repo = "github.com/MGYamada/ACMG.jl.git",
)
