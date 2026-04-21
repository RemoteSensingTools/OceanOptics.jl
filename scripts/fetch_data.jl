#!/usr/bin/env julia
#
# scripts/fetch_data.jl
#
# Regenerate bundled reference-data CSVs in `data/` from their original
# upstream sources. Run this script when:
#
#   • the package is first cloned (optional — the CSVs are committed),
#   • an upstream table has been corrected / updated, or
#   • you want to verify that the committed CSVs match the upstream files.
#
# Usage:  julia --project=.. scripts/fetch_data.jl [--verify]
#
# With `--verify`, the script downloads each file, regenerates the CSV
# in-memory, and diffs against the committed version. Exits non-zero on
# mismatch. Without the flag, it overwrites `data/*.csv` from upstream.
#
# Datasets
# --------
#   pope_fry_1997.csv           ← omlc.org/spectra/water/data/pope97.txt
#   smith_baker_1981.csv        ← omlc.org/spectra/water/data/smith81.txt
#
# Tables *not* fetched by this script (no stable redistributable URL; must
# be transcribed from the original publication):
#   mason_cone_fry_2016.csv     — Appl. Opt. 55, 7163, Table 2
#   bricaud_1995.csv            — JGR 100, 13321, Table 2
#   petzold_sdh.csv             — Mobley (1994) "Light and Water" Appendix A.4
#
# Unit conventions
# ----------------
# Every CSV stores absorption in `[1/m]`. OMLC publishes in `[1/cm]`;
# fetch_data.jl applies the factor-100 conversion at write time and
# records it in the CSV provenance header.

using Dates
using Downloads
using Printf

const DATA_DIR = abspath(joinpath(@__DIR__, "..", "data"))
const UPSTREAM = Dict(
    "pope_fry_1997"    => (
        url          = "https://omlc.org/spectra/water/data/pope97.txt",
        citation     = "R. M. Pope and E. S. Fry (1997), \"Absorption spectrum (380-700 nm) of pure water. II. Integrating cavity measurements,\" Appl. Opt. 36, 8710-8723.",
        doi          = "10.1364/AO.36.008710",
        dataset_desc = "Pure-water absorption coefficient (liquid water)",
        grid_desc    = "380.0 to 727.5 nm, 2.5 nm resolution",
    ),
    "smith_baker_1981" => (
        url          = "https://omlc.org/spectra/water/data/smith81.txt",
        citation     = "R. C. Smith and K. S. Baker (1981), \"Optical properties of the clearest natural waters (200-800 nm),\" Appl. Opt. 20, 177-184.",
        doi          = "10.1364/AO.20.000177",
        dataset_desc = "Pure-water absorption coefficient (\"clearest natural waters\")",
        grid_desc    = "200 to 800 nm, 10 nm resolution",
    ),
)

"""
Parse the OMLC plain-text format. Files have ~3 lines of free-text
citation, a `lambda / absorption` column-name line, a units line, and
then rows of `λ[nm]\\t a[1/cm]`. We look for the first row whose first
token parses as a float; everything after that is data.
"""
function parse_omlc(bytes::Vector{UInt8})
    # Upstream encodes citation text with ISO-8859-1 characters (en-dash);
    # decode permissively by treating non-ASCII bytes as-is (the citation
    # is discarded anyway — we only care about the numeric rows).
    text = String(filter(b -> b < 0x80, bytes))
    λs   = Float64[]
    ys   = Float64[]
    started = false
    for raw in eachline(IOBuffer(text))
        line = strip(raw)
        isempty(line) && continue
        parts = split(line)
        length(parts) < 2 && continue
        λ = tryparse(Float64, parts[1])
        y = tryparse(Float64, parts[2])
        if λ === nothing || y === nothing
            started && break
            continue
        end
        started = true
        push!(λs, λ)
        push!(ys, y)
    end
    return (λ_nm = λs, absorption_inv_cm = ys)
end

"Format the value with 7 significant digits, plain decimal notation."
_fmt(x) = @sprintf("%.7g", x)

function write_absorption_csv(path, meta, λs, as_inv_m)
    today = Dates.format(Dates.today(), dateformat"yyyy-mm-dd")
    open(path, "w") do io
        println(io, "# Dataset:      ", meta.dataset_desc)
        println(io, "# Source:       ", meta.citation)
        println(io, "# DOI:          ", meta.doi)
        println(io, "# Retrieved:    ", meta.url, " (", today, ")")
        println(io, "# Range:        ", meta.grid_desc)
        println(io, "# Unit change:  Source 1/cm -> stored 1/m (factor 100)")
        println(io, "# Columns:      lambda_nm, absorption_m_inv")
        println(io, "lambda_nm,absorption_m_inv")
        for (λ, a) in zip(λs, as_inv_m)
            println(io, @sprintf("%.2f", λ), ",", _fmt(a))
        end
    end
    return path
end

function fetch_and_build(name; verify::Bool = false)
    meta = UPSTREAM[name]
    @info "$name: downloading $(meta.url)"
    io = IOBuffer()
    Downloads.download(meta.url, io)
    parsed = parse_omlc(take!(io))

    # Convert 1/cm -> 1/m.
    as_inv_m = parsed.absorption_inv_cm .* 100

    out_path = joinpath(DATA_DIR, string(name, ".csv"))
    if verify
        tmp = tempname() * ".csv"
        write_absorption_csv(tmp, meta, parsed.λ_nm, as_inv_m)
        diff = read(`diff $out_path $tmp`, String)
        if !isempty(diff)
            @warn "$name: committed CSV differs from upstream regeneration"
            println(diff)
            return false
        else
            @info "$name: verified identical to upstream"
            return true
        end
    else
        write_absorption_csv(out_path, meta, parsed.λ_nm, as_inv_m)
        @info "$name: wrote $out_path ($(length(parsed.λ_nm)) rows)"
        return true
    end
end

function main(args)
    verify = "--verify" in args
    isdir(DATA_DIR) || mkpath(DATA_DIR)
    ok = true
    for name in sort!(collect(keys(UPSTREAM)))
        ok &= fetch_and_build(name; verify)
    end
    exit(ok ? 0 : 1)
end

abspath(PROGRAM_FILE) == @__FILE__ && main(ARGS)
