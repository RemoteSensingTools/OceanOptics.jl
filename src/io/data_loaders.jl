# =============================================================================
# Spectral data loading and interpolation
# =============================================================================
#
# All tabulated optical data in OceanOptics.jl is loaded from the `data/`
# directory as provenance-commented CSV files (see DESIGN §3.4). Every file
# opens with a comment block giving the source citation, DOI, download URL
# and date, the spectral range and resolution, and the units of the stored
# values. `load_spectral_csv` strips those comments and returns a typed
# `SpectralTable`; `linterp` performs piecewise-linear interpolation via
# `Interpolations.jl`'s `linear_interpolation` with `Flat()` boundary
# extrapolation (AD-transparent: ForwardDiff `Dual` inputs propagate
# through the interior and receive zero partials outside the grid).
# =============================================================================

"""
$(TYPEDEF)

Tabulated spectral quantity `y(λ)` with monotone `λ` axis, backed by an
`Interpolations.linear_interpolation` object for efficient repeated
queries. Columns are always named `λ` and `y` regardless of the physical
quantity stored; the unit interpretation is fixed by the loader's
provenance header and by the consumer that dispatches on it.

# Fields
$(TYPEDFIELDS)
"""
struct SpectralTable{T<:AbstractFloat, I}
    "Wavelengths `[nm]`, sorted ascending."
    λ::Vector{T}
    "Tabulated value at each λ (units fixed by consumer)."
    y::Vector{T}
    "Cached `Interpolations.jl` linear interpolant (queries are range-checked by `linterp`)."
    itp::I
end

function SpectralTable(λ::Vector{T}, y::Vector{T}) where {T<:AbstractFloat}
    length(λ) == length(y) ||
        throw(ArgumentError("SpectralTable: λ and y must have equal length"))
    length(λ) ≥ 2 ||
        throw(ArgumentError("SpectralTable: need at least two points for interpolation"))
    issorted(λ) ||
        throw(ArgumentError("SpectralTable: λ must be sorted ascending"))
    itp = linear_interpolation(λ, y; extrapolation_bc = Flat())
    return SpectralTable{T, typeof(itp)}(λ, y, itp)
end

Base.length(t::SpectralTable) = length(t.λ)

function Base.show(io::IO, t::SpectralTable{T}) where {T}
    print(io, "SpectralTable{$T}(", length(t), " pts, λ ∈ [",
          first(t.λ), ", ", last(t.λ), "])")
end

"""
    linterp(tbl::SpectralTable, x) -> value

Piecewise-linear interpolation of `tbl.y` at wavelength `x`, delegating
to `Interpolations.jl`. Throws a `DomainError` outside
`[tbl.λ[1], tbl.λ[end]]`: a table never supplies values beyond its data.
ForwardDiff `Dual` inputs propagate through the interpolation.
"""
function linterp(tbl::SpectralTable, x)
    first(tbl.λ) ≤ x ≤ last(tbl.λ) ||
        throw(DomainError(x, "outside the tabulated range $(first(tbl.λ))–$(last(tbl.λ)) nm"))
    return tbl.itp(x)
end

# =============================================================================
# Data directory and provenance-CSV loader
# =============================================================================

"""
    data_path(filename) -> String

Absolute path to `filename` under the package's bundled `data/` directory.
Used both by the materials code (to load shipped tables at module init)
and by the `scripts/fetch_data.jl` reproducibility script.
"""
data_path(filename::AbstractString) =
    normpath(joinpath(@__DIR__, "..", "..", "data", filename))

"""
    load_spectral_csv(path; col_λ=1, col_y=2) -> SpectralTable{Float64}

Parse a provenance-commented CSV and return a sorted `SpectralTable`.

Recognized structure:
  * Comment lines begin with `#` and are discarded.
  * The first non-comment, non-numeric line (the column header) is also
    discarded via `tryparse`.
  * Numeric rows are parsed as `Float64`; `col_λ` and `col_y` (1-based)
    select the wavelength column and the value column. Extra columns are
    ignored so multi-column files (e.g. Bricaud A(λ)/E(λ)) can share one
    loader.

Throws `ArgumentError` when `path` does not exist or when no rows parse.
"""
function load_spectral_csv(path::AbstractString; col_λ::Int = 1, col_y::Int = 2)
    isfile(path) ||
        throw(ArgumentError("Data file not found: $path"))
    col_λ ≥ 1 && col_y ≥ 1 ||
        throw(ArgumentError("Column indices must be positive (got col_λ=$col_λ, col_y=$col_y)"))

    λs, ys = Float64[], Float64[]
    for raw in eachline(path)
        line = strip(raw)
        (isempty(line) || startswith(line, '#')) && continue
        parts = split(line, ',')
        length(parts) < max(col_λ, col_y) && continue
        λ = tryparse(Float64, strip(parts[col_λ]))
        y = tryparse(Float64, strip(parts[col_y]))
        (λ === nothing || y === nothing) && continue
        push!(λs, λ)
        push!(ys, y)
    end
    isempty(λs) &&
        throw(ArgumentError("No numeric rows parsed from $path"))

    # Guarantee sorted output for downstream searchsortedlast.
    if !issorted(λs)
        perm = sortperm(λs)
        λs  .= λs[perm]
        ys  .= ys[perm]
    end
    return SpectralTable(λs, ys)
end

"""
    _read_csv_columns(path) -> Dict{String, Vector{Float64}}

Read a provenance-commented, multi-column numeric CSV (as in `data/`): `#`
comment lines are skipped, the first remaining line names the columns, and
every following row must parse as numbers. Used for reference datasets
such as `data/reference/loisel2023_hydrolight_iops.csv`.
"""
function _read_csv_columns(path::AbstractString)
    lines  = filter(l -> !isempty(l) && !startswith(l, '#'), strip.(readlines(path)))
    isempty(lines) && throw(ArgumentError("No header or data rows in $path"))
    names  = String.(strip.(split(first(lines), ',')))
    values = [Float64[] for _ in names]
    for (k, line) in enumerate(lines[2:end])
        fields = split(line, ',')
        length(fields) == length(names) ||
            throw(ArgumentError("$path: data row $k has $(length(fields)) fields, " *
                                "expected $(length(names))"))
        for (v, f) in zip(values, fields)
            push!(v, parse(Float64, strip(f)))
        end
    end
    return Dict(zip(names, values))
end

# -----------------------------------------------------------------------------
# Exact-wavelength tables (no interpolation)
# -----------------------------------------------------------------------------
#
# Tabulated benchmark inputs are defined only at their published
# wavelengths. These helpers validate such grids and look wavelengths up
# exactly — to a few ulps of the less precise of table and query, so that
# range-generated grids such as `400:0.1:500` and Float32 queries of a
# Float64 table still match. Any other wavelength is an error rather than
# an interpolation.

function _check_exact_λ_grid(λ::AbstractVector{<:Real}, who::AbstractString)
    isempty(λ) && throw(ArgumentError("$who: empty wavelength grid"))
    all(x -> isfinite(x) && x > zero(x), λ) ||
        throw(ArgumentError("$who: wavelengths must be finite and positive"))
    _is_strictly_increasing(λ) ||
        throw(ArgumentError("$who: wavelengths must be strictly increasing"))
    return nothing
end

function _exact_λ_index(λ_tab::AbstractVector{FT}, λ::Real,
                        who::AbstractString) where {FT}
    i   = searchsortedfirst(λ_tab, λ)
    tol = 8 * _λ_match_eps(FT, λ)
    for k in (i - 1, i)
        checkbounds(Bool, λ_tab, k) && abs(λ_tab[k] - λ) ≤ tol * λ_tab[k] && return k
    end
    shown = join(first(λ_tab, 8), ", ") * (length(λ_tab) > 8 ? ", …" : "")
    throw(ArgumentError("$who: λ = $λ nm is not a tabulated wavelength " *
                        "(tabulated: $shown); exact-wavelength tables are not interpolated"))
end

_λ_match_eps(::Type{FT}, λ::AbstractFloat) where {FT} = max(eps(FT), eps(typeof(λ)))
_λ_match_eps(::Type{FT}, ::Real) where {FT}          = eps(FT)
