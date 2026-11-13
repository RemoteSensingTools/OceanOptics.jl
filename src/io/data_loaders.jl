# =============================================================================
# Spectral data loading and interpolation
# =============================================================================
#
# All tabulated optical data in OceanOptics.jl is loaded from the `data/`
# directory as provenance-commented CSV files. Every file
# opens with a comment block giving the source citation, DOI, download URL
# and date, the spectral range and resolution, and the units of the stored
# values. `load_spectral_csv` strips those comments and returns a typed
# `SpectralTable`; `linterp` interpolates it linearly and never beyond its
# range.
# =============================================================================

"""
$(TYPEDEF)

Tabulated spectral quantity `y(λ)` on a strictly increasing wavelength
grid, for piecewise-linear interpolation with [`linterp`](@ref). Columns
are always named `λ` and `y` regardless of the physical quantity stored;
the unit interpretation is fixed by the loader's provenance header and by
the consumer that dispatches on it.

`SpectralTable(λ, y)` checks the grid and stores dense vectors of the
promoted floating-point type. The two vectors share the array type `V`,
and `Adapt.adapt(to, tbl)` moves them to a device; the range limits `lo`
and `hi` are stored as numbers, so a range check never reads a device
array.

# Fields
$(TYPEDFIELDS)
"""
struct SpectralTable{T<:AbstractFloat, V<:AbstractVector{T}}
    "Wavelengths `[nm]`, strictly increasing."
    λ::V
    "Tabulated value at each λ (units fixed by consumer)."
    y::V
    "First tabulated wavelength `[nm]`."
    lo::T
    "Last tabulated wavelength `[nm]`."
    hi::T

    # Stores its arguments as they are; `SpectralTable(λ, y)` validates.
    SpectralTable{T, V}(λ, y, lo, hi) where {T<:AbstractFloat, V<:AbstractVector{T}} =
        new{T, V}(λ, y, lo, hi)
end

function SpectralTable(λ::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    T = float(promote_type(eltype(λ), eltype(y)))
    length(λ) == length(y) ||
        throw(DimensionMismatch("SpectralTable: $(length(λ)) wavelengths for $(length(y)) values"))
    length(λ) ≥ 2 ||
        throw(ArgumentError("SpectralTable: need at least two points for interpolation"))
    x, v = _dense(T, λ), _dense(T, y)
    Base.require_one_based_indexing(x, v)
    all(isfinite, x) || throw(ArgumentError("SpectralTable: wavelengths must be finite"))
    all(isfinite, v) || throw(ArgumentError("SpectralTable: values must be finite"))
    _is_strictly_increasing(x) ||
        throw(ArgumentError("SpectralTable: λ must be strictly increasing; sort the table " *
                            "and merge repeated wavelengths"))
    lo, hi = extrema(x)
    return SpectralTable{T, typeof(x)}(x, v, lo, hi)
end

function Adapt.adapt_structure(to, t::SpectralTable{T}) where {T}
    λ = adapt(to, t.λ)
    _check_adapt_eltype(to, T, λ)
    return SpectralTable{T, typeof(λ)}(λ, adapt(to, t.y), t.lo, t.hi)
end

Base.length(t::SpectralTable) = length(t.λ)

function Base.show(io::IO, t::SpectralTable{T}) where {T}
    print(io, "SpectralTable{$T}(", length(t), " pts, λ ∈ [", t.lo, ", ", t.hi, "])")
end

# Linear interpolant of the table at `x`, without a range check: outside
# [lo, hi] the first or last segment is continued. A tabulated wavelength
# other than the first belongs to the segment that ends there; both
# segments give the tabulated value. Reads four elements, so it runs per
# element inside a broadcast on the device that holds the table.
function _interpolate(t::SpectralTable, x)
    λ, y = t.λ, t.y
    i = clamp(searchsortedfirst(λ, x) - 1, firstindex(λ), lastindex(λ) - 1)
    f = (x - λ[i]) / (λ[i+1] - λ[i])
    return (1 - f) * y[i] + f * y[i+1]
end

_throw_outside_table(t::SpectralTable, x) =
    throw(DomainError(x, "outside the tabulated range $(t.lo)–$(t.hi) nm"))

"""
    linterp(tbl::SpectralTable, x)

Return the piecewise-linear interpolant of `tbl.y` at the wavelength or vector of
wavelengths `x` `[nm]`.

Throws a `DomainError` outside `[tbl.lo, tbl.hi]`: a table never supplies values
beyond its data. A `ForwardDiff.Dual` wavelength gives the slope of the
segment as derivative. At a tabulated wavelength the interpolant has a
kink, and the derivative is the slope of one of the two segments that
meet there, whichever the number type's ordering against the tabulated
value selects.

The vector method checks the range of `x` once and interpolates in one
broadcast; `tbl` and `x` must live on the same device.
"""
function linterp(tbl::SpectralTable, x::Real)
    _within(x, tbl.lo, tbl.hi) || _throw_outside_table(tbl, x)
    return _interpolate(tbl, x)
end

function linterp(tbl::SpectralTable, x::AbstractVector{<:Real})
    isempty(x) && return _interpolate.(Ref(tbl), x)
    xmin, xmax = extrema(x)
    _within(xmin, tbl.lo, tbl.hi) || _throw_outside_table(tbl, xmin)
    _within(xmax, tbl.lo, tbl.hi) || _throw_outside_table(tbl, xmax)
    return _interpolate.(Ref(tbl), x)
end

# =============================================================================
# Data directory and provenance-CSV loader
# =============================================================================

"""
    data_path(filename) -> String

Return the absolute path of `filename` in the `data/` directory bundled with the
package.
"""
data_path(filename::AbstractString) =
    normpath(joinpath(@__DIR__, "..", "..", "data", filename))

"""
    load_spectral_csv(path; col_λ=1, col_y=2) -> SpectralTable{Float64}

Parse a provenance-commented CSV and return a [`SpectralTable`](@ref) sorted by
wavelength.

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

    # `SpectralTable` needs an increasing grid.
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
# wavelengths. These helpers validate such grids and match a query λ to the
# tabulated wavelength λ_tab with |λ − λ_tab| ≤ `_λ_MATCH_RTOL`·λ_tab. The
# tolerance covers the round-off of single-precision wavelengths (a few
# 1e-7), of grids computed as 1e7/ν from wavenumbers, and of range-generated
# grids such as `400:0.1:500`; it is far below the spacing of any published
# table. Any other wavelength is an error unless the table's owner
# interpolates explicitly.

"Relative tolerance of the match between a query wavelength and a tabulated one."
const _λ_MATCH_RTOL = 1e-6

function _check_exact_λ_grid(λ::AbstractVector{<:Real}, who::AbstractString)
    isempty(λ) && throw(ArgumentError("$who: empty wavelength grid"))
    all(x -> isfinite(x) && x > zero(x), λ) ||
        throw(ArgumentError("$who: wavelengths must be finite and positive"))
    _is_strictly_increasing(λ) ||
        throw(ArgumentError("$who: wavelengths must be strictly increasing"))
    # Farther apart than twice the match tolerance: a query matches at most one.
    all(i -> λ[i+1] - λ[i] > 2 * _λ_MATCH_RTOL * λ[i+1], firstindex(λ):lastindex(λ)-1) ||
        throw(ArgumentError("$who: neighbouring wavelengths must differ by more than " *
                            "$(2 * _λ_MATCH_RTOL) of their value (the wavelength-match " *
                            "tolerance); check the table for duplicated rows"))
    return nothing
end

# Index of the tabulated wavelength that matches `λ`, or `nothing`.
function _matching_λ_index(λ_tab::AbstractVector, λ::Real)
    i = searchsortedfirst(λ_tab, λ)
    for k in (i - 1, i)
        checkbounds(Bool, λ_tab, k) && abs(λ_tab[k] - λ) ≤ _λ_MATCH_RTOL * λ_tab[k] && return k
    end
    return nothing
end

function _exact_λ_index(λ_tab::AbstractVector, λ::Real, who::AbstractString)
    k = _matching_λ_index(λ_tab, λ)
    k === nothing || return k
    shown = string(join(first(λ_tab, 8), ", "), length(λ_tab) > 8 ? ", …" : "")
    throw(ArgumentError("$who: λ = $λ nm is not a tabulated wavelength (tabulated: $shown; " *
                        "match tolerance $(_λ_MATCH_RTOL) relative); this table is not " *
                        "interpolated"))
end

"""
$(TYPEDEF)

How a table defined at discrete wavelengths answers a query between them.

[`ExactLookup`](@ref) admits tabulated wavelengths only, [`LinearLookup`](@ref)
interpolates linearly between neighbours and [`NearestLookup`](@ref) takes the nearest
entry. Constructors accept the names `:exact`, `:linear` and `:nearest` for them.
"""
abstract type AbstractTableLookup end

"""
$(TYPEDEF)

Table lookup that admits tabulated wavelengths only (to a relative `1e-6`); any other
wavelength throws.
"""
struct ExactLookup <: AbstractTableLookup end

"""
$(TYPEDEF)

Table lookup that interpolates linearly between neighbouring tabulated wavelengths,
never beyond the table.
"""
struct LinearLookup <: AbstractTableLookup end

"""
$(TYPEDEF)

Table lookup that uses the entry of the nearest tabulated wavelength.
"""
struct NearestLookup <: AbstractTableLookup end

_lookup_name(::ExactLookup)   = :exact
_lookup_name(::LinearLookup)  = :linear
_lookup_name(::NearestLookup) = :nearest

# The lookup named by `lookup` (a name or a lookup object), which must be one of
# the instances in `allowed`.
function _table_lookup(lookup, allowed::Tuple, who::AbstractString)
    for l in allowed
        (lookup === l || lookup === _lookup_name(l)) && return l
    end
    throw(ArgumentError("$who: lookup must be " *
                        join((":$(_lookup_name(l))" for l in allowed), " or ") *
                        ", got $(repr(lookup))"))
end

# Value of the table `(λ_tab, y)` at `λ`. `ExactLookup`: the entry of the
# matching wavelength, in the table's precision. `LinearLookup`: that entry,
# or else the linear interpolant between the two neighbouring entries (never
# beyond the table), in the type that `λ` and the table promote to.
_tabulated_value(::ExactLookup, λ_tab::AbstractVector{FT}, y::AbstractVector{FT}, λ::Real,
                 who::AbstractString) where {FT} = y[_exact_λ_index(λ_tab, λ, who)]

function _tabulated_value(::LinearLookup, λ_tab::AbstractVector{FT}, y::AbstractVector{FT},
                          λ::Real, who::AbstractString) where {FT}
    k = _matching_λ_index(λ_tab, λ)
    k === nothing || return convert(_result_type(FT, λ), y[k])
    first(λ_tab) < λ < last(λ_tab) ||
        throw(DomainError(λ, "$who: outside the tabulated range $(first(λ_tab))–" *
                             "$(last(λ_tab)) nm; tables are interpolated, never extrapolated"))
    i = searchsortedlast(λ_tab, λ)
    t = (λ - λ_tab[i]) / (λ_tab[i+1] - λ_tab[i])
    return convert(_result_type(FT, λ), (1 - t) * y[i] + t * y[i+1])
end
