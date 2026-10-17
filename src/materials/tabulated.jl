# =============================================================================
# Tabulated constituents — exact published IOPs, no empirical model
# =============================================================================
#
# Benchmarks and reference datasets prescribe IOPs directly: absorption and
# scattering coefficients at a handful of wavelengths, plus a phase function
# or phase matrix. These two constituents carry such values into
# `layer_optics` through the ordinary absorber / scatterer roles, so that a
# published case runs without any of the package's bio-optical defaults.
#
# Values are defined only at their tabulated wavelengths; any other
# wavelength is an error, never an interpolation.
#
# Absorption scope
# ----------------
# Published absorption is often a *bulk* value that already contains every
# absorber — e.g. a_blk = a_p + a_y + a_w in Chowdhary et al. (2020)
# Eq. (11b), or a_DP in Ottaviani & Chowdhary (2025) Table A1. Adding
# `PureWater()` or `CDOM()` on top would count water and CDOM twice. A
# `TabulatedAbsorber(...; total = true)` marks such a value, and
# `layer_optics` then rejects any other constituent that absorbs.
# =============================================================================

"""
$(TYPEDEF)

Absorber defined by tabulated absorption coefficients at exact
wavelengths.

# Fields
$(TYPEDFIELDS)

# Construction

    TabulatedAbsorber(λ, a; total = false, source = "")

`λ` `[nm]` must be strictly increasing; `a` `[1/m]` finite and
non-negative. Set `total = true` when `a` is the bulk absorption of the
whole layer (all absorbers, including water): [`layer_optics`](@ref) then
requires every other constituent in the layer to have zero absorption.

# Example

```julia
# Chowdhary et al. (2020) Table 2, AOS-II/III: a_blk = a_w
a_blk = TabulatedAbsorber([350.0, 450.0, 550.0, 650.0],
                          [0.0204, 0.0092, 0.0565, 0.3400];
                          total = true, source = "Chowdhary et al. 2020, Table 2")
absorption(a_blk, 450.0)   # 0.0092
absorption(a_blk, 500.0)   # ArgumentError: not a tabulated wavelength
```
"""
struct TabulatedAbsorber{FT} <: AbstractAbsorber{FT}
    "Wavelengths `[nm]`, strictly increasing"
    λ::Vector{FT}
    "Absorption coefficient `[1/m]` at each `λ`"
    a::Vector{FT}
    "Whether `a` is the bulk absorption of the whole layer"
    total::Bool
    "Provenance: citation, table, and any transformation applied"
    source::String

    function TabulatedAbsorber{FT}(λ, a, total, source) where {FT<:AbstractFloat}
        c = new{FT}(λ, a, total, source)
        _check_tabulated_coefficients(c.λ, c.a, "TabulatedAbsorber", "a")
        return c
    end
end

function TabulatedAbsorber(λ::AbstractVector, a::AbstractVector;
                           total::Bool = false, source::AbstractString = "")
    FT = float(promote_type(eltype(λ), eltype(a)))
    return TabulatedAbsorber{FT}(collect(FT, λ), collect(FT, a), total, String(source))
end

"""
$(TYPEDEF)

Non-absorbing scatterer defined by tabulated scattering coefficients at
exact wavelengths and a phase function (or phase matrix).

# Fields
$(TYPEDFIELDS)

# Construction

    TabulatedScatterer(λ, b, phase; source = "")

`λ` `[nm]` must be strictly increasing; `b` `[1/m]` finite and
non-negative. `phase` is any [`AbstractOceanPhaseFunction`](@ref); a
wavelength-dependent [`TabulatedGreekPhase`](@ref) must tabulate every
wavelength in `λ`. The backscattering coefficient is derived from the
phase function, `b_b = b · B`, so the two can never disagree; for a
`TabulatedGreekPhase`, `B` is that of the tabulated (possibly truncated)
expansion. Model an
absorbing component as a `TabulatedAbsorber` plus a `TabulatedScatterer`.

# Example

```julia
# Chowdhary et al. (2020) AOS-II/III: pure-seawater scattering, δ_w = 0
water_scat = TabulatedScatterer([350.0, 450.0, 550.0, 650.0],
                                [0.0134, 0.0045, 0.0019, 0.0010],
                                RayleighWaterPhase(depolarization = 0.0);
                                source = "Chowdhary et al. 2020, Table 2")
```
"""
struct TabulatedScatterer{FT, P<:AbstractOceanPhaseFunction{FT}} <: AbstractScatterer{FT}
    "Wavelengths `[nm]`, strictly increasing"
    λ::Vector{FT}
    "Scattering coefficient `[1/m]` at each `λ`"
    b::Vector{FT}
    "Phase function or phase matrix"
    phase::P
    "Provenance: citation, table, and any transformation applied"
    source::String

    function TabulatedScatterer{FT, P}(λ, b, phase, source) where {FT<:AbstractFloat,
                                                                  P<:AbstractOceanPhaseFunction{FT}}
        c = new{FT, P}(λ, b, phase, source)
        _check_tabulated_coefficients(c.λ, c.b, "TabulatedScatterer", "b")
        if is_wavelength_dependent(phase)
            for λi in c.λ   # every tabulated wavelength must have a phase column
                phase_function_moments(phase, 0, λi)
            end
        end
        return c
    end
end

function TabulatedScatterer(λ::AbstractVector, b::AbstractVector,
                            phase::AbstractOceanPhaseFunction{FT};
                            source::AbstractString = "") where {FT}
    return TabulatedScatterer{FT, typeof(phase)}(collect(FT, λ), collect(FT, b), phase,
                                                 String(source))
end

function _check_tabulated_coefficients(λ::AbstractVector, x::AbstractVector,
                                       who::AbstractString, name::AbstractString)
    _check_exact_λ_grid(λ, who)
    length(x) == length(λ) ||
        throw(DimensionMismatch("$who: $name must have one value per wavelength " *
                                "($(length(λ))), got $(length(x))"))
    all(v -> isfinite(v) && v ≥ zero(v), x) ||
        throw(ArgumentError("$who: $name must be finite and non-negative"))
    return nothing
end

absorption(c::TabulatedAbsorber, λ::Real) = c.a[_exact_λ_index(c.λ, λ, "TabulatedAbsorber")]
absorption(c::TabulatedAbsorber, λ::AbstractVector) = [absorption(c, λi) for λi in λ]

scattering(c::TabulatedScatterer, λ::Real) =
    c.b[_exact_λ_index(c.λ, λ, "TabulatedScatterer")]
scattering(c::TabulatedScatterer, λ::AbstractVector) = [scattering(c, λi) for λi in λ]

backscattering(c::TabulatedScatterer, λ::Real) =
    scattering(c, λ) * backscatter_fraction(c.phase, λ)
backscattering(c::TabulatedScatterer, λ::AbstractVector) =
    [backscattering(c, λi) for λi in λ]

phase_function(c::TabulatedScatterer) = c.phase

"""
    is_total_absorption(c::AbstractOceanConstituent) -> Bool

Whether `c`'s absorption is the bulk absorption of its whole layer, so
that no other constituent in the layer may absorb. Only a
[`TabulatedAbsorber`](@ref) built with `total = true` returns `true`.
"""
is_total_absorption(::AbstractOceanConstituent) = false
is_total_absorption(c::TabulatedAbsorber)       = c.total

function Base.show(io::IO, c::TabulatedAbsorber{FT}) where {FT}
    print(io, "TabulatedAbsorber{$FT}(", length(c.λ), " λ",
          c.total ? ", total" : "",
          isempty(c.source) ? "" : ", \"$(c.source)\"", ")")
end

function Base.show(io::IO, c::TabulatedScatterer{FT}) where {FT}
    print(io, "TabulatedScatterer{$FT}(", length(c.λ), " λ, ", c.phase,
          isempty(c.source) ? "" : ", \"$(c.source)\"", ")")
end
