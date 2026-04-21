# =============================================================================
# Spectral redistribution functions f^I(λ', λ)
# =============================================================================
#
# In Fell's (1997) factorization of an inelastic source operator (DESIGN
# §7.1),
#
#     b^I(τ; λ', λ) = a^I(τ; λ') · f^I(λ', λ)
#
# the `f^I(λ', λ)` density tells us how excitation at wavelength `λ'`
# is redistributed over emission wavelength `λ`. It is dimensioned
# `[1/nm]` and normalized `∫ f(λ', λ) dλ = 1` for every excitation `λ'`.
#
# Two shapes cover the ocean cases of interest:
#
# * `GaussianEmission(center, σ)`    — Gaussian in *wavelength* space,
#                                      centered at a fixed `center`
#                                      regardless of `λ'`. This is the
#                                      chlorophyll-a fluorescence profile
#                                      (685 nm, σ = 10.6 nm per Fell 1997
#                                      §2.2.9) and the standard CDOM
#                                      fluorescence model.
# * `WavenumberShiftEmission(Δν, σ_ν)` — Gaussian in *wavenumber* space,
#                                      centered at `ν' - Δν` where
#                                      `ν = 10⁷/λ` [1/cm]. Water Raman
#                                      emission (Haltrin & Kattawar 1991;
#                                      Δν ≈ 3400 cm⁻¹, σ_ν ≈ 85 cm⁻¹).
#                                      The loader applies the Jacobian
#                                      `|dν/dλ| = 10⁷/λ²` to convert
#                                      `g(ν) dν` into `f(λ) dλ`,
#                                      preserving `∫ f dλ = 1`.
# =============================================================================

"""
$(TYPEDEF)

Spectral redistribution `f^I(λ', λ)` `[1/nm]` for an inelastic process:
the probability density of emission wavelength `λ` given excitation at
`λ'`. Concrete subtypes implement [`redistribution`](@ref).
"""
abstract type AbstractSpectralRedistribution{FT<:AbstractFloat} end

"""
    redistribution(r, λ_prime, λ) -> FT

Evaluate the redistribution density `f(λ', λ)` `[1/nm]`. Must be
normalized so `∫ f(λ', λ) dλ = 1` for every `λ'` inside the
physically meaningful range of `r`.
"""
function redistribution end

# -----------------------------------------------------------------------------
# Fixed-wavelength Gaussian (chlorophyll / CDOM fluorescence)
# -----------------------------------------------------------------------------

"""
$(TYPEDEF)

Gaussian emission band in wavelength space centered at a fixed
wavelength, independent of excitation. Canonical example: chlorophyll-a
fluorescence,

```julia
GaussianEmission{Float64}(center_nm = 685.0, σ_nm = 10.6)
```

which matches the `685 nm` peak and `σ = 10.6 nm` (FWHM ≈ 24.9 nm) used
in Fell (1997) §2.2.9 and Gordon (1979) for the oceanic
chlorophyll-fluorescence profile.

# Fields
$(TYPEDFIELDS)
"""
struct GaussianEmission{FT<:AbstractFloat} <: AbstractSpectralRedistribution{FT}
    "Centre of the emission Gaussian `[nm]`"
    center_nm::FT
    "Standard deviation of the emission Gaussian `[nm]`"
    σ_nm::FT

    function GaussianEmission{FT}(center_nm::FT, σ_nm::FT) where {FT<:AbstractFloat}
        center_nm > zero(FT) ||
            throw(ArgumentError("GaussianEmission: center_nm must be positive, got $center_nm"))
        σ_nm > zero(FT) ||
            throw(ArgumentError("GaussianEmission: σ_nm must be positive, got $σ_nm"))
        new{FT}(center_nm, σ_nm)
    end
end

GaussianEmission{FT}(; center_nm = FT(685.0), σ_nm = FT(10.6)) where {FT<:AbstractFloat} =
    GaussianEmission{FT}(FT(center_nm), FT(σ_nm))
GaussianEmission(; kwargs...) = GaussianEmission{Float64}(; kwargs...)

function redistribution(r::GaussianEmission, λ_prime::Real, λ::Real)
    # `λ_prime` is unused — the emission shape is excitation-independent for
    # this redistribution class. Keeping it in the signature matches the
    # abstract contract so callers can dispatch uniformly.
    σ = r.σ_nm
    μ = r.center_nm
    return inv(σ * sqrt(2π)) * exp(-((λ - μ) / σ)^2 / 2)
end

# -----------------------------------------------------------------------------
# Wavenumber-shifted Gaussian (water Raman)
# -----------------------------------------------------------------------------

"""
$(TYPEDEF)

Gaussian emission band in *wavenumber* space, centered at a fixed offset
from the excitation wavenumber. Used for water Raman scattering, where
emission wavenumber `ν = ν' − Δν` with `Δν ≈ 3400 cm⁻¹` (Haltrin &
Kattawar 1991/1993).

Construction:

```julia
# Water Raman: 3400 cm⁻¹ Stokes shift, 85 cm⁻¹ Gaussian σ
WavenumberShiftEmission{Float64}(shift_cm_inv = 3400.0, σ_cm_inv = 85.0)
```

# Fields
$(TYPEDFIELDS)

# Mathematical note
The stored Gaussian is normalized in wavenumber (`∫ g(ν) dν = 1`);
[`redistribution`](@ref) applies the Jacobian `|dν/dλ| = 10⁷/λ²` so the
emitted density `f(λ', λ)` is also normalized *in wavelength*, satisfying
`∫ f dλ = 1`.
"""
struct WavenumberShiftEmission{FT<:AbstractFloat} <: AbstractSpectralRedistribution{FT}
    "Stokes shift `Δν = ν_excitation − ν_emission` `[1/cm]`"
    shift_cm_inv::FT
    "Gaussian standard deviation in wavenumber `[1/cm]`"
    σ_cm_inv::FT

    function WavenumberShiftEmission{FT}(shift_cm_inv::FT, σ_cm_inv::FT) where {FT<:AbstractFloat}
        σ_cm_inv > zero(FT) ||
            throw(ArgumentError("WavenumberShiftEmission: σ_cm_inv must be positive, got $σ_cm_inv"))
        new{FT}(shift_cm_inv, σ_cm_inv)
    end
end

WavenumberShiftEmission{FT}(; shift_cm_inv = FT(3400.0),
                              σ_cm_inv     = FT(85.0)) where {FT<:AbstractFloat} =
    WavenumberShiftEmission{FT}(FT(shift_cm_inv), FT(σ_cm_inv))
WavenumberShiftEmission(; kwargs...) = WavenumberShiftEmission{Float64}(; kwargs...)

function redistribution(r::WavenumberShiftEmission, λ_prime::Real, λ::Real)
    # Convert wavelengths to wavenumbers, evaluate Gaussian, then apply the
    # |dν/dλ| = 10⁷/λ² Jacobian so ∫ f(λ',λ) dλ = 1 holds in wavelength space.
    ν_prime  = 1e7 / λ_prime
    ν        = 1e7 / λ
    ν_center = ν_prime - r.shift_cm_inv
    σ        = r.σ_cm_inv
    g_ν      = inv(σ * sqrt(2π)) * exp(-((ν - ν_center) / σ)^2 / 2)
    return g_ν * (1e7 / λ^2)
end

# -----------------------------------------------------------------------------
# Convenience: FWHM ↔ σ conversion (useful for user-facing constructors)
# -----------------------------------------------------------------------------

"Convert Gaussian FWHM `[nm or 1/cm]` to standard deviation in the same unit."
fwhm_to_sigma(fwhm::Real) = fwhm / (2 * sqrt(2 * log(2)))
