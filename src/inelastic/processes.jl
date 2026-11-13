# =============================================================================
# Abstract inelastic-process interface (Fell 1997 §2.2.9)
# =============================================================================
#
# An inelastic source term couples *excitation* irradiance at wavelength
# λ' into *emission* radiance at wavelength λ. In the Fell factorization
# (Fell 1997, Eqs. 2.86–2.94),
#
#     J^I(τ, μ, φ, λ) = (1/c(τ)) · ∫_{λ'} b^I(λ', λ) · ⟨L(λ')⟩_Ω
#                                           · p^I(cos θ; λ', λ) dλ'
#
#     b^I(λ', λ) = a^I(λ') · f^I(λ', λ)
#
# factorizes every inelastic source (chlorophyll fluorescence, water
# Raman, CDOM fluorescence, phycobilipigments, …) into three pieces:
#
#   a^I(λ')   — excitation-side absorption coefficient `[1/m]`
#   f^I(λ',λ) — spectral redistribution density `[1/nm]`, ∫f dλ = 1
#   p^I(θ)    — angular redistribution (isotropic ⇒ Fourier mode m=0 only)
#
# Concrete process types implement the three trait functions exported here.
#
# Conversion units (the contract a solver relies on):
#   * Each process declares, through `conversion_units(p)`, whether a^I(λ')
#     converts photons (`PhotonConversion`) or energy (`EnergyConversion`),
#     and f^I(λ', λ) redistributes that same quantity per nm of emission
#     wavelength, ∫ f^I dλ = 1.
#   * Fluorescence is a photon conversion: a^I = φ·a_absorber already
#     includes the quantum yield (photons per photon). Do not apply φ again.
#   * Water Raman is an energy conversion by default: b_R(λ') E₀(λ') is the
#     power emitted over the Raman band per unit volume and per nm of
#     excitation wavelength (Mobley et al. 1993, Problem 7 and Appendix A).
#   * `radiance_source_kernel` is the energy kernel and
#     `inelastic_coefficient` the photon kernel; they differ by the
#     photon-energy ratio λ'/λ, applied on the side that the convention of
#     the process requires.
#   * The bundled water absorption models do not add a separate Raman loss;
#     whether b_R(λ') is added to the attenuation at λ' is the solver's
#     convention and must match the reference it is compared with.
#
# Note the scope: this file is *data and kernels only*. The two-pass
# elastic→inelastic solve itself lives in the radiative-transfer
# machinery (vSmartMOM.jl), which calls the kernels of this file.
# =============================================================================

# AbstractOceanInelasticProcess is declared in `src/types/abstract_types.jl`
# so that `OceanLayer` can reference it before this file is included.
# This file provides the kernel interface only.

"""
    AbstractConversionUnits

Quantity that the conversion coefficient of an inelastic process counts:
[`PhotonConversion`](@ref) or [`EnergyConversion`](@ref).

See [`conversion_units`](@ref).
"""
abstract type AbstractConversionUnits end

"""
    PhotonConversion()

Conversion units of a process whose coefficient `a^I(λ′)` `[1/m]` counts photons:
`a^I(λ′) f^I(λ′, λ)` is the number of photons emitted per nm at `λ` per unit photon
fluence at `λ′`.

The emitted energy carries the photon-energy ratio `λ′/λ`. Fluorescence, defined by a
quantum yield in photons per photon, is a photon conversion.
"""
struct PhotonConversion <: AbstractConversionUnits end

"""
    EnergyConversion()

Conversion units of a process whose coefficient `a^I(λ′)` `[1/m]` counts energy:
`a^I(λ′) f^I(λ′, λ) E₀(λ′)` is the power emitted per unit volume and per nm at `λ` for
the spectral scalar irradiance `E₀(λ′)` `[W m⁻² nm⁻¹]`, with no photon-energy ratio.

This is the definition of the Raman scattering coefficient `b_R` of Mobley et al.
(1993, Appl. Opt. 32, 7484; Problem 7 and Eqs. A3–A4), in which the radiance source is
`b_R E₀ / (4π)` times the angular and spectral distributions.
"""
struct EnergyConversion <: AbstractConversionUnits end

"""
    conversion_units(p::AbstractOceanInelasticProcess) -> AbstractConversionUnits

Return whether the coefficient [`excitation_absorption`](@ref) of `p` converts photons
([`PhotonConversion`](@ref)) or energy ([`EnergyConversion`](@ref)).

With `κ = excitation_absorption(p, λ′)` and `f = emission(p, λ′, λ)`, `∫ f dλ = 1`:

| units | energy kernel [`radiance_source_kernel`](@ref) | photon kernel [`inelastic_coefficient`](@ref) |
|:------|:-----------|:-----------|
| `PhotonConversion` | `κ f λ′/λ` | `κ f` |
| `EnergyConversion` | `κ f` | `κ f λ/λ′` |

# Implementation
Every concrete process defines a method; there is no default.
"""
function conversion_units end

# Energy and photon kernels from the conversion kernel κ·f of a process.
_energy_kernel(::PhotonConversion, κf, λ_prime, λ) = κf * λ_prime / λ
_energy_kernel(::EnergyConversion, κf, λ_prime, λ) = κf
_photon_kernel(::PhotonConversion, κf, λ_prime, λ) = κf
_photon_kernel(::EnergyConversion, κf, λ_prime, λ) = κf * λ / λ_prime

"""
    excitation_absorption(p::AbstractOceanInelasticProcess, λ_prime) -> FT

Return the absorption-side coefficient `a^I(λ')` `[1/m]` of the factorization
`b^I(λ', λ) = a^I(λ') · f^I(λ', λ)` of Fell (1997).

For fluorescence this is the quantum-yield-weighted absorption of the
fluorescing material (`φ · a_absorber(λ')`), a photon conversion; for Raman it is the
Raman scattering coefficient `b_R(λ')`. [`conversion_units`](@ref) says whether the
coefficient counts photons or energy.
"""
function excitation_absorption end

"""
    emission(p::AbstractOceanInelasticProcess, λ_prime, λ) -> FT

Return the spectral redistribution density `f^I(λ', λ)` `[1/nm]` (per nm of emission
wavelength) of the quantity that `p` converts: photons for a
[`PhotonConversion`](@ref), energy for an [`EnergyConversion`](@ref).

It is normalized so that `∫ f(λ', λ) dλ = 1` over all emission wavelengths for every
`λ'`. A finite emission grid captures only part of this integral; the solver accounts
for the remainder.
"""
function emission end

"""
    excitation_range(p::AbstractOceanInelasticProcess) -> (λ_min, λ_max)

Return the wavelength window `[nm]` over which the excitation integral is taken.

Outside this window `excitation_absorption` is treated as negligible, so integrators can
skip those grid points.
"""
function excitation_range end

"""
    is_isotropic(p::AbstractOceanInelasticProcess) -> Bool

Return whether the process emits isotropically in angle.

Isotropic emission projects onto the `m = 0` azimuthal Fourier mode only. Fluorescence
of solution-phase molecules (chlorophyll, CDOM, most pigments) is isotropic;
polarization-preserving Raman is not.
"""
is_isotropic(::AbstractOceanInelasticProcess) = true

"""
    inelastic_coefficient(p, λ_prime, λ) -> FT

Return the photon kernel `K_N(λ', λ)` `[1/m/nm]`: photons emitted per nm at `λ` per unit
photon fluence at `λ'`.

With `κ = excitation_absorption(p, λ')` and `f = emission(p, λ', λ)` it is `κ f` for a
[`PhotonConversion`](@ref) and `κ f λ/λ'` for an [`EnergyConversion`](@ref)
(see [`conversion_units`](@ref)). For fluorescence its integral over `λ` is
`φ Σ_k η_k a_k(λ')`.
"""
inelastic_coefficient(p::AbstractOceanInelasticProcess, λ_prime, λ) =
    _photon_kernel(conversion_units(p),
                   excitation_absorption(p, λ_prime) * emission(p, λ_prime, λ), λ_prime, λ)

"""
    radiance_source_kernel(p, λ_prime, λ) -> FT

Return the energy kernel `K_E(λ', λ)` `[1/m/nm]` that converts excitation spectral
irradiance into emitted energy.

The power emitted per unit volume and per nm at `λ` is `∫ K_E(λ', λ) E₀(λ') dλ'`, with
`E₀` the spectral scalar irradiance `[W m⁻² nm⁻¹]` (total intensity integrated over
solid angle, direct and diffuse). For an isotropic process the volume source of
spectral radiance is

```math
J(\\lambda) = \\frac{1}{4\\pi} \\int K_E(\\lambda', \\lambda) E_0(\\lambda')\\, d\\lambda'
```

With `κ = excitation_absorption(p, λ')` and `f = emission(p, λ', λ)`, `K_E = κ f λ'/λ`
for a [`PhotonConversion`](@ref) (fluorescence; `λ'/λ` is the ratio of photon energies)
and `K_E = κ f` for an [`EnergyConversion`](@ref) (see [`conversion_units`](@ref)). For
the default [`WaterRaman`](@ref), an energy conversion, `∫ K_E dλ = b_R(λ')`.
[`inelastic_coefficient`](@ref) is the photon kernel `K_N = K_E λ/λ'`.

On a discrete grid the solver integrates with its own bin weights and,
in wavenumber, the Jacobian `|dλ/dν| = 10⁷/ν²`.
"""
radiance_source_kernel(p::AbstractOceanInelasticProcess, λ_prime, λ) =
    _energy_kernel(conversion_units(p),
                   excitation_absorption(p, λ_prime) * emission(p, λ_prime, λ), λ_prime, λ)
