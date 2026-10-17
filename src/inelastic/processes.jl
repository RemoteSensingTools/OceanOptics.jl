# =============================================================================
# Abstract inelastic-process interface (Fell 1997 §2.2.9)
# =============================================================================
#
# An inelastic source term couples *excitation* irradiance at wavelength
# λ' into *emission* radiance at wavelength λ. In the Fell factorization
# (DESIGN §7.1, eq. 2.86–2.94),
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
# Photon conventions (the contract a solver relies on):
#   * a^I(λ') is a *photon* conversion rate: for fluorescence it already
#     includes the quantum yield, a^I = φ·a_absorber; for Raman it is the
#     Raman scattering coefficient b_R(λ'). Do not apply φ again.
#   * f^I(λ', λ) redistributes photons, per nm of emission wavelength.
#   * An energy (radiance) source therefore needs the photon-energy ratio
#     λ'/λ: `radiance_source_kernel(p, λ', λ) = a^I f^I λ'/λ`.
#   * The bundled water absorption models do not add a separate Raman loss;
#     whether b_R(λ') is added to the attenuation at λ' is the solver's
#     convention and must match the reference it is compared with.
#
# Note the scope: this file is *data and kernels only*. The two-pass
# elastic→inelastic solve itself lives in the radiative-transfer
# machinery (vSmartMOM.jl; see DESIGN §2.2.9 and the `OceanRS` extension
# planned in §7.3).
# =============================================================================

# AbstractOceanInelasticProcess is declared in `src/types/abstract_types.jl`
# so that `OceanLayer` can reference it before this file is included.
# This file provides the kernel interface only.

"""
    excitation_absorption(p::AbstractOceanInelasticProcess, λ_prime) -> FT

Absorption-side coefficient `a^I(λ')` `[1/m]` in Fell's factorization
`b^I(λ', λ) = a^I(λ') · f^I(λ', λ)`.

For fluorescence this is the quantum-yield-weighted absorption of the
fluorescing material (`φ · a_absorber(λ')`); for Raman it is the
Raman-scattering cross-section `b_w^R(λ')` taken in the "absorption
role" of the inelastic-source operator.
"""
function excitation_absorption end

"""
    emission(p::AbstractOceanInelasticProcess, λ_prime, λ) -> FT

Spectral redistribution density of emitted *photons* `f^I(λ', λ)`
`[1/nm]` (per nm of emission wavelength), normalized so that
`∫ f(λ', λ) dλ = 1` over all emission wavelengths for every `λ'`. A finite
emission grid captures only part of this integral; the solver accounts for
the remainder.
"""
function emission end

"""
    excitation_range(p::AbstractOceanInelasticProcess) -> (λ_min, λ_max)

Wavelength window `[nm]` over which the excitation integral should be
taken. Outside this window `excitation_absorption` is understood to be
negligible, allowing integrators to skip unnecessary grid points.
"""
function excitation_range end

"""
    is_isotropic(p::AbstractOceanInelasticProcess) -> Bool

Whether the process emits isotropically in angle. Isotropic emission
projects onto the `m = 0` azimuthal Fourier mode only, which the Fell
machinery exploits to halve the work. Fluorescence of solution-phase
molecules (chlorophyll, CDOM, most pigments) is isotropic;
polarization-preserving Raman is not.
"""
is_isotropic(::AbstractOceanInelasticProcess) = true

"""
    inelastic_coefficient(p, λ_prime, λ) -> FT

Combined `b^I(λ', λ) = a^I(λ') · f^I(λ', λ)` `[1/m/nm]`, the kernel of
the excitation integral that an RT solver convolves against the
pass-1 scalar irradiance `E°(τ, λ')`. Default: the product of
[`excitation_absorption`](@ref) and [`emission`](@ref). Override only
if a process supplies a joint closed form cheaper than the product.
"""
inelastic_coefficient(p::AbstractOceanInelasticProcess, λ_prime, λ) =
    excitation_absorption(p, λ_prime) * emission(p, λ_prime, λ)

"""
    radiance_source_kernel(p, λ_prime, λ) -> FT

Kernel `b^I(λ', λ) · λ'/λ` `[1/m/nm]` that converts excitation spectral
irradiance into emitted *energy*: for an isotropic process the volume
source of spectral radiance is

```math
J(\\lambda) = \\frac{1}{4\\pi} \\int b^I(\\lambda', \\lambda) \\frac{\\lambda'}{\\lambda} E_0(\\lambda')\\, d\\lambda'
```

with `E₀` the spectral scalar irradiance `[W m⁻² nm⁻¹]` (total intensity
integrated over solid angle, direct and diffuse). [`inelastic_coefficient`](@ref)
is the photon kernel; the factor `λ'/λ` is the ratio of photon energies.
On a discrete grid the solver integrates with its own bin weights and,
in wavenumber, the Jacobian `|dλ/dν| = 10⁷/ν²`.
"""
radiance_source_kernel(p::AbstractOceanInelasticProcess, λ_prime, λ) =
    inelastic_coefficient(p, λ_prime, λ) * λ_prime / λ
