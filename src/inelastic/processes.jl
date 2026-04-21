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
# Package-downstream code (vSmartMOM's `AbstractRamanType` adapter, users'
# own solvers) consumes those kernels against an excitation-field
# `E°(τ, λ')` produced by a first-pass elastic solve.
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

Spectral redistribution density `f^I(λ', λ)` `[1/nm]`, normalized so
`∫ f(λ', λ) dλ = 1` for every `λ'`.
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
