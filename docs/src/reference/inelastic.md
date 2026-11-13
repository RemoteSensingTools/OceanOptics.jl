# Inelastic processes

Fluorescence, water Raman, and CDOM fluorescence packaged as kernels
consumed by a two-pass RT solver (Fell 1997 §2.2.9), and the light
absorbed per layer and component that drives them.

## Conventions

**Conversion units.** Every process factorizes into a conversion coefficient
`κ(λ′) = `[`excitation_absorption`](@ref) `[1/m]` and a redistribution density
`f(λ′, λ) = `[`emission`](@ref) `[1/nm]` with `∫ f dλ = 1`. [`conversion_units`](@ref)
states whether `κ` counts photons or energy, and the two kernels follow by dispatch:

| `conversion_units(p)` | energy kernel `K_E` ([`radiance_source_kernel`](@ref)) | photon kernel `K_N` ([`inelastic_coefficient`](@ref)) | processes |
|:--|:--|:--|:--|
| [`PhotonConversion`](@ref) | `κ f λ′/λ` | `κ f` | [`IsotropicFluorescence`](@ref), [`CDOMFluorescence`](@ref); `WaterRaman(conversion = PhotonConversion())` |
| [`EnergyConversion`](@ref) | `κ f` | `κ f λ/λ′` | [`WaterRaman`](@ref) (default) |

In both cases `K_N = K_E λ/λ′`, and the power emitted per unit volume and per nm of
emission wavelength is

```math
4\pi \bar S(\lambda) = \int K_E(\lambda', \lambda)\, E_0(\lambda')\, d\lambda'
```

with `E₀` the spectral scalar irradiance `[W m⁻² nm⁻¹]`. [`inelastic_emission`](@ref)
integrates the same kernel over a layer.

**Fluorescence** is a photon conversion: `κ = φ Σ_k η_k a_k(λ′)` with the quantum yield
`φ` in photons emitted per photon delivered to the fluorescing pigment, so
`∫ K_N dλ = κ` and the emitted energy carries the photon-energy ratio `λ′/λ`.

**Water Raman** is an energy conversion by default: the Raman scattering coefficient
`b_R(λ′) = b_R(λ_ref) (λ_ref/λ′)^n` (defaults `2.7·10⁻⁴ m⁻¹` at 488 nm, `n = 5.5`) is
defined so that `b_R(λ′) E₀(λ′)` is the power emitted over the whole Raman band per unit
volume and per nm of excitation wavelength,

```math
\int K_E(\lambda', \lambda)\, d\lambda = b_R(\lambda') ,
```

with no factor `λ′/λ`. This is the definition of Mobley et al. (1993, Appl. Opt. 32,
7484; Problem 7 and Eqs. A3–A4): in Problem 7 all light Raman scattered at 417 nm
appears at 486 nm and the source of radiance is `b_R E₀/(4π)` times the phase function.
`WaterRaman(conversion = PhotonConversion())` treats the same number as a photon
coefficient instead (`∫ K_N dλ = b_R`); its energy source is lower by the band mean of
`λ′/λ`: 0.881 at 350 nm excitation, 0.858 at 417 nm, 0.834 at 488 nm, 0.796 at 600 nm.

**Raman angular distribution and depolarization.** `WaterRaman` stores the
depolarization ratio `ρ` for linearly polarized excitation (`linear_depolarization`,
default 0.17 from Mobley et al. 1993, Eq. 9). The phase function is

```math
p_R(\cos\Theta) = \frac{1}{4\pi}\left[1 + \beta_2 P_2(\cos\Theta)\right],
\qquad \beta_2 = \frac{1-\rho}{2(1+2\rho)} = 0.3097 \ \text{for}\ \rho = 0.17 ,
```

identical to Mobley's Eq. 9, `(3/16π) (1+3ρ)/(1+2ρ) [1 + (1−ρ)/(1+3ρ) cos²Θ]`.
The ratio for unpolarized (natural) excitation is `ρ_n = 2ρ/(1+ρ)` (0.2906;
[`natural_depolarization`](@ref)), and `β₂ = (1−ρ_n)/(2+ρ_n)` is the same number.
[`raman_phase`](@ref) returns the [`RayleighWaterPhase`](@ref) with depolarization
`ρ_n`, whose polarized Greek expansion is the Raman phase matrix.

**Raman loss.** The bundled water absorption models do not add a separate Raman loss;
whether `b_R(λ′)` is added to the attenuation at the excitation wavelength is the
solver's convention and must match the reference it is compared with (Mobley et al.
1993 take it to be contained in the Smith & Baker absorption).

**Pigment partitions outside their tables.** The pigment spectra of
[`ProspectPigments`](@ref) (400–800 nm) and [`Bidigare1990Pigments`](@ref) (400–750 nm)
do not cover the UV. By default (`outside = `[`HoldEdgeShares`](@ref)`()`) the shares of
the nearest tabulated wavelength apply outside the table, so absorption below 400 nm is
attributed to the pigment groups as at 400 nm and [`chlorophyll_fluorescence`](@ref)
with `excitation_range = (350, 700)` is excited continuously across 400 nm; with all
transfer efficiencies equal to 1 its excitation equals that of the unpartitioned
phytoplankton. The UV shares are an extrapolation, not data. A partition built with
`outside = `[`UnassignedOutsideTable`](@ref)`()` puts all absorption outside its table
in `:unassigned` (η = 0); the excitation then steps to zero just below 400 nm, which
removes 21.9 % of the 350–700 nm excitation integral at Chl = 0.5 mg m⁻³. Under both
rules `Bidigare1990Pigments` assigns nothing above about 732 nm, where all its group
spectra vanish. The default `excitation_range` of `chlorophyll_fluorescence` is
(400, 700) nm.

**Absorption that excites chlorophyll fluorescence.** With [`Bricaud1998`](@ref) the
absorption of `Phytoplankton` is the particulate absorption `a_p` of Case-1 water,
detritus included; the fluorescence excitation and the pigment shares are fractions
of that total.

## Spectral redistribution

```@docs
AbstractSpectralRedistribution
TabulatedEmission
chloroplast_emission
GaussianEmission
WavenumberShiftEmission
redistribution
fwhm_to_sigma
```

## Concrete processes

```@docs
ExcitationPathway
IsotropicFluorescence
CDOMFluorescence
WaterRaman
chlorophyll_fluorescence
raman_peak_wavelength
raman_phase
natural_depolarization
```

## Trait interface

```@docs
AbstractOceanInelasticProcess
AbstractConversionUnits
PhotonConversion
EnergyConversion
conversion_units
excitation_absorption
emission
excitation_range
excitation_pathways
is_isotropic
quantum_yield
inelastic_coefficient
radiance_source_kernel
```

## Light absorbed per layer

```@docs
absorbed_power
absorbed_photons
absorbed_par
inelastic_emission
excitation_optical_depth
absorber_index
```
