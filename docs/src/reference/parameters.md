# Physical parameters

Every physical parameter of the package is a field of a model type or a keyword
argument with a default. This page lists them per model with the default, the unit,
where the number comes from, and how to set it. Coefficient sets that define a
published model as a whole (the Quan–Fry polynomial, the Zhang–Hu–He formulas, the
tabulated data files) are named constants or data files; the choice between such
models is made by type.

"Package default" marks a value that is a modelling choice of this package and not
taken from a paper that was checked for this release.

## Pure water: [`PureWater`](@ref)

| parameter | symbol | default | unit | source | set with |
|:--|:--|:--|:--|:--|:--|
| temperature | `T` | 20 | °C | package default | `PureWater(temperature = 5)` |
| salinity | `S` | 35 | PSU | package default | `PureWater(salinity = 0)` |
| depolarization ratio (unpolarized light, 90°) | `δ` | 0.039 | – | Farinato & Rowell (1976), as used by Zhang, Hu & He (2009) | `PureWater(depolarization = 0.09)` |
| absorption table | – | [`IOCCG2018`](@ref) | – | see the model types | `PureWater(absorption_model = PopeFry1997())` |
| scattering model | – | [`ZhangHu2009`](@ref) | – | see the model types | `PureWater(scattering_model = Morel1974())` |

`δ` is one stored value: it sets the phase matrix ([`RayleighWaterPhase`](@ref) from
`scattering_phase`) and the magnitude of the `ZhangHu2009` scattering,
`b ∝ (6 + 6δ)/(6 − 7δ) · (2 + δ)/(1 + δ)`. `Morel1974` scattering is an empirical
seawater law and uses neither `δ`, `T` nor `S`. Temperature enters the absorption
through the coefficient `Ψ_T` of Röttgers et al. (2014), applied from 550 nm upward
relative to the 22 °C of the Pope & Fry and Kou tables; salinity does not enter the
absorption.

Fixed coefficient sets: `Morel1974` (`0.00193 m⁻¹` at 550 nm, exponent 4.32); the
Zhang, Hu & He (2009) formulas with their auxiliary fits; the absorption tables in
`data/`.

## Rayleigh phase function: [`RayleighWaterPhase`](@ref)

| parameter | symbol | default | unit | source | set with |
|:--|:--|:--|:--|:--|:--|
| depolarization ratio (unpolarized light, 90°) | `ρ` | 0.039 | – | Farinato & Rowell (1976) | `RayleighWaterPhase(depolarization = 0.09)` |

`β₂ = (1 − ρ)/(2 + ρ)`. The default is the same constant as that of `PureWater`.

## Refractive index: [`FixedRefractiveIndex`](@ref), [`QuanFry1995`](@ref)

| parameter | symbol | default | unit | source | set with |
|:--|:--|:--|:--|:--|:--|
| index relative to air | `n` | none (required) | – | – | `FixedRefractiveIndex(1.34)` |
| temperature | `T` | none (required) | °C | – | `QuanFry1995(temperature = 20, salinity = 35)` |
| salinity | `S` | none (required) | PSU | – | as above |

`column_optics` has no default refractive index. The ten Quan & Fry (1995)
coefficients are a fixed set.

## CDOM: [`CDOM`](@ref)

`a_g(λ) = a_ref · exp(−S (λ − λ_ref))`

| parameter | symbol | default | unit | source | set with |
|:--|:--|:--|:--|:--|:--|
| reference absorption | `a_ref` | 0.03 | m⁻¹ | package default (open ocean) | `CDOM(a_ref = 0.5)` |
| reference wavelength | `λ_ref` | 440 | nm | package default | `CDOM(λ_ref = 400)` |
| spectral slope | `S` | 0.014 | nm⁻¹ | Bricaud et al. (1981), as given by Fell (1997), Eq. 5.4 | `CDOM(slope = 0.017)` |

## Particles of Case-1 water: [`Phytoplankton`](@ref), [`GordonMorel1983`](@ref)

`b_p(λ) = b_ref · Chl^m · (λ_ref/λ)^γ`; absorption `a_p = A(λ) Chl^{E(λ)}` from the
[`Bricaud1998`](@ref) tables (fixed data).

| parameter | symbol | default | unit | source | set with |
|:--|:--|:--|:--|:--|:--|
| chlorophyll-a concentration | `Chl` | none (required) | mg m⁻³ | – | `Phytoplankton(Chl = 0.5)` |
| scattering at `λ_ref` for `Chl = 1` | `b_ref` | 0.30 | m⁻¹ | Gordon & Morel (1983) as cited by Fell (1997), Eq. 5.5 (`0.30 ± 0.15`) | `Phytoplankton(Chl = 0.5, scattering_model = GordonMorel1983(b_ref = 0.45))` |
| chlorophyll exponent | `m` | 0.62 | – | same | `GordonMorel1983(chl_exponent = 0.766)` |
| reference wavelength | `λ_ref` | 550 | nm | same | `GordonMorel1983(λ_ref = 660)` |
| spectral exponent | `γ` | 1 | – | same | `GordonMorel1983(spectral_exponent = 0.5)` |
| backscatter fraction of the Fournier–Forand phase function | `B` | 0.01 | – | package default; equals `b_bph/b_ph` of the database of Loisel et al. (2023) | `Phytoplankton(Chl = 0.5, B = 0.015)` |
| carotenoid-to-chlorophyll ratio | – | none (required) | – | – | `ProspectPigments(carotenoid_ratio = 0.3)` |
| pigment concentrations | – | `chl_a` required, others 0 | mg m⁻³ | – | `Bidigare1990Pigments(chl_a = 0.5, ppc = 0.1)` |
| shares outside the pigment tables | – | [`HoldEdgeShares`](@ref) | – | package default | `ProspectPigments(carotenoid_ratio = 0.3, outside = UnassignedOutsideTable())` |

## Non-algal particles: [`NonAlgalParticles`](@ref)

[`BabinBricaud2003`](@ref): `a(λ) = a_ref · exp(−S (λ − λ_ref))`,
`b(λ) = b_ref · (λ_ref/λ)^γ`.

| parameter | symbol | default | unit | source | set with |
|:--|:--|:--|:--|:--|:--|
| reference absorption | `a_ref` | 0.05 | m⁻¹ | package default | `NonAlgalParticles(a_ref = 0.1)` |
| reference scattering | `b_ref` | 0.45 | m⁻¹ | package default | `NonAlgalParticles(b_ref = 0.6)` |
| reference wavelength | `λ_ref` | 443 | nm | package default | `NonAlgalParticles(λ_ref = 440)` |
| absorption slope | `S` | 0.011 | nm⁻¹ | package default (Babin et al. 2003 was not available for checking) | `NonAlgalParticles(slope = 0.0123)` |
| spectral exponent of scattering | `γ` | 1 | – | package default | `NonAlgalParticles(spectral_exponent = 0.5)` |
| backscatter fraction | `B` | 0.02 | – | package default | `NonAlgalParticles(B = 0.03)` |

[`FellMineral`](@ref): `c(λ) = c_ref · (λ_ref/λ)^γ`, `a = (1 − ω) c`, `b = ω c`.

| parameter | symbol | default | unit | source | set with |
|:--|:--|:--|:--|:--|:--|
| beam attenuation at `λ_ref` | `c_ref` | 0.40 | m⁻¹ | package default | `NonAlgalParticles{Float64, FellMineral}(c_ref = 1.0)` |
| single-scattering albedo | `ω₀` | 0.95 | – | Fell (1997), Eq. 5.11 | `…(ω_nap = 0.9)` |
| reference wavelength | `λ_ref` | 440 | nm | Fell (1997), Eq. 5.10 | `…(λ_ref = 550)` |
| spectral exponent of attenuation | `γ` | 1 | – | Fell (1997), Eq. 5.10 | `…(spectral_exponent = 0.5)` |
| backscatter fraction | `B` | 0.036 | – | Fell (1997), §5.1.4 | `…(B = 0.02)` |

## Phase-function parameters

| model | parameter | symbol | default | source | set with |
|:--|:--|:--|:--|:--|:--|
| [`HenyeyGreensteinPhase`](@ref) | asymmetry parameter | `g` | 0.85 | package default | `HenyeyGreensteinPhase(g = 0.9)` |
| [`FournierForandPhase`](@ref) | real refractive index of the particles | `n` | from `B` | Mobley, Sundman & Boss (2002) | `FournierForandPhase{Float64}(1.10, 3.5835)` |
| [`FournierForandPhase`](@ref) | Junge slope | `μ_j` | from `B` | same | as above |
| [`FournierForandPhase`](@ref) | backscatter fraction | `B` | none (required in this form) | – | `FournierForandPhase(backscatter_fraction = 0.0183)` |

`backscatter_fraction = B` selects `(n, μ_j)` on the family
`n = 1.01 + 0.1542 (μ_j − 3)` of Mobley, Sundman & Boss (2002), a fixed coefficient
set. [`PetzoldPhase`](@ref) is tabulated data.

## Water Raman scattering: [`WaterRaman`](@ref)

`b_R(λ′) = b_ref (λ_ref/λ′)^n`; emission band Gaussian in wavenumber at `ν′ − Δν`.

| parameter | symbol | default | unit | source | set with |
|:--|:--|:--|:--|:--|:--|
| Stokes shift | `Δν` | 3400 | cm⁻¹ | attributed to Haltrin & Kattawar (1993); not in a paper checked for this release | `WaterRaman(shift_cm_inv = 3357)` |
| band standard deviation | `σ_ν` | 85 | cm⁻¹ | same | `WaterRaman(σ_cm_inv = 90)` |
| scattering coefficient at `λ_ref` | `b_ref` | 2.7·10⁻⁴ | m⁻¹ | attributed to Bartlett et al. (1998); not checked. The HydroLight runs of Loisel et al. (2023) are reproduced with 2.6·10⁻⁴ | `WaterRaman(cross_section_ref = 2.6e-4)` |
| reference wavelength | `λ_ref` | 488 | nm | same | `WaterRaman(cross_section_λ_ref = 488)` |
| exponent | `n` | 5.5 | – | attributed to Bartlett et al. (1998), energy units; not checked | `WaterRaman(cross_section_n = 5.3)` |
| depolarization ratio, linearly polarized excitation | `ρ` | 0.17 | – | Mobley et al. (1993), Eq. 9 | `WaterRaman(linear_depolarization = 0.2)` |
| excitation window | – | (350, 700) | nm | package default | `WaterRaman(excitation_range = (300, 700))` |
| convention of `b_R` | – | [`EnergyConversion`](@ref) | – | Mobley et al. (1993), Problem 7 and Appendix A | `WaterRaman(conversion = PhotonConversion())` |

The shift and width defaults are shared with [`WavenumberShiftEmission`](@ref)
(`WavenumberShiftEmission(shift_cm_inv = …, σ_cm_inv = …)`). Mobley's Problem 7 is
`WaterRaman(cross_section_ref = 0.0063, cross_section_λ_ref = 417)`.

## Fluorescence parameters

| model | parameter | symbol | default | unit | source | set with |
|:--|:--|:--|:--|:--|:--|:--|
| [`chlorophyll_fluorescence`](@ref) | quantum yield | `φ` | none (required) | – | 0.01–0.05 in situ; HydroLight 0.02 (Loisel et al. 2023) | `chlorophyll_fluorescence(phyto; quantum_yield = 0.02)` |
| | excitation window | – | (400, 700) | nm | PAR limits | `…; excitation_range = (350, 700))` |
| | carotenoid transfer efficiency | `η` | none (required with carotenoids) | – | – | `…; carotenoid_transfer = 0.5)` |
| | phycobilin transfer efficiency | `η` | 1 | – | package default | `…; phycobilin_transfer = 0.8)` |
| | photoprotective-carotenoid transfer efficiency | `η` | 0 | – | photoprotective pigments dissipate absorbed energy | `…; photoprotective_transfer = 0.1)` |
| | emission spectrum | `f` | [`chloroplast_emission`](@ref) | nm⁻¹ | Fluspect (data) | `…; emission = GaussianEmission())` |
| [`GaussianEmission`](@ref) | band centre | – | 685 | nm | Fell (1997), §2.2.9 | `GaussianEmission(center_nm = 683)` |
| | band standard deviation | `σ` | 10.6 | nm | Fell (1997), §2.2.9 | `GaussianEmission(σ_nm = 12)` |
| [`chloroplast_emission`](@ref) | photosystem I share | – | combined spectrum | – | Fluspect | `chloroplast_emission(psi_share = 0.3)` |
| [`CDOMFluorescence`](@ref) | quantum yield | `φ` | none (required) | – | – | `CDOMFluorescence(cdom; emission, quantum_yield = 0.01)` |
| | excitation window | – | (300, 400) | nm | package default | `…; excitation_range = (310, 390))` |
| [`ExcitationPathway`](@ref) | transfer efficiency | `η` | 1 | – | – | `ExcitationPathway(phyto, :carotenoids, 0.5)` |

## Absorbed-light parameters

| function | parameter | default | unit | source | set with |
|:--|:--|:--|:--|:--|:--|
| [`absorbed_par`](@ref) | integration limits | (400, 700) | nm | PAR limits; the same constant as the default excitation window of `chlorophyll_fluorescence` | `absorbed_par(o, E₀; range = (400, 690))` |
| [`inelastic_emission`](@ref) | excitation limits | `excitation_range(p)` | nm | the process | `inelastic_emission(o, j, E₀, λ; range = (370, 690))` |
| [`absorbed_power`](@ref), [`absorbed_photons`](@ref), [`absorbed_par`](@ref), [`inelastic_emission`](@ref) | negative irradiance treated as round-off, relative to the largest entry | `√eps` (1.5·10⁻⁸ in `Float64`) | – | package default | `absorbed_power(o, E₀; roundoff_rtol = 0)` |

## Numerical settings that affect results

| setting | default | set with |
|:--|:--|:--|
| number of expansion moments | `ℓ_max = 64` | `layer_optics(layer, λ; ℓ_max = 128)` |
| wavelength match of [`ExactLookup`](@ref) tables | relative `1e-6` | fixed constant of the package |
| normalization tolerance of [`TabulatedMuellerPhase`](@ref) | `1e-3` | `TabulatedMuellerPhase(…; normalization_tolerance = 1e-4)` |

Universal constants (`h`, `c`, `N_A`, SI 2019 values) are named constants.
