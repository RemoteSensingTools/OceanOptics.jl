# OceanOptics.jl

[![CI](https://github.com/RemoteSensingTools/OceanOptics.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/RemoteSensingTools/OceanOptics.jl/actions/workflows/CI.yml)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://RemoteSensingTools.github.io/OceanOptics.jl/)

Inherent optical properties (IOPs) and phase-matrix expansions of natural
waters, for plane-parallel radiative-transfer solvers — primarily the
coupled atmosphere–ocean solver of
[vSmartMOM.jl](https://github.com/RemoteSensingTools/vSmartMOM.jl), with
no dependency on it. It is the ocean counterpart of
[CanopyOptics.jl](https://github.com/RemoteSensingTools/CanopyOptics.jl).

## What it provides

- **Constituents:** `PureWater` (absorption `IOCCG2018` (default),
  `PopeFry1997`, `Kou1993`, `SmithBaker1981`; scattering `ZhangHu2009`
  (default) or `Morel1974`), `CDOM`, `Phytoplankton{Bricaud1998}`
  (optionally split into pigment groups: `ProspectPigments`,
  `Bidigare1990Pigments`), `NonAlgalParticles`, and exact published
  values (`TabulatedAbsorber`, with a `total` flag for bulk absorption;
  `TabulatedScatterer`).
- **Phase functions and matrices**, expanded as `GreekCoefficients` in
  vSmartMOM's convention: `RayleighWaterPhase` (phase matrix),
  `FournierForandPhase` (Mobley et al. 2002 B → (n, μ) family),
  `HenyeyGreensteinPhase`, `PetzoldPhase`, `TabulatedGreekPhase`,
  `TabulatedMuellerPhase` (angular F(Θ) tables, e.g. benchmark
  hydrosols) and `WavelengthDependentPhase`. Expansions are integrated in
  scattering angle until converged; nothing is renormalized silently.
- **Solver output:** `layer_optics` → `OceanLayerOptics` (absorption
  components per constituent or pigment group; scattering components,
  each with its phase-matrix expansion and a polarization label; inelastic
  processes linked to the absorption that excites them) and
  `column_optics` → `OceanColumnOptics` (layers, interface refractive
  index, black or Lambertian floor).
- **Seawater refractive index:** `FixedRefractiveIndex`, `QuanFry1995`.
- **Fluorescence and Raman:** chlorophyll fluorescence with the
  chloroplast-level emission of the leaf model Fluspect, excited by
  chlorophyll and a chosen share of carotenoid absorption; water Raman;
  an explicit photon/energy contract (`radiance_source_kernel`).
- **Absorbed light:** absorbed power, photons and PAR per layer and
  component (e.g. by chlorophyll), and the fluorescence emission it
  drives, from the solver's scalar irradiance.

Every table throws outside its data range, scalar-only phase functions are
never given an invented phase matrix, and published bulk absorption cannot
be combined with separately modelled water or CDOM absorption.

## Install

```julia
pkg> add OceanOptics
```

Julia 1.10 or later.

## New in 0.2.0

- Water Raman is an energy conversion by default (no λ′/λ factor); its phase function follows Mobley et al. (1993) Eq. 9 with ρ = 0.17.
- Every physical parameter is a field or keyword with a documented default ([parameter reference](https://RemoteSensingTools.github.io/OceanOptics.jl/stable/reference/parameters/)).
- Derivatives with respect to model parameters (Chl, CDOM slope, temperature, salinity, Fournier-Forand B, …) by ForwardDiff, through `layer_optics` and `column_optics`.
- Pigment partitions hold their edge shares below 400 nm; photoprotective carotenoids no longer excite chlorophyll fluorescence by default.
- Several breaking changes (`WaterRaman`, `Phytoplankton.scattering_model`, `Gordon1992` removed): see the CHANGELOG and the "Migrating from 0.1" section of the Consumer API page.

## Quick start

```julia
using OceanOptics

water = PureWater{Float64}(temperature = 18.0, salinity = 35.0)
cdom  = CDOM{Float64}(a_ref = 0.03, λ_ref = 440.0, slope = 0.014)
phyto = Phytoplankton{Float64, Bricaud1998}(Chl = 0.5)
column = uniform_column(AbstractOceanConstituent{Float64}[water, cdom, phyto];
                        n_layers = 4, depth = 40.0)

col = column_optics(column, collect(400.0:5.0:700.0);
                    refractive_index = QuanFry1995(temperature = 18.0, salinity = 35.0),
                    bottom = BlackBottom(), ℓ_max = 256, polarized = true,
                    scalar_phase_approximation = :depolarizing)

layer = col[1]
optical_depth(layer), single_scattering_albedo(layer)
[(c.label, c.polarization) for c in layer.scatterers]   # PureWater :full, Phytoplankton :depolarizing
greek_expansion(layer.scatterers[2], 1)                # GreekCoefficients of the particles
```

Chlorophyll fluorescence, excited by chlorophyll and half of the
carotenoid absorption, and the PAR each pigment group absorbs:

```julia
phyto = Phytoplankton(Chl = 0.5, pigments = ProspectPigments(carotenoid_ratio = 0.3))
sif   = chlorophyll_fluorescence(phyto; quantum_yield = 0.02, carotenoid_transfer = 0.5)
layer = OceanLayer(0.0, 5.0, AbstractOceanConstituent{Float64}[water, phyto];
                   inelastic = AbstractOceanInelasticProcess{Float64}[sif])
o  = layer_optics(layer, collect(400.0:5.0:750.0))
E₀ = ones(n_wavelengths(o))   # stand-in for the solver's layer-mean scalar irradiance [W m⁻² nm⁻¹]
absorbed_par(o, E₀)           # µmol photons m⁻² s⁻¹ per absorption component
```

## Precision and devices

The results have the element type of the layer or column (`Float64` above);
pass another as the first argument. They are host arrays, and `Adapt.adapt`
moves their numeric arrays to a device without changing the element type:

```julia
col32 = column_optics(Float32, column, collect(400.0:5.0:700.0);
                      refractive_index = FixedRefractiveIndex(1.34),
                      bottom = BlackBottom())
# using Adapt, CUDA; adapt(CuArray, col32)
c = col32[1].scatterers[1]
phase_id(c, 1) == phase_id(greek_expansion(c, 1))   # true
```

`phase_id` is a content key of an expansion (a 64-bit hash of its element
type and coefficients): equal keys mean equal coefficients in any layer,
column or session, so a solver can build phase operators once per key.

## Documentation

The [documentation](https://RemoteSensingTools.github.io/OceanOptics.jl/)
covers the consumer contract (units, conventions, validation), the
[parameter reference](https://RemoteSensingTools.github.io/OceanOptics.jl/stable/reference/parameters/)
(every physical parameter with its default, unit and source), the
vSmartMOM mapping with a tested adapter sketch, a comparison with the
HydroLight IOP database of Loisel et al. (2023), and interactive figures.

## Testing

```shell
julia --project=. -e 'using Pkg; Pkg.test()'
```

Tests check against independent references where they exist: closed
forms (Rayleigh, Henyey-Greenstein, Fournier-Forand normalization and
backscatter), the published Mobley et al. (2002) Fournier-Forand pair, a
vSmartMOM-generated Siewert (2000) scattering matrix, Chowdhary et al.
(2020) Table 2 ocean bodies, and the Loisel et al. (2023) HydroLight
inputs.

## Data

Spectral tables live in `data/` as provenance-commented CSVs (citation,
DOI, URL, retrieval date, unit changes):

| File | Source (license) | What |
|---|---|---|
| `ioccg_2018_pure_water.csv` | Loisel et al. 2023 HydroLight database (CC0) | IOCCG (2018) recommended pure-water absorption, 350–750 nm |
| `kou_1993.csv` | omlc.org | Liquid-water absorption, 667–2500 nm (Kou et al. 1993) |
| `pope_fry_1997.csv` | omlc.org | Pure-water absorption, 380–727.5 nm |
| `smith_baker_1981.csv` | omlc.org | Clearest natural waters, 200–800 nm |
| `rottgers_2014_temperature.csv` | Röttgers et al. 2014, figshare (CC BY 4.0) | Temperature coefficient Ψ_T, 550–2500 nm |
| `bricaud_1998_A.csv`, `bricaud_1998_E.csv` | Ocean Optics Web Book (CC-BY) | Particulate absorption `a_p = A(λ)·Chl^E(λ)` (phytoplankton plus detritus; HydroLight Case-1 tables), 300–1000 nm |
| `petzold_sdh.csv` | Ocean Optics Web Book (CC-BY) | Petzold San Diego Harbor VSF, 0.1°–180° |
| `bidigare_1990_pigments.csv` | Ocean Optics Web Book (CC-BY) | In vivo specific absorption of algal pigment groups, 400–750 nm (Bidigare et al. 1990) |
| `fluspect_chloroplast_emission.csv` | SCOPE parameter file `Optipar2017_ProspectD` (used there by Fluspect, Vilfan et al. 2016, 2018; original measurements not traced; data only, no SCOPE code) | Chloroplast-level fluorescence emission (combined, PSI, PSII), 640–851 nm |
| `prospect_d_pigment_absorption.csv` | same SCOPE parameter file, attributed there to PROSPECT-D (Féret et al. 2017); not compared with the authors' own distribution; data only | PROSPECT-D in vivo chlorophyll and carotenoid absorption, 400–800 nm |
| `reference/loisel2023_hydrolight_iops.csv` | Loisel et al. 2023 (CC0) | HydroLight input IOPs: pure seawater and three scenarios |

Cite the original papers in derived work.

## License

Apache License, Version 2.0 — see [`LICENSE`](LICENSE). Bundled
reference data keep their upstream licenses (above).
