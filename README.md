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
  (default) or `Morel1974`), `CDOM`, `Phytoplankton{Bricaud1998}`,
  `NonAlgalParticles`, and exact published values (`TabulatedAbsorber`,
  with a `total` flag for bulk absorption; `TabulatedScatterer`).
- **Phase functions and matrices**, expanded as `GreekCoefficients` in
  vSmartMOM's convention: `RayleighWaterPhase` (phase matrix),
  `FournierForandPhase` (Mobley et al. 2002 B → (n, μ) family),
  `HenyeyGreensteinPhase`, `PetzoldPhase`, `TabulatedGreekPhase`,
  `TabulatedMuellerPhase` (angular F(Θ) tables, e.g. benchmark
  hydrosols) and `WavelengthDependentPhase`. Expansions are integrated in
  scattering angle until converged; nothing is renormalized silently.
- **Solver output:** `layer_optics` → `OceanLayerOptics` (absorption
  optical depth plus one scattering component per constituent, each with
  its phase-matrix expansion and a polarization label) and
  `column_optics` → `OceanColumnOptics` (layers, interface refractive
  index, bottom).
- **Seawater refractive index:** `FixedRefractiveIndex`, `QuanFry1995`.
- **Inelastic kernels** for fluorescence and water Raman, with an
  explicit photon/energy contract (`radiance_source_kernel`).

Every table throws outside its data range, scalar-only phase functions are
never given an invented phase matrix, and published bulk absorption cannot
be combined with separately modelled water or CDOM absorption.

## Install

```julia
pkg> add OceanOptics
```

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
[(c.label, c.polarization) for c in layer.components]   # PureWater :full, Phytoplankton :depolarizing
greek_expansion(layer.components[2], 1)                # GreekCoefficients of the particles
```

## Documentation

The [documentation](https://RemoteSensingTools.github.io/OceanOptics.jl/)
covers the consumer contract (units, conventions, validation), the
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
| `bricaud_1998_A.csv`, `bricaud_1998_E.csv` | Ocean Optics Web Book (CC-BY) | Phytoplankton `A(λ)`, `E(λ)`, 300–1000 nm |
| `petzold_sdh.csv` | Ocean Optics Web Book (CC-BY) | Petzold San Diego Harbor VSF, 0.1°–180° |
| `reference/loisel2023_hydrolight_iops.csv` | Loisel et al. 2023 (CC0) | HydroLight input IOPs: pure seawater and three scenarios |

Cite the original papers in derived work.

## License

Apache License, Version 2.0 — see [`LICENSE`](LICENSE). Bundled
reference data keep their upstream licenses (above).
