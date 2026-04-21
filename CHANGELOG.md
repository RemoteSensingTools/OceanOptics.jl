# Changelog

All notable changes to OceanOptics.jl are recorded here. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and
this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added — Phase 4 (in progress)

- `OceanColumn{FT}` — ordered, surface-anchored stack of `OceanLayer`s
  with contiguity validation (DESIGN §4.5). `n_layers`, `thickness`,
  `bottom_depth`, `midpoint_depths`, `layer_thicknesses` accessors plus
  `Base.length` / `getindex` / `iterate`.
- `uniform_column(constituents; n_layers, depth, fluorophores)` —
  equal-thickness stack.
- `fell_column(constituents; depth, fluorophores)` — Fell (1997) §5.2.1
  grid (1 m in 0–10 m, 2 m in 10–20 m, 5 m below), with graceful
  truncation for `depth < 20 m`.
- Apache License, Version 2.0 (matches sibling `vSmartMOM.jl` and
  `CanopyOptics.jl`).
- CHANGELOG.

### Added — Phase 1 polish

- Faithful Zhang, Hu & He (2009) seawater-scattering port of the
  reference MATLAB `betasw_ZHH2009.m` — density + concentration
  fluctuations with UNESCO 1981 EOS, Millero (1980) secant-bulk
  compressibility, Quan-Fry 1995 refractive index, Ciddor 1996 air
  correction, Millero-Leung 1976 water activity, PMH density
  derivative, and the Cabannes `b = (8π/3)·β(90°)·(2+δ)/(1+δ)` closure.
  Replaces the earlier calibrated power-law approximation.
- `data/petzold_sdh.csv` — Petzold "San Diego Harbor" VSF, CC-BY via
  Ocean Optics Web Book (from Mobley 1994 Appendix A.4 / Table 3.10).
  `PetzoldPhase{FT}()` now constructs and normalizes without error.
- `data/bricaud_1998_A.csv`, `data/bricaud_1998_E.csv` — Bricaud (1998)
  `A(λ)` / `E(λ)` phytoplankton spectra with Morrison & Nelson (2004)
  and Vasilkov et al. (2005) UV extensions, 300–1000 nm. CC-BY via
  Ocean Optics Web Book.

### Changed

- Renamed `Bricaud1995` → `Bricaud1998` to match the actually-bundled
  data. No back-compat alias (pre-1.0 package).
- `Phytoplankton{FT, Bricaud1998}.absorption(λ)` now returns real
  numerical values from the bundled spectra (previously errored with a
  "data not bundled" message).

### Removed

- `to_core_scattering_properties` stub from `src/mixing/layer_assembly.jl`.
  Consumer-side adapters live in the consumer — DESIGN §3.1 rewritten
  accordingly, `[weakdeps]` on vSmartMOM dropped. See
  [No solver deps feedback](#) for the architectural decision.

## [0.1.0] — initial development

### Added — Phase 2 (inelastic framework)

- `AbstractOceanInelasticProcess{FT}` hierarchy with kernel traits
  `excitation_absorption(p, λ')`, `emission(p, λ', λ)`,
  `excitation_range(p)`, `is_isotropic(p)`,
  `inelastic_coefficient(p, λ', λ)` — the Fell (1997) §2.2.9
  factorization primitives that a two-pass RT solver consumes against
  a pass-1 scalar irradiance field.
- `AbstractSpectralRedistribution{FT}` with `GaussianEmission`
  (fixed-wavelength Gaussian — chlorophyll / CDOM) and
  `WavenumberShiftEmission` (Stokes-shifted Gaussian in wavenumber
  with `|dν/dλ|` Jacobian — water Raman).
- `IsotropicFluorescence{FT, A, E}`, `CDOMFluorescence{FT, E}`,
  `WaterRaman{FT}` concrete processes. `chlorophyll_fluorescence(phyto)`
  convenience constructor with Fell (1997) defaults (685 nm peak,
  σ = 10.6 nm, φ = 0.003).
- `OceanLayer.fluorophores::Vector{<:AbstractOceanInelasticProcess}`
  field (backward-compatible default empty).

### Added — Phase 1 (elastic path)

- Type spine: `AbstractOceanConstituent`, `AbstractAbsorber`,
  `AbstractScatterer`, `AbstractAbsorbingScatterer`, plus
  model-tag hierarchies (`AbstractPureWaterAbsorptionModel`,
  `AbstractPureWaterScatteringModel`, `AbstractPhytoplanktonModel`,
  `AbstractCDOMModel`, `AbstractNAPModel`).
- `IOP{FT, A}` with `+`/`*` algebra and `zero(iop)`.
- `OceanLayer{FT}` with contiguity-validated constructor.
- `OceanLayerOptics{FT, V, M}` — the stable RT-solver-neutral output
  bundle of `(λ, Δz, τ, ϖ, B, β)` plus optional polarized Greek
  elements `(α, γ, δ, ϵ, ζ)`.
- `layer_optics(layer, λ_grid; ℓ_max)` — the single forward entry
  point producing an `OceanLayerOptics`.
- Materials: `PureWater` (Pope-Fry / Smith-Baker / Mason-Cone-Fry
  stub / `CombinedPopeFryMason` piecewise, plus Morel 1974 or
  Zhang-Hu 2009 scattering), `CDOM` (Bricaud-Morel-Prieur 1981
  exponential), `Phytoplankton`, `NonAlgalParticles` with
  `BabinBricaud2003` and `FellMineral` parameterizations.
- Phase functions: `RayleighWaterPhase` (closed-form β₂ = (1-σ)/(5(2+σ))),
  `FournierForandPhase` (Mobley 2002 B → n,μ_j inversion + adaptive
  renormalization of the forward peak),
  `HenyeyGreensteinPhase` (closed-form β_ℓ = g^ℓ and analytic
  backscatter fraction), `PetzoldPhase` (stub on initial ship, data
  bundled in Phase 1 polish).
- Spectral data loader (`SpectralTable`, `linterp`, `load_spectral_csv`,
  `data_path`) backed by `Interpolations.linear_interpolation` with
  `Flat()` extrapolation — AD-transparent via `clamp`.
- Bundled reference data with provenance headers:
  `pope_fry_1997.csv`, `smith_baker_1981.csv`,
  `pegau_1997_temperature.csv`. Reproducible re-fetch via
  `scripts/fetch_data.jl`.
- Tests: 270 passing across IOP algebra, material spectra,
  phase-function normalization and moments, layer assembly, inelastic
  kernels, and forward-mode AD (`test_ad.jl`, DESIGN §12.5).
- CI: GitHub Actions, Julia 1.9 + stable, Linux + macOS, x64 + aarch64.

[Unreleased]: https://github.com/RemoteSensingTools/OceanOptics.jl/compare/main
[0.1.0]: https://github.com/RemoteSensingTools/OceanOptics.jl/commit/3a9e2d3
