# Changelog

All notable changes to OceanOptics.jl are recorded here. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and
this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0] — 2026-10-16

First registered release. The entries below cover the development since
the initial phases.

### Packaging

- Requires Julia ≥ 1.10. Runtime dependencies are DocStringExtensions,
  FastGaussQuadrature and Interpolations (QuadGK, LinearAlgebra and Printf
  were unused and are dropped; QuadGK remains a test dependency). Compat
  bounds are tested at both ends (lowest versions on Julia 1.10; Julia
  1.13), and Interpolations starts at 0.14.2, the first release with
  `linear_interpolation`. TagBot and CompatHelper workflows added.

### Changed — solver-facing interface (breaking)

The interface was redesigned for direct use by vSmartMOM's coupled
atmosphere–ocean solver; there are no compatibility shims.

- `layer_optics` returns a component-based `OceanLayerOptics`: absorption
  optical depth `τ_abs` plus one `ScatteringComponent` per scattering
  constituent (label, scattering optical depth, `GreekCoefficients` once
  or per wavelength, backscatter fraction, polarization label `:scalar`,
  `:full` or `:depolarizing`). Mixed quantities come from
  `optical_depth`, `single_scattering_albedo`, `backscatter_fraction` and
  `greek_expansion(o, iλ)`. Layers are validated at construction.
- New `column_optics` → `OceanColumnOptics`: contiguous layers on one grid,
  the interface refractive index and the bottom (`BlackBottom`,
  `LambertianBottom`); both must be chosen explicitly. `interface_depths`.
- Components store their distinct expansions plus a per-wavelength
  `phase_index`, and a `phase_id`; each distinct expansion is computed
  once, and `column_optics` shares expansions between layers.
  `WavelengthDependentPhase(...; lookup = :nearest)` gives
  band-representative phase matrices on fine spectral grids.
- `GreekCoefficients` replaces the NamedTuple and matrix layouts: the
  fields, order and signs of vSmartMOM's `GreekCoefs`, validated (`β₀ = 1`,
  `|β_ℓ| ≤ 2ℓ+1`, zero `α, γ, ϵ, ζ` for `ℓ < 2`, all-or-none polarization).
  `resize_moments`, `depolarized`, `mix_greek`, `evaluate_scattering_matrix`.
- Phase interface: `greek_expansion(pf, ℓ_max[, λ]; polarized)` replaces
  `phase_matrix_moments`; scalar-only phase functions throw for
  `polarized = true` unless the completely depolarizing approximation is
  chosen (`depolarizing_greek_expansion`, or
  `layer_optics(...; scalar_phase_approximation = :depolarizing)`).
  `phase_function(c)` is now `scattering_phase(c)`. These names no
  longer collide with vSmartMOM's exports.
- Expansions are integrated in the scattering angle on graded
  Gauss-Legendre panels until converged, instead of a fixed rule in `μ`
  followed by renormalization; `phase_function_moments` throws for a
  phase function that is not normalized.
- Every interpolated table throws a `DomainError` outside its data range
  (previously clamped to the end values).
- Dependency bounds include the versions vSmartMOM resolves to
  (`FastGaussQuadrature` 0.4/0.5, `Interpolations` 0.16).

### Fixed — physics

- Fournier-Forand: the closed form has a removable 0/0 where δ = 1
  (θ ≈ 10° for n = 1.1); evaluated directly it returned values off by many
  orders of magnitude close to that angle. A series in 1 − δ is used there.
  Very small backscatter fractions now include the analytic forward tail
  below the quadrature floor; `n` is restricted to (1, 1.5).
- Fournier-Forand: the "renormalization" factor was a quadrature artifact
  (the closed form is exactly normalized) that inflated `p` and `B` by up
  to 0.3 % for small `B`; removed. The `B → (n, μ)` inversion now uses the
  family of Mobley, Sundman & Boss (2002), `n = 1.01 + 0.1542 (μ − 3)`
  (previously `μ = 3 + 20 (n − 1)`, a different forward-peak shape for the
  same `B`); `B = 0.0183` gives their `n = 1.10, μ = 3.5835`.
- Petzold: power-law (log–log) interpolation with power-law forward
  extrapolation and accurate normalization, instead of linear-in-μ
  interpolation, a flat forward peak and trapezoid normalization.
- Pure-water absorption: the default `CombinedPopeFryMason` fell back to
  Smith & Baker (4–6.5× too high below 380 nm) because the Mason, Cone &
  Fry table is not bundled; it and `MasonConeFry2016` are removed. New
  default `IOCCG2018` (the IOCCG 2018 composite used for HydroLight by
  Loisel et al. 2023, 350–750 nm, continued with Kou et al. 1993 to
  2500 nm) and `Kou1993`. The temperature correction is referenced to
  each table's measurement temperature (22 °C; none for Smith & Baker)
  instead of 20 °C.
- Temperature correction: the bundled five-point `pegau_1997_temperature.csv`
  ("condensed" digitization) missed the 740 nm O–H band, giving 0.0013
  instead of ≈ 0.016 m⁻¹ K⁻¹ there and the wrong sign near 700 nm; replaced
  by the Röttgers, McKee & Utschig (2014) coefficients (2 nm, 550–2500 nm,
  CC BY 4.0).

### Added

- Tabulated inputs: `TabulatedAbsorber` (`total = true` for bulk
  absorption; `layer_optics` then rejects other absorbers),
  `TabulatedScatterer`, `TabulatedGreekPhase`, `TabulatedMuellerPhase`
  (angular F(Θ) tables in the standard convention, monotone-cubic
  interpolation, recorded normalization) and `WavelengthDependentPhase`.
  The Chowdhary et al. (2020) AOS-IV hydrosol matrices expand to the
  paper's asymmetry parameters and backscatter ratios.
- Seawater refractive index: `FixedRefractiveIndex`, `QuanFry1995`
  (relative to air, domain-checked).
- `radiance_source_kernel` and an explicit photon/energy contract for
  inelastic processes.
- References and tests: Siewert (2000) scattering matrix reconstructed by
  vSmartMOM (`test/data/`), Chowdhary et al. (2020) Table 2, and the
  Loisel et al. (2023) HydroLight input IOPs
  (`data/reference/loisel2023_hydrolight_iops.csv`, CC0), with a docs page
  comparing them.
- Interactive Chart.js documentation figures with data tables
  (`docs/DocCharts.jl`, `docs/src/assets/oo_charts.js`); Plots.jl is no
  longer a docs dependency. Consumer API and vSmartMOM integration pages
  rewritten (with a tested adapter sketch).

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

## Initial development (not released)

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
