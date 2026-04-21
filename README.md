# OceanOptics.jl

[![CI](https://github.com/RemoteSensingTools/OceanOptics.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/RemoteSensingTools/OceanOptics.jl/actions/workflows/CI.yml)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://RemoteSensingTools.github.io/OceanOptics.jl/)

Inherent optical properties (IOPs), phase-function Legendre expansions, and
inelastic source terms of natural waters from biogeochemical state — in
Julia, AD-transparent, and designed to feed
[vSmartMOM.jl](https://github.com/RemoteSensingTools/vSmartMOM.jl)'s
forthcoming `OceanSurface` boundary.

`OceanOptics.jl` is the ocean counterpart to
[CanopyOptics.jl](https://github.com/RemoteSensingTools/CanopyOptics.jl) in
the three-package vSmartMOM ecosystem: atmospheric RT couples to either
canopy or ocean as a lower boundary through a small, RT-solver-neutral
`OceanLayerOptics` bundle. The full design rationale, type-system
philosophy, and implementation roadmap live in [`DESIGN.md`](DESIGN.md).

## Status

**Phases 1 (elastic IOPs) and 2 (inelastic kernels) are complete.** The
current release ships:

- `PureWater` with choice of absorption model (Pope-Fry 1997, Smith-Baker
  1981, Mason-Cone-Fry 2016 stub, piecewise combined) and scattering
  model (Morel 1974 or Zhang-Hu 2009 — the latter a faithful port of the
  `betasw_ZHH2009.m` reference MATLAB, matching Zhang et al.'s published
  agreement with Morel measurements to 1 %).
- `CDOM` (Bricaud, Morel & Prieur 1981 exponential form).
- `Phytoplankton{Bricaud1998}` — power-law `a_φ*(λ) = A(λ)·Chl^(-E(λ))`
  with bundled 1998 Bricaud + Morrison-Nelson 2004 + Vasilkov 2005
  spectra (300–1000 nm), plus Morel & Maritorena 2001 Case-1 scattering.
- `NonAlgalParticles` with Babin et al. (2003) exponential and
  Fell (1997) mineral parameterizations.
- Phase functions: `RayleighWaterPhase`, `FournierForandPhase`
  (Mobley 2002 B → n,μ_j inversion + forward-peak renormalization),
  `HenyeyGreensteinPhase` (closed-form moments), and `PetzoldPhase`
  (San Diego Harbor VSF from the Ocean Optics Web Book).
- `layer_optics(layer, λ_grid)` — the main entry point, producing a
  complete `OceanLayerOptics` bundle of τ, ϖ, B, and β_ℓ(λ).
- `AbstractOceanInelasticProcess` + trait functions for
  `IsotropicFluorescence`, `CDOMFluorescence`, `WaterRaman`, plus
  `chlorophyll_fluorescence(phyto)` with Fell (1997) defaults.

Phase 3 (coupled RT with `OceanSurface <: AbstractSurfaceType`) lives
in `vSmartMOM.jl`, not here. Phase 4 (polarized Greek-coefficient
expansion, `OceanColumn` + Fell-grid helper, artifact-based data
loading, JOSS paper) is the remaining in-repo scope. See DESIGN §10.

## Install

```julia
julia> ]
(@v1.x) pkg> add https://github.com/RemoteSensingTools/OceanOptics.jl
```

## Quick start

```julia
using OceanOptics

# Open-ocean Case-1 water, moderate chlorophyll.
water = PureWater{Float64}(
    absorption_model = CombinedPopeFryMason(),   # HydroLight-style piecewise
    scattering_model = ZhangHu2009(),            # faithful ZHH-2009 port
    temperature      = 20.0,
    salinity         = 35.0,
)
cdom  = CDOM{Float64}(a_ref = 0.03, λ_ref = 440.0, slope = 0.014)
phyto = Phytoplankton{Float64, Bricaud1998}(Chl = 0.5)

layer = OceanLayer(0.0, 5.0,
    AbstractOceanConstituent{Float64}[water, cdom, phyto])

λ   = collect(400.0:5.0:700.0)
opt = layer_optics(layer, λ; ℓ_max = 64)

opt.τ        # per-λ optical depth           Vector{Float64}
opt.ϖ        # per-λ single-scattering albedo
opt.β        # Legendre moments β_ℓ(λ)       Matrix (ℓ_max+1, nλ)

# Inelastic kernels for a downstream two-pass RT solver
sif = chlorophyll_fluorescence(phyto)          # Fell 1997 defaults
excitation_absorption(sif, 440.0)              # φ · a_φ(440)
emission(sif, 440.0, 685.0)                    # 1/nm, ~0.0376
```

`OceanLayerOptics` is the RT-solver-neutral interface: any plane-parallel
ocean RT code (vSmartMOM, HydroLight-style ODE, Monte Carlo, SOS) can
consume the same output.

## Consumer API (what an RT solver consumes)

OceanOptics.jl publishes a stable, RT-agnostic interface. Downstream
solvers (vSmartMOM, HydroLight-style ODE codes, Monte-Carlo,
SOS, …) consume it *without* OceanOptics having any knowledge of
the solver. `Project.toml` therefore lists no RT-solver package
under `[deps]` or `[weakdeps]`; adapters live on the solver side.

The stable contract is:

* **`OceanLayerOptics{FT}`** — per-layer bundle with fields `λ`
  (`[nm]`), `Δz` (`[m]`), `τ`, `ϖ`, `B`, and `β` (Legendre moments,
  indexed `[ℓ+1, λ]`). Optional polarized Greek elements α, γ, δ,
  ϵ, ζ follow the same `[ℓ+1, λ]` layout when populated.
* **`layer_optics(layer::OceanLayer, λ_grid; ℓ_max)`** — the
  single entry point that produces an `OceanLayerOptics` from a
  biogeochemical state.
* **`AbstractOceanInelasticProcess{FT}`** hierarchy and its trait
  functions `excitation_absorption(p, λ')`, `emission(p, λ', λ)`,
  `excitation_range(p)`, `is_isotropic(p)`,
  `inelastic_coefficient(p, λ', λ)` (DESIGN §7.1). These are the
  Fell (1997) factorization kernels a two-pass RT solver multiplies
  against its pass-1 scalar irradiance `E°(τ, λ')`.
* **`OceanLayer{FT}.fluorophores::Vector{<:AbstractOceanInelasticProcess{FT}}`**
  — per-layer list of active inelastic processes.

Any API surface outside the above list is internal and may change.

## Automatic differentiation

Every forward function on the elastic path is AD-transparent via
ForwardDiff (see DESIGN §12.5 and `test/test_ad.jl`):

```julia
using ForwardDiff, OceanOptics

w = PureWater{Float64}(absorption_model = PopeFry1997())
ForwardDiff.derivative(λ -> absorption(w, λ), 440.0)   # ∂a_w/∂λ at 440 nm
```

Boundary handling uses `Interpolations.jl`'s `Flat()` mode, whose
`clamp(x, x_min, x_max)` mechanism keeps `Dual` inputs intact with zero
partial derivative outside the tabulated range — matching the
derivative of a constant extrapolation.

## Data bundling

Spectral reference tables live in `data/` as provenance-commented CSVs.
Every file records its upstream citation, DOI, URL, retrieval date, and
any unit conversion applied at load time (DESIGN §3.4). Currently bundled:

| File                           | Source (license)                             | What                                      |
| ------------------------------ | -------------------------------------------- | ----------------------------------------- |
| `pope_fry_1997.csv`            | omlc.org (public)                            | Pure water absorption, 380–727.5 nm       |
| `smith_baker_1981.csv`         | omlc.org (public)                            | Pure water absorption, 200–800 nm         |
| `pegau_1997_temperature.csv`   | Pegau et al. 1997, digitized                 | Temperature derivative ∂a/∂T, 600–800 nm  |
| `bricaud_1998_A.csv`           | Ocean Optics Web Book (CC-BY), Bricaud 1998 + Morrison & Nelson 2004 + Vasilkov 2005 | Phytoplankton `A(λ)`, 300–1000 nm |
| `bricaud_1998_E.csv`           | Ocean Optics Web Book (CC-BY), Bricaud 1998 + UV extensions | Phytoplankton `E(λ)`, 300–1000 nm |
| `petzold_sdh.csv`              | Ocean Optics Web Book (CC-BY) / Mobley 1994 Table 3.10 | Petzold "San Diego Harbor" VSF, 0.1°–180° |

Re-fetch from upstream at any time with:

```shell
julia --project=. scripts/fetch_data.jl
```

### Still stubbed (paywalled sources)

- `mason_cone_fry_2016.csv` — UV-extension to the Pope-Fry table,
  Appl. Opt. 55, 7163 Table 2. No redistributable source found.
  Selecting `MasonConeFry2016` raises a helpful error pointing you at
  the paper; drop a transcribed CSV into `data/` to activate.
  `CombinedPopeFryMason` falls back to Smith-Baker below 380 nm in the
  meantime.

The Bricaud bundle is the **1998 update** to the original 1995 JGR fit;
the 1995 table is paywalled and the 1998 + UV-extended values are the
community-standard replacement (they are what HydroLight 5 and NASA
OCSSW ship).

## Testing

```shell
julia --project=. -e 'using Pkg; Pkg.test()'
```

276 tests across IOP algebra, material spectra, phase-function
normalization and moments, layer assembly, inelastic kernels
(fluorescence + Raman), and forward-mode AD.

## Related packages

- [vSmartMOM.jl](https://github.com/RemoteSensingTools/vSmartMOM.jl) —
  atmospheric RT solver (Caltech/JPL). The primary consumer of
  OceanOptics.jl's output.
- [CanopyOptics.jl](https://github.com/RemoteSensingTools/CanopyOptics.jl)
  — land-surface IOPs (sibling package in the vSmartMOM ecosystem).

## License

Apache License, Version 2.0 — see [`LICENSE`](LICENSE). Matches the
sibling packages in the RemoteSensingTools ecosystem
(`vSmartMOM.jl`, `CanopyOptics.jl`).

### Third-party reference data

Bundled reference tables in `data/` retain their upstream licenses:

- `pope_fry_1997.csv`, `smith_baker_1981.csv` — public-domain data
  re-hosted by OMLC.
- `bricaud_1998_A.csv`, `bricaud_1998_E.csv`, `petzold_sdh.csv` —
  redistributed from the Ocean Optics Web Book (Mobley, ongoing)
  under a Creative Commons Attribution license; cite the original
  papers (Bricaud et al. 1998, Morrison & Nelson 2004, Vasilkov et
  al. 2005, Petzold 1972, Mobley 1994) in derived work.
- `pegau_1997_temperature.csv` — digitized from Pegau, Gray &
  Zaneveld (1997) Appl. Opt. 36, 6035.
