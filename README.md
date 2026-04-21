# OceanOptics.jl

[![CI](https://github.com/RemoteSensingTools/OceanOptics.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/RemoteSensingTools/OceanOptics.jl/actions/workflows/CI.yml)

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

**Phase 1 (elastic path, single layer) — under active development.** The
current release ships:

- `PureWater` with choice of absorption model (Pope-Fry 1997,
  Smith-Baker 1981, Mason-Cone-Fry 2016, piecewise combined) and
  scattering model (Morel 1974, Zhang-Hu-like 2009 calibrated approximation).
- `CDOM` (Bricaud, Morel & Prieur 1981 exponential form).
- `Phytoplankton` with Bricaud (1995) absorption (data-dependent; see
  [data bundling notes](#data-bundling)) and Morel & Maritorena (2001)
  Case-1 scattering.
- `NonAlgalParticles` with both Babin et al. (2003) exponential and
  Fell (1997) mineral parameterizations.
- Phase functions: `RayleighWaterPhase`, `FournierForandPhase` with
  Mobley (2002) B→(n,μ_j) inversion and forward-peak renormalization,
  `HenyeyGreensteinPhase` (closed-form moments), and `PetzoldPhase`
  (data-dependent).
- `layer_optics(layer, λ_grid)` — the main entry point, producing a
  complete `OceanLayerOptics` bundle of τ, ϖ, B, and β_ℓ(λ).

Phases 2 (inelastic fluorescence / Raman), 3 (coupled RT through a
`vSmartMOM.OceanSurface`), and 4 (polarized Greek expansion, artifact
data loading, JOSS paper) are planned. See DESIGN §10.

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
    scattering_model = ZhangHu2009(),            # calibrated approximation
    temperature      = 20.0,
    salinity         = 35.0,
)
cdom  = CDOM{Float64}(a_ref = 0.03, λ_ref = 440.0, slope = 0.014)
# Bricaud 1995 absorption requires the data CSVs (see "Data bundling").
# Until those are provided, scattering-only phytoplankton still works:
phyto = Phytoplankton{Float64, Bricaud1995}(Chl = 0.5)

layer = OceanLayer(0.0, 5.0,
    AbstractOceanConstituent{Float64}[water, cdom, phyto])

λ   = collect(400.0:5.0:700.0)
opt = layer_optics(layer, λ; ℓ_max = 64)

opt.τ        # per-λ optical depth           Vector{Float64}
opt.ϖ        # per-λ single-scattering albedo
opt.β        # Legendre moments β_ℓ(λ)       Matrix (ℓ_max+1, nλ)
```

`OceanLayerOptics` is the RT-solver-neutral interface: any plane-parallel
ocean RT code (vSmartMOM, HydroLight-style ODE, Monte Carlo, SOS) can
consume the same output.

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

| File                           | Source                          | What                                      |
| ------------------------------ | ------------------------------- | ----------------------------------------- |
| `pope_fry_1997.csv`            | omlc.org                        | Pure water absorption, 380–727.5 nm       |
| `smith_baker_1981.csv`         | omlc.org                        | Pure water absorption, 200–800 nm         |
| `pegau_1997_temperature.csv`   | Pegau et al. 1997, digitized    | Temperature derivative ∂a/∂T, 600–800 nm  |

Re-fetch from upstream at any time with:

```shell
julia --project=. scripts/fetch_data.jl
```

### Not yet bundled (transcription required)

- `mason_cone_fry_2016.csv` — UV extension, Appl. Opt. 55, 7163 Table 2.
- `bricaud_1995_A.csv`, `bricaud_1995_E.csv` — phytoplankton `A(λ)` and
  `E(λ)` spectra, JGR 100, 13321 Table 2.
- `petzold_sdh.csv` — Petzold "San Diego Harbor" VSF, Mobley (1994)
  *Light and Water* Appendix A.4.

Selecting a model that depends on one of these files (e.g.
`MasonConeFry2016`, `Phytoplankton{Bricaud1995}.absorption`) raises an
error that points you to the paper; drop the CSV into `data/` to
activate. Mason, Bricaud, and Petzold are on the Phase 4 to-do list.

## Testing

```shell
julia --project=. -e 'using Pkg; Pkg.test()'
```

199 tests across IOP algebra, material spectra, phase-function
normalization and moments, layer assembly, and forward-mode AD.

## Related packages

- [vSmartMOM.jl](https://github.com/RemoteSensingTools/vSmartMOM.jl) —
  atmospheric RT solver (Caltech/JPL). The primary consumer of
  OceanOptics.jl's output.
- [CanopyOptics.jl](https://github.com/RemoteSensingTools/CanopyOptics.jl)
  — land-surface IOPs (sibling package in the vSmartMOM ecosystem).

## License

[Add license before registering.]
