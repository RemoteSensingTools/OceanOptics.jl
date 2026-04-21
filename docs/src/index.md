# OceanOptics.jl

Inherent optical properties (IOPs), phase-function Legendre expansions,
and inelastic source terms of natural waters from biogeochemical state —
in Julia, automatic-differentiation-transparent, and designed to feed
[vSmartMOM.jl](https://github.com/RemoteSensingTools/vSmartMOM.jl)'s
forthcoming `OceanSurface` boundary.

`OceanOptics.jl` is the ocean counterpart to
[`CanopyOptics.jl`](https://github.com/RemoteSensingTools/CanopyOptics.jl):
both are pure data-producing packages feeding a single shared RT solver.
`OceanOptics.jl` **declares no solver dependency** — the interface it
publishes is stable and solver-agnostic.

## Design principles

The full architecture memo lives in
[`DESIGN.md`](https://github.com/RemoteSensingTools/OceanOptics.jl/blob/main/DESIGN.md)
and is the authoritative long-form document. The four load-bearing
principles are:

1. **RT-solver neutrality.** The primary output is an
   [`OceanLayerOptics`](@ref) bundle consumable by any plane-parallel
   solver. No RT code lives here.
2. **Three orthogonal dispatch axes.** `(optical role) × (biogeochemical
   identity) × (parameterization tag)` each get their own slot in the
   type system so new parameterizations are additive changes.
3. **Physical quantities with algebra.** `IOP{FT}` mixes by `+`;
   phase-function moments mix by scattering-weighted averages; layer
   optics compose vertically.
4. **Data provenance as a first-class concern.** Every spectral table
   in `data/` carries its upstream citation, DOI, URL, retrieval date,
   and unit-conversion record as a header comment.

## Install

```julia
julia> ]
(@v1.x) pkg> add OceanOptics
```

Until registered in the General registry, install from source:

```julia
(@v1.x) pkg> add https://github.com/RemoteSensingTools/OceanOptics.jl
```

## 60-second tour

```julia
using OceanOptics

# 1. A Case-1 water layer at Chl = 0.5 mg/m³
water = PureWater{Float64}(scattering_model = ZhangHu2009())
cdom  = CDOM{Float64}(a_ref = 0.03)
phyto = Phytoplankton{Float64, Bricaud1998}(Chl = 0.5)
layer = OceanLayer(0.0, 5.0,
                   AbstractOceanConstituent{Float64}[water, cdom, phyto])

# 2. Compute the elastic bundle over a wavelength grid
opt = layer_optics(layer, collect(400.0:5.0:700.0); ℓ_max = 32)
opt.τ   # optical depth per λ
opt.ϖ   # single-scattering albedo per λ
opt.β   # Legendre moments β_ℓ(λ), size (ℓ_max+1, n_λ)

# 3. Inelastic kernels for a downstream two-pass RT solver
sif = chlorophyll_fluorescence(phyto)
excitation_absorption(sif, 440.0)   # φ · a_φ(440)
emission(sif, 440.0, 685.0)         # Gaussian at 685 nm, σ = 10.6 nm
```

For a longer walk-through with numerical output, see
[End-to-end example](@ref).

## What's covered

* **Materials**: [`PureWater`](@ref), [`CDOM`](@ref),
  [`Phytoplankton`](@ref), [`NonAlgalParticles`](@ref). See
  [Materials](@ref).
* **Phase functions**: [`RayleighWaterPhase`](@ref),
  [`FournierForandPhase`](@ref), [`HenyeyGreensteinPhase`](@ref),
  [`PetzoldPhase`](@ref). See [Phase functions](@ref).
* **Inelastic processes**: [`IsotropicFluorescence`](@ref),
  [`CDOMFluorescence`](@ref), [`WaterRaman`](@ref). See
  [Inelastic](@ref).
* **Layer + column**: [`OceanLayer`](@ref), [`OceanColumn`](@ref),
  `fell_column`, `uniform_column`. See [Types](@ref).
* **Data loading**: [`SpectralTable`](@ref), [`linterp`](@ref),
  [`load_spectral_csv`](@ref). See [Data loading](@ref).

## Scope and non-scope

OceanOptics.jl **does** produce IOPs, phase-function expansions, and
inelastic source kernels from biogeochemical state. It is AD-transparent
through ForwardDiff.

OceanOptics.jl **does not** solve radiative transfer. That work lives
in the consumer (e.g. `vSmartMOM.jl`) which imports our types and
builds its own `OceanSurface <: AbstractSurfaceType`.
