# OceanOptics.jl

Inherent optical properties (IOPs), phase-function Legendre expansions,
and inelastic source terms of natural waters from biogeochemical state —
in Julia, and designed to feed
[vSmartMOM.jl](https://github.com/RemoteSensingTools/vSmartMOM.jl)'s
forthcoming `OceanSurface` boundary.

`OceanOptics.jl` is the ocean counterpart to
[`CanopyOptics.jl`](https://github.com/RemoteSensingTools/CanopyOptics.jl):
both are pure data-producing packages feeding a single shared RT solver.
`OceanOptics.jl` **declares no solver dependency** — the interface it
publishes is stable and solver-agnostic.

## Design principles

The historical design memo
[`DESIGN.md`](https://github.com/RemoteSensingTools/OceanOptics.jl/blob/main/DESIGN.md)
records the original plans; where it differs from this manual, the manual is
authoritative. The four load-bearing principles are:

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

# 2. Optical properties on a wavelength grid: absorption and scattering
#    components per constituent
opt = layer_optics(layer, collect(400.0:5.0:700.0); ℓ_max = 32)
optical_depth(opt)              # τ per λ
single_scattering_albedo(opt)   # ϖ per λ
opt.scatterers[2].greek         # phytoplankton expansion (GreekCoefficients)
greek_expansion(opt, 1)         # mixed expansion at 400 nm

# 3. Chlorophyll fluorescence: chloroplast-level emission (Fluspect),
#    excitation by chlorophyll plus half of the carotenoid absorption
pig = Phytoplankton(Chl = 0.5, pigments = ProspectPigments(carotenoid_ratio = 0.3))
sif = chlorophyll_fluorescence(pig; quantum_yield = 0.02, carotenoid_transfer = 0.5)
excitation_absorption(sif, 440.0)   # φ · (a_chl + 0.5 a_car) at 440 nm
emission(sif, 440.0, 685.0)         # photons per nm at 685 nm
```

For a longer walk-through with numerical output, see
[End-to-end example](@ref).

## What's covered

* **Materials**: [`PureWater`](@ref), [`CDOM`](@ref),
  [`Phytoplankton`](@ref) (optionally with a pigment partition,
  [`ProspectPigments`](@ref) or [`Bidigare1990Pigments`](@ref)),
  [`NonAlgalParticles`](@ref). See [Materials](@ref).
* **Phase functions and matrices**: [`RayleighWaterPhase`](@ref),
  [`FournierForandPhase`](@ref), [`HenyeyGreensteinPhase`](@ref),
  [`PetzoldPhase`](@ref), [`TabulatedGreekPhase`](@ref),
  [`TabulatedMuellerPhase`](@ref), [`WavelengthDependentPhase`](@ref),
  expanded as [`GreekCoefficients`](@ref). See [Phase functions](@ref).
* **Solver output**: [`layer_optics`](@ref) → [`OceanLayerOptics`](@ref)
  (absorption, scattering and inelastic components),
  [`column_optics`](@ref) → [`OceanColumnOptics`](@ref) with the
  interface [`refractive_index`](@ref) and the floor (black or
  Lambertian). See [Consumer API](@ref) and [vSmartMOM integration](@ref).
* **Absorbed light**: [`absorbed_power`](@ref), [`absorbed_par`](@ref)
  and [`inelastic_emission`](@ref) per layer and component, from the
  solver's scalar irradiance.
* **Inelastic processes**: [`IsotropicFluorescence`](@ref) with
  [`ExcitationPathway`](@ref)s, [`chlorophyll_fluorescence`](@ref) with
  the chloroplast-level emission of Fluspect ([`chloroplast_emission`](@ref)),
  [`CDOMFluorescence`](@ref), [`WaterRaman`](@ref). See
  [Inelastic processes](@ref).
* **Layer + column**: [`OceanLayer`](@ref), [`OceanColumn`](@ref),
  `fell_column`, `uniform_column`. See [Types](@ref).
* **Data loading**: [`SpectralTable`](@ref), [`linterp`](@ref),
  [`load_spectral_csv`](@ref). See [Data loading](@ref).

## Scope and non-scope

OceanOptics.jl **does** produce IOPs, phase-function expansions, and
inelastic source kernels from biogeochemical state. The elastic path is
generic in its element type, so that ForwardDiff dual numbers can be
carried by the model parameters (chlorophyll, CDOM reference absorption
and slope, the non-algal particle parameters, temperature, salinity and
depolarization of water, the Fournier-Forand backscatter fraction, the
Henyey-Greenstein asymmetry, layer thickness) through `layer_optics` and
`column_optics`, in addition to derivatives with respect to wavelength,
scattering angle and irradiance. Inelastic processes and measured phase
tables are floating point only; see
[Forward-mode automatic differentiation](@ref) for what is and is not
differentiable. ForwardDiff is not a dependency of the package.

OceanOptics.jl **does not** solve radiative transfer. That work lives
in the consumer (e.g. `vSmartMOM.jl`), which maps the column optics onto
its own layer operators and refracting interface.
