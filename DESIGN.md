# OceanOptics.jl — Design Memo

Version 0.1 · Draft · April 2026

---

## 1. Scope and audience

This memo specifies the design of **OceanOptics.jl**, a Julia package
computing the inherent optical properties (IOPs), phase-function
expansions, and inelastic source terms of natural waters from
biogeochemical state.

The audience is anyone picking this up to implement, extend, or
review: a Julia developer familiar with multiple dispatch, the
general shape of atmospheric RT, and curious about ocean optics but
not necessarily fluent in it. All external references used below
are linked; no prior familiarity with the ocean-optics literature is
assumed.

The memo is deliberately organized as a design document in the
CliMA/Oceananigans tradition: the abstractions come first, the
implementation second, with extensive cross-references to the
ecosystem packages this one lives in.

---

## 2. Context

OceanOptics.jl occupies a specific niche in a three-package
ecosystem built around [**vSmartMOM.jl**](https://github.com/RemoteSensingTools/vSmartMOM.jl),
the Caltech/JPL radiative transfer solver:

- [**CanopyOptics.jl**](https://github.com/RemoteSensingTools/CanopyOptics.jl)
  provides leaf-angle distributions, PROSPECT leaf optics, and
  canopy-scattering kernels. It feeds vSmartMOM's `CanopySurface`.
- **OceanOptics.jl** (this package) provides pure-water, CDOM,
  phytoplankton, and non-algal-particle IOPs, phase functions, and
  inelastic sources. It feeds vSmartMOM's `OceanSurface`
  (to-be-implemented; see the plan in
  `vSmartMOM/docs/dev_notes/ocean.md` on the `unified-vsmartmom`
  branch).
- **vSmartMOM.jl** consumes both: atmospheric RT layers couple to
  either canopy or ocean as a lower boundary via the
  `AbstractSurfaceType` interface. The atmospheric solver, adding-
  doubling kernel, and polarization machinery are shared.

The design is deliberately parallel to CanopyOptics: Ocean and
Canopy are the two non-atmospheric physics the atmospheric solver
integrates. They share patterns (abstract material types,
spectral-grid conventions, data-artifact loading) where it makes
sense, and diverge where the physics differs (water scatters
everywhere; leaves don't; fluorescence/Raman emission from water
constituents has no canopy analog).

A longer-term goal is that design lessons from OceanOptics.jl —
particularly the orthogonal-axis type system described in §4 —
can be ported back into CanopyOptics.jl.

### 2.1 External reference codes

| Code | License | Language | Status |
|---|---|---|---|
| [HydroLight](https://www.numericaloptics.com/hydrolight.html) | Commercial (Numerical Optics Ltd.) | Fortran | Industry standard; not source-available |
| [OSOAA](https://github.com/CNES/RadiativeTransferCode-OSOAA) | GPLv3 | Fortran + C | Open; vector coupled atm-ocean |
| [HYDROPT](https://www.mdpi.com/2072-4292/13/15/3006) | Open | Python | Wraps HydroLight LUT results; not a solver |

HydroLight is the *de facto* validation standard in the field but is
commercial. We cannot compare OceanOptics.jl against its source, only
against published HydroLight outputs (including Mobley et al. 2002,
ESSD 2023 synthetic database, IOCCG reports). OSOAA is the open-source
target for numerical validation of the coupled atmosphere-ocean result;
it's the closest free analog to what vSmartMOM + OceanOptics will be.

### 2.2 Authoritative knowledge sources

- **[Ocean Optics Web Book](https://www.oceanopticsbook.info/)**
  (Mobley, continuously updated). The go-to online reference for
  everything from the definition of IOPs to modern bio-optical
  models. Chapters we return to throughout this memo: "Optical
  Constituents of the Ocean", "Scattering", "Radiative Transfer
  Theory".
- **Mobley (1994)** *Light and Water: Radiative Transfer in Natural
  Waters*, Academic Press. The canonical textbook. Out of print but
  widely available as PDF; Mobley has effectively released it online
  in the form of the Web Book.
- **Fell (1997)** PhD thesis (German), provided as
  `vSmartMOM/docs/papers/1997_Thesis_Fell.pdf`. Applies the Matrix
  Operator Method to the coupled atmosphere-ocean system. Primary
  reference for the fluorescence source operator (§2.2.9) and the
  Fresnel interface (§3).
- **[IOCCG Reports](https://ioccg.org/what-we-do/ioccg-reports/)**,
  particularly Report 5 (remote sensing of IOPs) and Report 14
  (phase function effects). Open-access consensus documents from
  the International Ocean Colour Coordinating Group.

---

## 3. Design principles

Four principles constrain the design. They are enumerated here so
that later decisions can be traced back to them.

### 3.1 RT-solver neutrality

The package's primary output is an `OceanLayerOptics` bundle: per-λ
optical depth τ, single-scattering albedo ϖ, and Legendre-expansion
coefficients β_ℓ of the mixed phase function. This is the smallest
stable interface the package can produce that every plane-parallel
ocean RT solver (vSmartMOM, HydroLight, OSOAA, Monte-Carlo codes)
can consume. The package itself contains no RT solver.

Translation from `OceanLayerOptics` into a specific solver's internal
representation — vSmartMOM's `CoreScatteringOpticalProperties` with
`Z`-matrices, HydroLight's tabular I/O, a Monte-Carlo code's
per-photon sampler — happens on the *solver* side. This matches the
CanopyOptics.jl / vSmartMOM precedent: the data package has zero
knowledge of which consumer imports it, and `Project.toml` lists no
solver packages under `[deps]` or `[weakdeps]`.

### 3.2 Multiple dispatch on three orthogonal axes

A biogeochemical constituent is described by three independent
attributes. The type system represents each as a separate axis so
that adding a new parameterization doesn't require modifying other
constituents.

| Axis | Examples | Implementation |
|---|---|---|
| **Optical role** | absorber, scatterer, absorbing-scatterer | abstract supertype |
| **Biogeochemical identity** | water, CDOM, phytoplankton, NAP | concrete struct |
| **Parameterization tag** | Pope-Fry 1997, Bricaud 1995, Mobley-Case-1 | singleton type parameter |

Concretely, `Phytoplankton{FT, Bricaud1995}` is a different type
from `Phytoplankton{FT, Gordon1992}`. Swapping parameterization is a
single-symbol change; dispatch resolves at compile time; no runtime
branches. This pattern was inspired by Oceananigans.jl's
`AbstractArchitecture` (CPU/GPU as type parameters driving compile-
time dispatch) and by CliMA's physics-closure hierarchies.

CanopyOptics.jl currently uses a flatter hierarchy (PROSPECT is one
function with branches on model version). A future CanopyOptics
refactor could adopt the three-axis pattern described here and gain
the same compile-time dispatch benefits — this is noted for the
cross-package learning loop mentioned in §2.

### 3.3 Physical quantities as values with algebra

The package carries three structurally distinct physical objects:

- `IOP{FT}`: a bundle `(a, b, bb)` of absorption, scattering, and
  backscattering coefficients. Overloads `+` for constituent mixing
  using the extensive-quantity sum (each coefficient adds).
- `OceanLayerOptics{FT}`: the ready-for-solver bundle
  `(τ, ϖ, β_ℓ, …)`. Overloads `*` for vertical concatenation of
  layers.
- `AbstractOceanPhaseFunction{FT}`: a phase-function object (e.g.
  `FournierForandPhase`, `RayleighWaterPhase`). Knows how to
  evaluate `p(cos θ)` and how to produce its Legendre expansion.

The algebra is chosen to match vSmartMOM's `CoreScatteringOpticalProperties`
conventions (which also overload `+` for mixing and `*` for vertical
concatenation), so the adapter is a near-identity operation.

### 3.4 Data provenance as a first-class concern

Every tabulated spectral quantity in the package is traceable to a
specific, citable, redistributable source. Data files live in a
`data/` directory (for the first release) or under Julia's
`Artifacts.toml` system (for subsequent releases, matching
CanopyOptics.jl's PROSPECT data pattern). Each file carries its own
header naming the original publication, data range, and units. Unit
conversions happen at load time, not at evaluation time.

The rationale: ocean-optics users frequently want to know which
specific parameterization a number came from (Pope-Fry 1997 vs.
Mason-Cone-Fry 2016 matters in the UV). Hiding the provenance inside
a generic "pure water absorption" function is a disservice.

---

## 4. Type system

### 4.1 The role hierarchy

```julia
abstract type AbstractOceanConstituent{FT<:AbstractFloat} end

abstract type AbstractAbsorber{FT}            <: AbstractOceanConstituent{FT} end
abstract type AbstractScatterer{FT}           <: AbstractOceanConstituent{FT} end
abstract type AbstractAbsorbingScatterer{FT}  <: AbstractOceanConstituent{FT} end
```

The three concrete branches correspond to physically meaningful
categories. Each has its own required method contract:

| Branch | Required methods | Zero fallbacks |
|---|---|---|
| `AbstractAbsorber` | `absorption(c, λ)` | `scattering(c, λ) = 0`, `backscattering(c, λ) = 0` |
| `AbstractScatterer` | `scattering(c, λ)`, `phase_function(c)` | `absorption(c, λ) = 0` |
| `AbstractAbsorbingScatterer` | all four | none |

The `attenuation(c, λ) = a + b` default is provided. Subtypes
typically don't override it; the zero fallbacks keep the IOP
summation branch-free.

### 4.2 Parameterization tags

Each identity — water, CDOM, phytoplankton, non-algal particles —
has its own small hierarchy of parameterization-model tags:

```julia
abstract type AbstractPureWaterAbsorptionModel end
struct SmithBaker1981      <: AbstractPureWaterAbsorptionModel end
struct PopeFry1997         <: AbstractPureWaterAbsorptionModel end
struct MasonConeFry2016    <: AbstractPureWaterAbsorptionModel end
struct CombinedPopeFryMason <: AbstractPureWaterAbsorptionModel end

abstract type AbstractPureWaterScatteringModel end
struct Morel1974  <: AbstractPureWaterScatteringModel end
struct ZhangHu2009 <: AbstractPureWaterScatteringModel end
```

Pure water carries **two** independent model tags because absorption
and scattering come from separate literatures that advance at
different rates. HydroLight's default is `PopeFry1997 + ZhangHu2009`
and we match that.

The concrete struct:

```julia
struct PureWater{FT, A<:AbstractPureWaterAbsorptionModel,
                     S<:AbstractPureWaterScatteringModel} <: AbstractAbsorbingScatterer{FT}
    absorption_model::A
    scattering_model::S
    temperature::FT       # [°C]
    salinity::FT          # [PSU]
    depolarization::FT    # for the Rayleigh-water phase function
end
```

Dispatched methods:

```julia
absorption(w::PureWater{FT, PopeFry1997,  S},    λ) where {FT, S}    = ...
absorption(w::PureWater{FT, SmithBaker1981, S},  λ) where {FT, S}    = ...
scattering(w::PureWater{FT, A, ZhangHu2009},     λ) where {FT, A}    = ...
scattering(w::PureWater{FT, A, Morel1974},       λ) where {FT, A}    = ...
```

Zero runtime dispatch cost. Adding a new parameterization is a new
tag + one method; nothing else in the package changes.

### 4.3 Phase functions

```julia
abstract type AbstractOceanPhaseFunction{FT<:AbstractFloat} end

struct RayleighWaterPhase{FT}   <: AbstractOceanPhaseFunction{FT}
    depolarization::FT
end

struct FournierForandPhase{FT}  <: AbstractOceanPhaseFunction{FT}
    n::FT          # real refractive index
    μ_j::FT        # Junge slope
    renorm::FT     # normalization factor (see §6.2)
end

struct PetzoldPhase{FT}         <: AbstractOceanPhaseFunction{FT} end
struct HenyeyGreensteinPhase{FT} <: AbstractOceanPhaseFunction{FT}
    g::FT          # asymmetry parameter
end
```

Contract: a concrete subtype implements `phase_function_value(pf, μ)`
*or* `phase_function_moments(pf, ℓ_max)`. The other is derived
numerically via Gauss-Legendre quadrature. Closed-form moments
(Rayleigh, Henyey-Greenstein) override for speed.

This is the same trait-style pattern used in Oceananigans'
`AbstractAdvectionScheme`: a concrete scheme supplies either the
explicit flux function or the stencil; other methods are derived.

### 4.4 Inelastic sources (Phase 2)

```julia
abstract type AbstractOceanInelasticProcess{FT} end

struct IsotropicFluorescence{FT, A<:AbstractAbsorber{FT},
                                 E<:AbstractSpectralRedistribution{FT}} <: AbstractOceanInelasticProcess{FT}
    absorber::A           # e.g. Phytoplankton — determines a_φ(λ')
    emission::E           # e.g. GaussianEmission(685 nm, 10.6 nm)
    quantum_yield::FT     # φ; typical 0.003 for Chl-a
    excitation_range::Tuple{FT,FT}    # e.g. (400.0, 700.0)
end

struct WaterRaman{FT} <: AbstractOceanInelasticProcess{FT}
    shift_cm_inv::FT      # ~3400 cm⁻¹ for liquid water
    bandwidth_cm_inv::FT  # ~200 cm⁻¹ Gaussian half-width
    depolarization::FT    # Cabannes depolarization
end

struct CDOMFluorescence{FT, E<:AbstractSpectralRedistribution{FT}} <: AbstractOceanInelasticProcess{FT}
    cdom::CDOM{FT}
    emission::E
    quantum_yield::FT
end
```

An `AbstractOceanInelasticProcess` *describes* the physics. A
separate adapter in `ext/OceanOpticsVSmartMOMExt` constructs a
vSmartMOM `AbstractRamanType` subtype (`OceanRS`) from a list of
processes — this is where the physics description meets the RT
solver's expected data layout. See §8 for details.

### 4.5 Layers and columns

```julia
struct OceanLayer{FT}
    depth_top::FT
    depth_bottom::FT
    constituents::Vector{<:AbstractOceanConstituent{FT}}
    fluorophores::Vector{<:AbstractOceanInelasticProcess{FT}}
end

struct OceanColumn{FT}
    layers::Vector{OceanLayer{FT}}   # surface-first, contiguity validated
end
```

The column type deferred to Phase 2 — Phase 1 just needs a single
`OceanLayer` to test end-to-end. Column construction includes a
contiguity check (no gaps, no overlaps, starts at z=0).

Convenience constructors following Fell (1997) §5.2.1:

```julia
fell_column(constituents; depth=60.0) -> OceanColumn
# → 1m layers 0-10m, 2m layers 10-20m, 5m layers 20-60m
uniform_column(constituents; n_layers=10, depth=50.0) -> OceanColumn
```

---

## 5. Pipeline: from biogeochemistry to RT-solver input

```
     biogeochemical state
           │
           │  PureWater{FT, PopeFry1997, ZhangHu2009}(T=20, S=35)
           │  CDOM{FT}(a_ref=0.1, λ_ref=440, slope=0.014)
           │  Phytoplankton{FT, Bricaud1995}(Chl=0.5)
           ▼
     OceanLayer{FT}
           │
           │  layer_optics(layer, λ_grid; ℓ_max=64)
           │     for each λ:
           │        a = Σ absorption(c, λ)          for c in constituents
           │        b = Σ scattering(c, λ)          for c in scatterers
           │        bb = Σ backscattering(c, λ)     for c in scatterers
           │        β_ℓ = mix_phase_moments(
           │                 [scattering(c, λ) for c in scatterers],
           │                 [phase_function_moments(phase_function(c), ℓ_max)
           │                  for c in scatterers])
           │        τ(λ) = (a + b) * Δz
           │        ϖ(λ) = b / (a + b)
           ▼
     OceanLayerOptics{FT}
           │
           │  (solver-side, not in this repo)
           │  vSmartMOM builds CoreScatteringOpticalProperties from β:
           │     Z⁺⁺, Z⁻⁺ = vSmartMOM.Scattering.compute_Z_matrices(
           │                   β, quad_points, pol_type)
           ▼
     vSmartMOM.CoreRT.CoreScatteringOpticalProperties
```

Every arrow inside this repo is a function call; every box is a
value. No hidden state, no mutable workspaces, no implicit
initialization. The full forward path is side-effect-free through
`OceanLayerOptics`; the Z-matrix and solver-specific packing steps
happen on the consumer side (see §3.1).

---

## 6. Phase functions in depth

### 6.1 Why Fournier-Forand, not Petzold

Mobley et al. (2002) Appl. Opt. 41, 1035, established
Fournier-Forand (Fournier & Jonasz 1999) as the modern default for
ocean-particle phase functions. Its advantages over the tabulated
Petzold "San Diego Harbor" phase function:

- Closed analytic form: `p(cos θ) = f(n, μ_j; cos θ)` where `n` is
  the real refractive index and `μ_j` is the Junge slope.
- Parameterizable by a single scalar: backscatter fraction
  `B = b_b / b`, via the Mobley (2002) inversion.
- Matches measured phase functions across three orders of magnitude
  of forward-scattered intensity (better than Petzold).

We ship `FournierForandPhase` as the default particle phase
function. `PetzoldPhase` (tabulated, bundled as `data/petzold_sdh.csv`
from Mobley 1994 Appendix) is available for reproducing legacy
results but not recommended.

### 6.2 Renormalization

The Fournier & Jonasz (1999) closed form uses an approximate
back-hemisphere correction that breaks exact normalization by
~1% for typical `(n, μ_j)`. We compute a renormalization factor at
construction time using `QuadGK.quadgk` with forced subdivision near
`μ = 1` (where the forward peak is integrable but near-singular,
`p ~ 1/sin²(θ/2)`), and store it as a struct field. Every
`phase_function_value` call multiplies by `renorm`, making
`(1/2) ∫₋₁¹ p(μ) dμ ≡ 1` exactly.

This is a Julia-idiomatic solution: the normalization is baked into
the immutable struct, so there is no way to accidentally use an
un-normalized phase function.

### 6.3 Mixing

Two or more scatterers in a layer mix with a scattering-weighted
average of phase-function Legendre moments:

```
β_ℓ^mix = (b_1 β_ℓ^{(1)} + b_2 β_ℓ^{(2)} + …) / (b_1 + b_2 + …)
```

This is exact at the Legendre-moment level (where mixing is linear).
Doing it on `p(μ)` samples and then re-expanding would lose
precision and require a quadrature pass per mixture. The
`mix_phase_moments` utility in `src/phase/moments.jl` implements
this.

---

## 7. Inelastic sources and the vSmartMOM Raman framework

### 7.1 The Fell factorization

Fell (1997) §2.2.9 gives the clean mathematical structure for any
inelastic source in the Matrix Operator Method. Translated from the
thesis (German) equations 2.86, 2.88, 2.91, 2.94:

```
J^I(τ, μ, φ, λ)  =  (1/c(τ)) · ∫_{λ'} b^I(τ; λ', λ) · ⟨L(τ; λ')⟩_Ω · p^I(cos θ; λ', λ) dλ'
                                     ↑                ↑              ↑
                       coupling    scalar irradiance   angular pattern
                       strength    from pass-1 solve

     b^I(τ; λ', λ) = a^I(τ; λ') · f^I(λ', λ)
```

where `a^I(λ')` is the absorption coefficient of the fluorescing /
Raman-scattering material at the excitation wavelength λ', and
`f^I(λ', λ)` is the spectral redistribution function in energy units
(normalized `∫ f dλ = 1`).

For isotropic emission (chlorophyll fluorescence per Mobley 1994
p.308, CDOM fluorescence in solution, and approximately water Raman):

- The angular integral reduces to scalar irradiance `E°(τ, λ')`.
- The source contributes only to Fourier mode `m = 0`.

This is one unified structure; chlorophyll fluorescence, phycocyanin
fluorescence, CDOM fluorescence, and water Raman are all instances
of it with different `a^I(λ')` and `f^I(λ', λ)`.

### 7.2 Mapping to vSmartMOM's AbstractRamanType

vSmartMOM's `src/Inelastic/types.jl` on the **sanghavi branch**
already defines the data contract. Every `AbstractRamanType` subtype
(`RRS`, `VS_0to1`, `sol_RRS`, etc.) carries the same fields:

```julia
ϖ_λ₁λ₀::Vector{FT}      # inelastic coupling strengths
i_λ₁λ₀::Vector{Int}     # excitation-to-emission index map
Z⁻⁺_λ₁λ₀::Matrix{FT}    # backward phase-matrix block at couplings
Z⁺⁺_λ₁λ₀::Matrix{FT}    # forward phase-matrix block
greek_raman::GreekCoefs # Greek expansion
F₀::Matrix{FT}          # source spectrum (solar / stellar)
SIF₀::Matrix{FT}        # Solar-Induced Fluorescence slot — already present!
```

The presence of `SIF₀` in every existing subtype is deliberate — the
framework was built anticipating fluorescence as an extension. This
is a "don't build parallel machinery" signal: ocean fluorescence and
ocean Raman should plug in as a new subtype `OceanRS <: AbstractRamanType`,
populated by a function patterned on
`getRamanSSProp!(RS_type, depol, λ, grid_in)`.

### 7.3 The adapter pattern

```julia
# In OceanOptics.jl (data-producer, RT-solver-neutral)
struct IsotropicFluorescence{FT, A, E} <: AbstractOceanInelasticProcess{FT}
    absorber::A
    emission::E
    quantum_yield::FT
    excitation_range::Tuple{FT, FT}
end

# In vSmartMOM.jl (consumer)
function getOceanInelasticProp!(RS_type::OceanRS{FT},
                                processes::Vector{<:AbstractOceanInelasticProcess},
                                layer::OceanLayer,
                                λ_grid::AbstractVector{FT}) where {FT}
    # Fold each process into ϖ_λ₁λ₀, i_λ₁λ₀, Z_λ₁λ₀ fields of RS_type.
    # Fluorescence: isotropic, δ₁(0, m) → contributes only to m=0.
    # Water Raman: Cabannes phase function with depolarization δ.
    # CDOM fluorescence: isotropic, broader emission band.
    # ...
end
```

The existing `elemental_inelastic!`, `doubling_inelastic!`,
`interaction_inelastic!` kernels in
`vSmartMOM/src/CoreRT/CoreKernel/` consume the populated `OceanRS`
struct unchanged — zero new kernel code.

Note the direction: `getOceanInelasticProp!` lives in vSmartMOM,
not here, so that OceanOptics.jl declares no solver dependency
(see §3.1). The types `IsotropicFluorescence`, `WaterRaman`,
`CDOMFluorescence` and their `excitation_absorption` /
`emission` / `excitation_range` / `is_isotropic` kernel trait
functions are the stable public API vSmartMOM consumes.

### 7.4 Two-pass solve

Both fluorescence and Raman require a two-pass RT solve. Fell (1997)
§2.2.9 explains: Pass 1 solves the elastic RT problem to obtain the
scalar irradiance `E°(τ, λ')`. Pass 2 uses `E°` to evaluate the
activation function `A(τ) = ∫ a^I(λ') E°(τ, λ') λ' dλ'` and add the
source term. vSmartMOM's existing Raman pipeline already implements
two-pass for atmospheric rotational Raman; the ocean inherits this
machinery.

---

## 8. Data sources

All tables below are redistributable and available as plain-text
files. For the first release, download to `data/`; for later
releases, use `Artifacts.toml` (matching CanopyOptics.jl's PROSPECT
pattern).

### 8.1 Pure-water absorption

| Source | Range (nm) | Res. | URL | Units |
|---|---|---|---|---|
| Pope & Fry (1997) | 380–727.5 | 2.5 nm | <https://omlc.org/spectra/water/data/pope97.txt> | 1/cm |
| Smith & Baker (1981) | 200–800 | 10 nm | <https://omlc.org/spectra/water/data/smith81.txt> | 1/cm |
| Segelstein (1981) | 10 nm – 10 m | dense | <https://omlc.org/spectra/water/data/segelstein81.txt> | n, k |
| Hale & Querry (1973) | 200 nm – 200 µm | variable | <https://omlc.org/spectra/water/data/hale73.txt> | n, k |

**OMLC units are `1/cm`** — multiply by 100 to get `1/m`.

**Recommended default**: the piecewise `CombinedPopeFryMason` model —
Mason-Cone-Fry (2016) below 380 nm (UV), Pope-Fry in 380–700 nm, and
Smith-Baker above 700 nm (NIR). Matches HydroLight's default "aw"
table.

Mason-Cone-Fry (2016) data is not on OMLC; Table 2 of the paper
(Appl. Opt. 55, 7163) lists values at 250–550 nm which we transcribe
into `data/mason_cone_fry_2016.csv`.

Temperature corrections (Pegau, Gray & Zaneveld 1997, Appl. Opt. 36,
6035) apply above ~600 nm. Values in Table 2 of the paper;
transcribed directly.

### 8.2 Pure-seawater scattering — Zhang et al. (2009)

- Paper: Zhang, Hu & He, Opt. Express 17, 5698 (2009).
- **MATLAB reference implementation**:
  <https://www.seanoe.org/data/00318/42916/> (direct download).
- ~100 lines, straightforward Julia port.

The older Morel (1974) power law `b_w(λ) = 0.00193 · (550/λ)^4.32`
is kept as a fallback (tag: `Morel1974`) but Zhang-Hu is the
default. Linear extrapolation of Morel's pure-water result to
seawater overestimates `b` by up to 30%.

### 8.3 Phytoplankton absorption — Bricaud et al. (1995, 1998)

Empirical power-law form:

```
a_φ*(λ) = A(λ) · Chl^(-E(λ))    [m²/mg]
```

where `A(λ)` and `E(λ)` are tabulated spectra at 1 nm resolution,
400–700 nm, from fits to 815 globally-distributed spectra. Table 2
of:

- Bricaud, Babin, Morel & Claustre (1995), JGR 100, 13321,
  doi:10.1029/95JC00463
- Bricaud, Morel, Babin, Allali & Claustre (1998), JGR 103, 31033,
  doi:10.1029/98JC02712

No public data file. The table is transcribed into
`data/bricaud_1995.csv` (~300 lines, one-time transcription).

Alternative parameterizations for future implementation:

- Gordon (1992): simpler, less accurate, but historical baseline.
  `tag: Gordon1992`.
- Ciotti, Lewis & Cullen (2002): accounts for algal size structure.
  `tag: Ciotti2002`.
- Mobley "New Case 1 IOP Model": HydroLight 5 default combining
  several empirical pieces. See
  <https://www.oceanopticsbook.info/view/optical-constituents-of-the-ocean/level-2/new-iop-model-case-1-water>.

### 8.4 CDOM absorption

Bricaud, Morel & Prieur (1981), Limnol. Oceanogr. 26, 43:

```
a_g(λ) = a_g(λ_ref) · exp(-S · (λ - λ_ref))     [1/m]
```

Typical slope `S = 0.014 nm⁻¹`, `λ_ref = 440 nm`. Babin et al.
(2003), JGR 108, 3211 gives coastal slope range 0.011–0.020 nm⁻¹.

No data file needed — the model is fully analytic.

### 8.5 Non-algal particles

- Babin et al. (2003), JGR 108, 3211: specific absorption with
  exponential slope, similar to CDOM form.
- Morel & Maritorena (2001), JGR 106, 7163: scattering
  parameterization for Case-1 waters.
- For inorganic/mineral particles specifically, Fell (1997) §5.1.4
  gives a crude `c(λ) = c_ref · (440/λ)` extinction law, which we
  use as the simplest `NonAlgalParticles` model.

### 8.6 Fournier-Forand phase function

No data file — closed-form analytic.

- Fournier & Forand (1994), Proc. SPIE 2258, 194.
- Fournier & Jonasz (1999), Proc. SPIE 3761, 62. Latest form.
- Mobley, Sundman & Boss (2002), Appl. Opt. 41, 1035, Eqs. 10, 17,
  18 for B → (n, μ_j) inversion.
- Ocean Optics Web Book:
  <https://www.oceanopticsbook.info/view/scattering/the-fournier-forand-phase-function>.

### 8.7 Petzold "San Diego Harbor" phase function

Tabulated VSF from Petzold (1972), SIO Ref. 72-78. Mobley (1994)
*Light and Water* Appendix A.4 transcribes the table. Bundle as
`data/petzold_sdh.csv`.

### 8.8 Raman scattering

- **Haltrin & Kattawar (1991, 1993)**, Appl. Opt. 30, 630 and 32,
  5356: canonical references for water Raman in ocean RT. Gives
  Raman scattering coefficient `b_w^R(λ)` and the Stokes-shifted
  emission band.
- **Bartlett et al. (1998)**, JGR 103, 17875: more recent Raman
  measurements with temperature and salinity dependence.

Water Raman has a fixed wavenumber shift (~3400 cm⁻¹) and a
Gaussian-like emission band of width ~200 cm⁻¹. In terms of the §7
factorization: `a^R(λ') = b_w^R(λ')` (the Raman "absorption" is
actually a scattering cross-section, but in the inelastic-source
formalism it plays the absorption role) and `f^R(λ', λ)` is a
Gaussian centered at wavelength corresponding to λ' shifted by
3400 cm⁻¹.

Ocean Optics Web Book Raman page:
<https://www.oceanopticsbook.info/view/scattering/level-2/raman-scattering>.

### 8.9 Chlorophyll fluorescence

- **Fell (1997)** §2.2.9: Gaussian emission at 685 nm with σ = 10.6 nm
  (FWHM ≈ 25 nm), quantum yield φ_C = 0.003 (North Sea average).
- **Gordon (1979)**, Appl. Opt. 18, 1161: original Gaussian model.
- **Fischer & Kronfeld (1990)**, Deep-Sea Res. 37, 865: the study
  Fell cites for the 0.003 quantum yield.
- No data file — analytic model.

---

## 9. File layout

```
OceanOptics.jl/
├── Project.toml
├── Artifacts.toml              # later: auto-download data tables
├── DESIGN.md                   # this file
├── README.md                   # quick-start and package overview
│
├── src/
│   ├── OceanOptics.jl          # module entry point
│   │
│   ├── types/
│   │   ├── abstract_types.jl   # AbstractOceanConstituent + role subtypes
│   │   ├── iop.jl              # IOP{FT} with + algebra
│   │   ├── layer_optics.jl     # OceanLayerOptics (the RT-solver contract)
│   │   ├── ocean_layer.jl      # OceanLayer
│   │   └── ocean_column.jl     # OceanColumn + fell_column, uniform_column
│   │
│   ├── materials/
│   │   ├── pure_water.jl       # PureWater{FT, A, S} + all model dispatch
│   │   ├── cdom.jl             # CDOM (pure absorber)
│   │   ├── phytoplankton.jl    # Phytoplankton{FT, M} — Bricaud / Gordon / ...
│   │   ├── nap.jl              # NonAlgalParticles (detritus + minerals)
│   │   └── refractive_index.jl # seawater n(λ, T, S)
│   │
│   ├── phase/
│   │   ├── abstract_phase.jl   # contract + default fallbacks
│   │   ├── moments.jl          # numerical Legendre expansion + mixing
│   │   ├── rayleigh_water.jl   # closed-form depolarized Rayleigh
│   │   ├── fournier_forand.jl  # FJ1999 + Mobley 2002 inversion
│   │   ├── petzold.jl          # tabulated (loads data/petzold_sdh.csv)
│   │   └── henyey_greenstein.jl # analytical, legacy
│   │
│   ├── inelastic/              # Phase 2
│   │   ├── processes.jl        # AbstractOceanInelasticProcess hierarchy
│   │   ├── redistribution.jl   # GaussianEmission, RamanShift, Tabulated
│   │   ├── fluorescence.jl     # concrete fluorescence process types
│   │   └── raman.jl            # WaterRaman
│   │
│   ├── mixing/
│   │   └── layer_assembly.jl   # layer_optics(layer, λ, ℓ_max)
│   │
│   └── io/
│       └── data_loaders.jl     # parse OMLC txt, Bricaud csv, etc.
│
├── ext/                              # optional, non-solver extensions only
│   └── (future: OceanOpticsUnitfulExt/unitful_wrappers.jl
│        for u"nm" wavelength grids)
│
├── data/
│   ├── pope_fry_1997.csv
│   ├── smith_baker_1981.csv
│   ├── mason_cone_fry_2016.csv
│   ├── pegau_1997_temperature.csv
│   ├── bricaud_1995.csv
│   └── petzold_sdh.csv
│
└── test/
    ├── runtests.jl
    ├── test_pure_water.jl
    ├── test_fournier_forand.jl
    ├── test_layer_assembly.jl
    └── test_vsmartmom_integration.jl   # runs only if vSmartMOM loadable
```

The `ext/` directory is reserved for non-solver convenience
extensions (Unitful wavelength-grid support is the obvious
candidate). RT-solver integrations — vSmartMOM's `OceanSurface`,
`OceanRS`, `CoreScatteringOpticalProperties` adapter — live in
the *solver* package, not here, per the CanopyOptics.jl pattern
discussed in §3.1.

---

## 10. Implementation roadmap

### Phase 1 — elastic path, single layer

Target: `layer_optics(layer, λ_grid)` works for pure seawater, CDOM,
and chlorophyll. Result validates against published HydroLight
IOP tables.

| Deliverable | Notes |
|---|---|
| Type spine: `AbstractOceanConstituent` + roles | §4.1 |
| `IOP{FT}` + mixing algebra | §3.3 |
| `PureWater` with 4 × 2 model tags | §4.2 |
| `CDOM` analytic | §8.4 |
| `Phytoplankton{FT, Bricaud1995}` | §8.3 |
| `FournierForandPhase` with renorm | §6.2 |
| `RayleighWaterPhase` closed-form | §4.3 |
| `OceanLayer` + `layer_optics` | §5 |
| Data loaders for OMLC + Bricaud | §8 |
| Tests: absorption at 440 nm, normalization, layer assembly | — |

### Phase 2 — inelastic framework

Target: fluorescence and water Raman as `AbstractOceanInelasticProcess`
instances, buildable into an `OceanRS <: AbstractRamanType` in the
vSmartMOM extension.

| Deliverable | Notes |
|---|---|
| `AbstractOceanInelasticProcess` hierarchy | §4.4 |
| `GaussianEmission`, `WavenumberShiftEmission` | §7.1 |
| `IsotropicFluorescence`, `WaterRaman`, `CDOMFluorescence` | §7.1 |
| `OceanRS` in vSmartMOM extension | §7.2 |
| `getOceanInelasticProp!` | §7.3 |
| Validation: chlorophyll fluorescence magnitude (Fell Fig 5.5) | — |

### Phase 3 — coupled RT in vSmartMOM

`OceanSurface <: AbstractSurfaceType` with adding-doubling through
ocean sub-layers, Fresnel interface coupling, plus the
`OceanRS <: AbstractRamanType` adapter, plus the
`CoreScatteringOpticalProperties` builder from an `OceanLayerOptics`.
**All of this lives in vSmartMOM.jl, not in OceanOptics.jl** (see
§3.1 for the rationale). From this repo's perspective Phase 3 is a
no-op — the contract to honor is keeping `OceanLayerOptics`, the
`AbstractOceanInelasticProcess` kernel traits, and the
`AbstractOceanConstituent` hierarchy stable so the solver side can
code against them without churn. See
`vSmartMOM/docs/dev_notes/ocean.md` for the detailed solver-side
plan.

### Phase 4 — polish

- Full polarized Greek-coefficient expansions (not just β)
- `OceanColumn` with Fell-grid helper
- Artifact-based data loading
- JOSS paper

---

## 11. Validation strategy

Three tiers:

1. **Unit tests** against published reference values:
   - Pope-Fry 440 nm: `a_w = 0.00635 m⁻¹`
   - Morel 1974 scattering at 500 nm: `b_w = 0.00317 m⁻¹`
   - Zhang-Hu at 500 nm, S=35, T=20°C: `b_w ≈ 0.00194 m⁻¹`
   - Fournier-Forand normalization: `(1/2) ∫₋₁¹ p dμ = 1 ± 10⁻⁶`
   - Bricaud 1995 at Chl=1 mg/m³, 440 nm: `a_φ* ≈ 0.03 m⁻¹`

2. **Cross-model closure**: `PureWater{PopeFry1997} + CDOM + Phytoplankton{Bricaud1995}`
   for Case-1 water at Chl=1 mg/m³, compared against Mobley "New
   Case 1 IOP Model" values from the Ocean Optics Web Book.

3. **End-to-end RT closure** (Phase 3): run vSmartMOM +
   `OceanSurface` on the IOCCG Report 5 synthetic dataset and the
   ESSD 2023 database, compare to published HydroLight outputs. Also
   run against OSOAA for an open-source cross-check.

---

## 12. Open questions

These are design decisions deferred pending real use:

1. **Polarization**: scalar RT is the Phase-1 target. Full Greek
   coefficients (α, β, γ, δ, ε, ζ) for polarized Stokes-vector RT
   adds a ~6× expansion cost; add when vSmartMOM's `OceanSurface`
   actually uses polarization.

2. **Temperature-depth coupling**: T and S profiles affect Zhang-Hu
   scattering and Pegau temperature-corrected absorption. Phase 1
   treats them as per-layer constants; Phase 2 might need
   interpolation on a separate T(z), S(z) profile.

3. **Multi-species phytoplankton**: Bricaud 1995 is
   single-phytoplankton. Multi-species retrievals (e.g.
   diatom/cyanobacteria/prochlorococcus) need a phytoplankton
   composition type — deferred.

4. **APAR accumulation**: Absorbed Photosynthetically Active
   Radiation is a byproduct of the pass-1 solve. Adding an
   accumulator into the RT output requires vSmartMOM changes —
   coordinate with the solver side.

5. **AD compatibility**: the broader vSmartMOM refactor is moving to
   ForwardDiff-based Jacobians for retrieval. OceanOptics.jl
   functions should be AD-transparent. Early check: ensure no
   mutating operations on values derived from inputs, no
   branches on input values (use `ifelse` not `if`), no
   `Complex{FT}` paths for values that might be dualized.

---

## 13. References

**Package repositories**:

- OceanOptics.jl (this package): `github.com/RemoteSensingTools/OceanOptics.jl`
- CanopyOptics.jl: <https://github.com/RemoteSensingTools/CanopyOptics.jl>
- vSmartMOM.jl: <https://github.com/RemoteSensingTools/vSmartMOM.jl>
- Oceananigans.jl: <https://github.com/CliMA/Oceananigans.jl>
- OSOAA (GPLv3 reference): <https://github.com/CNES/RadiativeTransferCode-OSOAA>

**Authoritative references**:

- Ocean Optics Web Book: <https://www.oceanopticsbook.info/>
- Mobley (1994) *Light and Water*, Academic Press. Widely available
  as PDF through the author's affiliation at Sequoia Scientific.
- IOCCG Reports: <https://ioccg.org/what-we-do/ioccg-reports/>
- Fell (1997) PhD thesis, in `vSmartMOM/docs/papers/1997_Thesis_Fell.pdf`

**Key primary sources** (reference-page numbers are for Ocean Optics
Web Book "Optical Constituents of the Ocean" section):

- Pope & Fry (1997), Appl. Opt. 36, 8710 — pure water absorption,
  visible.
- Mason, Cone & Fry (2016), Appl. Opt. 55, 7163 — pure water
  absorption, UV.
- Smith & Baker (1981), Appl. Opt. 20, 177 — pure water absorption,
  200–800 nm.
- Zhang, Hu & He (2009), Opt. Express 17, 5698 — seawater
  scattering.
- Bricaud et al. (1995, 1998), JGR 100, 13321 / 103, 31033 —
  phytoplankton absorption.
- Bricaud, Morel & Prieur (1981), Limnol. Oceanogr. 26, 43 — CDOM
  absorption.
- Fournier & Jonasz (1999), Proc. SPIE 3761, 62 — phase function.
- Mobley, Sundman & Boss (2002), Appl. Opt. 41, 1035 — phase
  function parameterization.
- Haltrin & Kattawar (1991, 1993), Appl. Opt. 30, 630 / 32, 5356 —
  water Raman.
- Fell (1997), PhD thesis, FU Berlin — fluorescence source operator
  in MOM.
