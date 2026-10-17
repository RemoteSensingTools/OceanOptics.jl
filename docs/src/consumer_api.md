# Consumer API

OceanOptics.jl publishes an RT-agnostic interface. Downstream solvers
(vSmartMOM, HydroLight-style ODE codes, Monte-Carlo, SOS, …) consume it
*without* OceanOptics having any knowledge of the solver.
`Project.toml` therefore lists no RT-solver package under `[deps]` or
`[weakdeps]`; adapters live on the solver side.

This page is the consumer contract: what the numbers mean, in which
units and conventions, and which inputs are rejected rather than
approximated.

## Quantities and units

| Quantity | Where | Units | Notes |
|---|---|---|---|
| Wavelength `λ` | everywhere | nm | As tabulated by the source data; no air/vacuum conversion is applied. `OceanLayerOptics.λ` is strictly monotone, increasing *or* decreasing |
| Absorption `a`, scattering `b`, backscattering `b_b` | constituents, [`IOP`](@ref) | m⁻¹ | Volumetric coefficients, not spectral densities |
| Thickness `Δz` | [`OceanLayer`](@ref), [`OceanLayerOptics`](@ref) | m | Positive; depth increases downward from the surface |
| Optical depth `τ = (a + b) Δz` | `OceanLayerOptics.τ` | — | Per layer and wavelength |
| Single-scattering albedo `ϖ = b / (a + b)` | `OceanLayerOptics.ϖ` | — | In `[0, 1]` |
| Backscatter fraction `B = b_b / b` | `OceanLayerOptics.B` | — | Diagnostic only; in `[0, 1]` |
| Greek coefficients `α, β, γ, δ, ϵ, ζ` | `OceanLayerOptics` | — | `(ℓ_max + 1, nλ)` matrices, row `ℓ + 1`, column `iλ` |
| Fluorescence source `J_source` | `OceanLayerOptics.J_source` | W m⁻² sr⁻¹ nm⁻¹ | Optional; per-nm spectral density |
| Redistribution `emission(p, λ′, λ)` | inelastic processes | nm⁻¹ | Density per nm of emission wavelength |

Consumers working in wavenumber (vSmartMOM's core RT uses cm⁻¹) must
apply the spectral-density Jacobian to per-nm quantities, and permute
every per-wavelength array identically when reversing the order.

## Phase-function convention

* **Normalization:** `(1/2) ∫₋₁¹ p(μ) dμ = 1`, i.e. `∫ p dΩ / (4π) = 1`.
* **Expansion:** `p(μ) = Σ_ℓ β_ℓ P_ℓ(μ)` with
  `β_ℓ = (2ℓ+1)/2 ∫₋₁¹ p(μ) P_ℓ(μ) dμ`. The `(2ℓ+1)` factor is already
  included — adapters must not apply it again. Hence `β₀ = 1`, and
  `|β_ℓ| ≤ 2ℓ+1` for every non-negative phase function; both are checked
  by the `OceanLayerOptics` constructor. The bound catches a doubled
  `(2ℓ+1)` factor only for strongly peaked phase functions and never a
  missing one (e.g. DISORT-style `g^ℓ` moments), so compare `β₁/3` with
  the published asymmetry parameter when importing data.
* **Phase matrix:** the Greek coefficients define, per order `ℓ`,

  ```
  B_ℓ = [β  γ   0   0
         γ  α   0   0
         0  0   ζ   ϵ
         0  0  -ϵ   δ]
  ```

  (Sanghavi 2014 Eq. 16; de Haan et al. 1987), with the Hovenier sign of
  `γ` used by vSmartMOM: `γ₂ = +√(3/2)·dpl_p > 0` for Rayleigh
  scattering. Tables in other conventions (e.g. VLIDORT/Mishchenko-format
  expansion coefficients) differ at least in the sign of `γ`; convert them
  before use — vSmartMOM's `docs/src/pages/conventions.md` documents the
  known flips. The coefficients describe scattering relative to the
  scattering plane; rotation into the meridian-plane Stokes basis is the
  solver's job. `α, γ, ϵ, ζ` vanish at `ℓ < 2`, where their generalized
  spherical functions are identically zero; the constructors check this.
* **Truncation:** [`layer_optics`](@ref) returns the first `ℓ_max + 1`
  coefficients with no δ-M or other truncation correction; the solver
  applies its own. Numerically integrated moments (Fournier-Forand,
  Petzold) are rescaled so that `β₀ = 1`; this fixes the normalization
  but does not by itself establish the accuracy of the higher moments of
  strongly forward-peaked functions.

## Polarization capability

| Phase function | `is_polarizable` | Polarized coefficients |
|---|---|---|
| [`RayleighWaterPhase`](@ref) | yes | Closed form (Hansen & Travis 1974) |
| [`TabulatedGreekPhase`](@ref) | if all six are given | As tabulated |
| [`FournierForandPhase`](@ref), [`HenyeyGreensteinPhase`](@ref), [`PetzoldPhase`](@ref) | no | none — scalar only |

A scalar `β_ℓ` sequence does not determine the other five coefficients,
so OceanOptics never invents them:

* [`phase_matrix_moments`](@ref) throws for scalar-only phase functions.
* `layer_optics(...; polarized = true)` throws if any scatterer in the
  layer is scalar-only, unless the caller passes
  `scalar_phase_approximation = :depolarizing`. That option selects the
  completely depolarizing matrix `diag(p, 0, 0, 0)` for those scatterers
  (see [`depolarizing_phase_matrix_moments`](@ref)): a modelling choice,
  not real particle polarization.
* Every `OceanLayerOptics` is either scalar (`α, γ, δ, ϵ, ζ` all empty)
  or polarized (all five present); [`is_polarized`](@ref) tells which.

## Zero scattering

At a wavelength where a layer does not scatter, `ϖ = 0`, `B = 0`,
`β = [1, 0, 0, …]` and all polarized coefficients are zero. A layer
without any absorption or scattering has `τ = 0`.

## Tabulated inputs

Benchmarks prescribe optical properties directly. Three types carry such
values through the ordinary constituent and phase-function interfaces,
without any of the package's empirical bio-optical defaults:

* [`TabulatedAbsorber`](@ref)`(λ, a; total = false, source = "")`
* [`TabulatedScatterer`](@ref)`(λ, b, phase; source = "")` — non-absorbing;
  `b_b` is derived from the phase function, so the two cannot disagree.
* [`TabulatedGreekPhase`](@ref)`(β; α, γ, δ, ϵ, ζ, source, normalize)` or
  `TabulatedGreekPhase(λ, β; …)` for a wavelength-dependent table.

Values are defined only at their tabulated wavelengths; any other
wavelength throws instead of interpolating. Each object records its
`source`.

**Total versus constituent absorption.** Published absorption is often
the bulk value of the whole water body, already containing water and
CDOM — e.g. `a_blk = a_p + a_y + a_w` in Chowdhary et al. (2020)
Eq. (11b). Mark such data with `total = true`; `layer_optics` then
rejects any other absorbing constituent in the same layer, so water or
CDOM cannot be counted twice.

```julia
# Chowdhary et al. (2020) AOS-II/III ocean body, Table 2
λ     = [350.0, 450.0, 550.0, 650.0]
a_blk = TabulatedAbsorber(λ, [0.0204, 0.0092, 0.0565, 0.3400]; total = true,
                          source = "Chowdhary et al. 2020, Table 2")
water = TabulatedScatterer(λ, [0.0134, 0.0045, 0.0019, 0.0010],
                           RayleighWaterPhase(depolarization = 0.0);
                           source = "Chowdhary et al. 2020, Table 2")
body  = OceanLayer(0.0, 100.0, AbstractOceanConstituent{Float64}[a_blk, water])
opt   = layer_optics(body, λ; ℓ_max = 16, polarized = true)
opt.ϖ    # ≈ [0.39644970, 0.32846715, 0.03253425, 0.002932551] (Table 2)
```

Not yet supported: angular scattering-matrix tables such as
`F₁₁(Θ), F₁₂(Θ), F₃₃(Θ), F₃₄(Θ)`, which need a reviewed
generalized-spherical-function expansion before they can become Greek
coefficients, and automatic conversion between Greek-coefficient
conventions.

## Validity and provenance

Each empirical model documents its source and intended range on the
[Materials](@ref) page; bundled tables in `data/` carry provenance
headers. Interpolated tables are clamped to their end values outside the
tabulated range (`linterp`), so a query beyond a model's validity range
returns an edge value rather than an error — check the documented range
before using a model there. The `MasonConeFry2016` absorption table is
not bundled yet; its first absorption query throws.

## The stable contract

A downstream solver should consume only the symbols below. Anything
else in the package is internal and can change in any minor version.

* **Output bundle:** [`OceanLayerOptics`](@ref), [`layer_optics`](@ref),
  [`is_polarized`](@ref), [`n_moments`](@ref), [`n_wavelengths`](@ref).
* **Constituent hierarchy:** [`AbstractOceanConstituent`](@ref),
  [`AbstractAbsorber`](@ref), [`AbstractScatterer`](@ref),
  [`AbstractAbsorbingScatterer`](@ref), and the trait functions
  [`absorption`](@ref), [`scattering`](@ref), [`backscattering`](@ref),
  [`attenuation`](@ref), [`phase_function`](@ref),
  [`is_total_absorption`](@ref).
* **Phase functions:** [`AbstractOceanPhaseFunction`](@ref),
  [`phase_function_moments`](@ref), [`phase_matrix_moments`](@ref),
  [`is_polarizable`](@ref), [`is_wavelength_dependent`](@ref).
* **Tabulated inputs:** [`TabulatedAbsorber`](@ref),
  [`TabulatedScatterer`](@ref), [`TabulatedGreekPhase`](@ref).
* **Inelastic-process hierarchy:** [`AbstractOceanInelasticProcess`](@ref)
  and trait functions [`excitation_absorption`](@ref), [`emission`](@ref),
  [`excitation_range`](@ref), [`is_isotropic`](@ref),
  [`inelastic_coefficient`](@ref). Canonical docstrings in [Inelastic processes](@ref).
* **Layer / column:** [`OceanLayer`](@ref) and [`OceanColumn`](@ref).

## The two-pass solve convention (Fell 1997 §2.2.9)

Inelastic sources require two passes:

1. **Pass 1** — solve the elastic RT problem using `OceanLayerOptics` to
   get the internal radiance field at the excitation wavelengths `λ'`.

2. **Pass 2** — for each inelastic process `p` in
   `layer.fluorophores`, evaluate the source kernel

   ```math
   J^I(\tau, \lambda) = \int_{\lambda'} b^I(\lambda', \lambda) \,
                                       \langle L(\tau, \lambda') \rangle_\Omega \,
                                       p^I(\cos\theta; \lambda', \lambda) \, d\lambda'
   ```

   with `b^I(λ', λ) = excitation_absorption(p, λ') × emission(p, λ', λ)`.
   For isotropic processes (`is_isotropic(p) == true`) the angular factor
   `p^I` reduces to a scalar and the source contributes only to the
   `m = 0` Fourier mode; for isotropic emission the excitation enters
   through the scalar irradiance (total intensity integrated over solid
   angle, direct and diffuse, up and down), not the downwelling plane
   irradiance. Anisotropic processes such as [`WaterRaman`](@ref) need
   the angular excitation field.

The adapter code that packs these kernels into a solver's internal
format lives in the solver package.
