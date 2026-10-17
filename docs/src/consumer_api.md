# Consumer API

OceanOptics.jl publishes an RT-agnostic interface. Downstream solvers
(vSmartMOM, HydroLight-style codes, Monte Carlo, …) consume it *without*
OceanOptics having any knowledge of the solver: `Project.toml` lists no
RT package, and adapters live on the solver side.

This page is the contract: what the numbers mean, in which units and
conventions, and which inputs are rejected rather than approximated.

## What a solver receives

```julia
col = column_optics(column, λ; refractive_index, bottom, ℓ_max, polarized,
                    scalar_phase_approximation)
```

returns an [`OceanColumnOptics`](@ref): contiguous layers from the
surface down, the refractive index of the surface water relative to air
per wavelength, and the bottom ([`BlackBottom`](@ref) or
[`LambertianBottom`](@ref)). Each layer is an [`OceanLayerOptics`](@ref):

| Field / accessor | Meaning |
|---|---|
| `λ` | wavelengths `[nm]`, strictly increasing or decreasing |
| `depth_top`, `depth_bottom`, [`thickness`](@ref) | depths below the surface `[m]` |
| `τ_abs`, [`absorption_optical_depth`](@ref) | absorption optical depth of all absorbers |
| `components` | one [`ScatteringComponent`](@ref) per scattering constituent |
| [`optical_depth`](@ref), [`scattering_optical_depth`](@ref), [`single_scattering_albedo`](@ref) | derived totals |
| [`greek_expansion`](@ref)`(o, iλ)` | mixed expansion at a wavelength (scattering-weighted) |
| [`backscatter_fraction`](@ref)`(o)` | scattering-weighted backscatter fraction of the phase functions |

A [`ScatteringComponent`](@ref) carries a `label` (constituent type, or
the `source` of tabulated data), its scattering optical depth `τ` per
wavelength, its distinct expansions `greek` ([`GreekCoefficients`](@ref);
one unless the phase function depends on wavelength, in which case
`phase_index[iλ]` selects the entry at each wavelength), the backscatter
fraction of its phase function per entry, a `phase_id`, and a
`polarization` label:

| `polarization` | Meaning |
|---|---|
| `:scalar` | phase function only (`polarized = false`) |
| `:full` | the constituent's actual phase matrix |
| `:depolarizing` | its phase function with the completely depolarizing matrix `diag(p, 0, 0, 0)`, selected explicitly |

Keeping components separate lets a solver build each component's phase
operators once per band and truncate them individually (as it does for
aerosol types), and keeps every phase matrix's provenance visible.
Components with equal `phase_id` — in any layer of a column — have
identical expansions (`column_optics` shares the objects), so their phase
operators can be reused. Layers are ordered from the surface down;
[`interface_depths`](@ref)`(col)` gives the layer boundaries.

**Derivatives.** Every optical depth is linear in its constituent's
concentration (`τ_c = b_c Δz`, `τ_abs = Σ a_c Δz`) and the layer mixture is
linear in the component optical depths, so derivatives with respect to
component amounts follow directly from the component split; the
assembly itself is evaluated in floating point (`FT`).

**Deep water.** Layers have finite thickness. An optically deep ocean is
represented by a bottom layer thick enough for the bottom not to matter
(checked by convergence in depth), with `BlackBottom()`.

## Quantities and units

| Quantity | Units | Notes |
|---|---|---|
| Wavelength `λ` | nm | As tabulated by the source data; no air/vacuum conversion is applied |
| `a`, `b`, `b_b` | m⁻¹ | Volumetric coefficients, not spectral densities |
| Depth, thickness | m | Positive downward from the surface |
| `τ_abs`, component `τ` | — | `a Δz`, `b Δz` |
| Refractive index | — | Real, seawater relative to air |
| Redistribution `emission(p, λ′, λ)` | nm⁻¹ | Photon density per nm of emission wavelength |

A solver working in wavenumber evaluates OceanOptics at `λ = 10⁷/ν`;
per-nm densities additionally need `|dλ/dν| = 10⁷/ν²` and the solver's bin
weights.

## Phase-matrix convention

* **Normalization:** `(1/2) ∫₋₁¹ p(μ) dμ = 1`, i.e. `∫ p dΩ / (4π) = 1`.
* **Expansion:** `p(μ) = Σ_ℓ β_ℓ P_ℓ(μ)` with the `(2ℓ+1)` factor
  included, so `β₀ = 1`.
* **Layout and signs:** [`GreekCoefficients`](@ref) uses the fields and
  convention of vSmartMOM's `GreekCoefs` (Sanghavi 2014 Eq. 16; Hovenier
  signs: `γ₂ = +√(3/2)` for Rayleigh scattering). Relative to
  Mishchenko/VLIDORT-format coefficients, `γ = −β₁` and `ϵ = −β₂`.
* **Angular matrices:** [`evaluate_scattering_matrix`](@ref) and
  [`TabulatedMuellerPhase`](@ref) use the standard convention of Hansen &
  Travis (1974) and Mishchenko et al. (2002), `F₁₂ < 0` for Rayleigh
  scattering at 90°. (vSmartMOM's internal `reconstruct_phase` returns
  `f₁₂` and `f₃₄` with the opposite sign.)
* **Accuracy:** expansions are converged numerically in the scattering
  angle (graded toward the forward peak) to `|Δβ_ℓ|/(2ℓ+1) < 10⁻¹⁰`, and
  [`phase_function_moments`](@ref) throws if a phase function is not
  normalized — nothing is renormalized silently. Analytic forms
  (Rayleigh, Henyey-Greenstein) are exact; Fournier-Forand is
  analytically normalized.
* **Truncation:** expansions are returned through `ℓ_max` without δ-M or
  other truncation; the solver truncates and must request `ℓ_max` at
  least as large as its fitting order.
* **Checks:** [`GreekCoefficients`](@ref) enforce `β₀ = 1`, `|β_ℓ| ≤
  2ℓ+1` and zero `α, γ, ϵ, ζ` for `ℓ < 2`. These catch a doubled `(2ℓ+1)`
  factor and row-shifted tables, not a missing `(2ℓ+1)` factor: compare
  [`asymmetry_parameter`](@ref) with published values when importing.

## Polarization

| Phase function | Phase matrix |
|---|---|
| [`RayleighWaterPhase`](@ref) | closed form (Hansen & Travis 1974) |
| [`TabulatedMuellerPhase`](@ref) | expanded from the tabulated `F(Θ)` |
| [`TabulatedGreekPhase`](@ref) | as tabulated, if all six coefficients are given |
| [`WavelengthDependentPhase`](@ref) | if every entry has one |
| [`FournierForandPhase`](@ref), [`HenyeyGreensteinPhase`](@ref), [`PetzoldPhase`](@ref) | none — scalar only |

A scalar phase function does not determine a phase matrix, so
OceanOptics never invents one: `greek_expansion(pf, ℓ_max; polarized =
true)` and `layer_optics(...; polarized = true)` throw for scalar-only
scatterers unless `scalar_phase_approximation = :depolarizing` is passed,
and the resulting components are labelled `:depolarizing`.

## Zero scattering

Where a layer does not scatter, `ϖ = 0`, the backscatter fraction is 0 and
the mixed expansion is `β = [1, 0, …]` with zero polarization terms.

## Tabulated inputs

Benchmarks prescribe optical properties directly, without the package's
bio-optical models:

* [`TabulatedAbsorber`](@ref)`(λ, a; total = false, source = "")`
* [`TabulatedScatterer`](@ref)`(λ, b, phase; source = "")` — `b_b` follows
  from the phase function, so the two cannot disagree.
* [`TabulatedGreekPhase`](@ref) — given Greek coefficients.
* [`TabulatedMuellerPhase`](@ref) — given `F₁₁, F₁₂, F₂₂, F₃₃, F₃₄, F₄₄`
  on an angular grid covering 0–180°, by keyword; the normalization
  integral is recorded, and rescaling beyond `normalization_tolerance`
  needs `normalize = true`.
* [`WavelengthDependentPhase`](@ref)`(λ, phases; lookup = :exact)` — one
  phase per tabulated wavelength; `lookup = :nearest` uses the nearest
  entry, e.g. one band-representative phase matrix for a fine spectral
  grid. Each distinct entry is expanded once.

Values are defined only at their tabulated wavelengths; any other
wavelength throws. **Total absorption:** published bulk absorption that
already contains water and CDOM (e.g. `a_blk = a_p + a_y + a_w` in
Chowdhary et al. 2020, Eq. 11b) is passed with `total = true`, and
`layer_optics` then rejects any other absorbing constituent.

```julia
# Chowdhary et al. (2020) AOS-II/III ocean body (Table 2)
λ     = [350.0, 450.0, 550.0, 650.0]
a_blk = TabulatedAbsorber(λ, [0.0204, 0.0092, 0.0565, 0.3400]; total = true,
                          source = "Chowdhary et al. 2020, Table 2")
water = TabulatedScatterer(λ, [0.0134, 0.0045, 0.0019, 0.0010],
                           RayleighWaterPhase(depolarization = 0.0);
                           source = "Chowdhary et al. 2020, Table 2")
body  = OceanLayer(0.0, 100.0, AbstractOceanConstituent{Float64}[a_blk, water])
col   = column_optics(OceanColumn([body]), λ; refractive_index = FixedRefractiveIndex(1.34),
                      bottom = BlackBottom(), ℓ_max = 16, polarized = true)
single_scattering_albedo(col[1])  # ≈ [0.39644970, 0.32846715, 0.03253425, 0.002932551]
```

For AOS-IV, the hydrosol matrices of the supplement (Table S1, tabulated
as `F/(4π)` with `F21`, `F43`) enter as `TabulatedMuellerPhase(θ; F11 = 4π .* f11,
F12 = 4π .* f21, F33 = 4π .* f33, F34 = -4π .* f43)`, one per [Chl] value,
combined with `WavelengthDependentPhase`. Their expansions reproduce the
paper's stated asymmetry parameters (0.95, 0.97) and backscatter ratios
(0.0108, 0.0058).

## Refractive index and bottom

`column_optics` requires both choices explicitly:

* `refractive_index = FixedRefractiveIndex(1.34)` (benchmarks) or
  `QuanFry1995(; temperature, salinity)` — Quan & Fry (1995), relative to
  air, fitted to 400–700 nm, 0–30 °C, 0–35 PSU; outside that domain it
  throws unless `extrapolate = true`.
* `bottom = BlackBottom()` or `LambertianBottom(albedo)` (one value or one
  per wavelength). The solver applies the bottom; OceanOptics only
  describes it.

## Validity and provenance

Every interpolated table throws a `DomainError` outside its data range —
nothing is extrapolated silently. Each model documents its source and
range on the [Materials](@ref) page; bundled tables in `data/` carry
provenance headers. The default water absorption, [`IOCCG2018`](@ref),
covers 350–2500 nm; [`SmithBaker1981`](@ref) reaches down to 200 nm but
is several times too high below ~420 nm.

## Inelastic source kernels

Inelastic sources are not part of `layer_optics`; they are kernels the
solver evaluates against its own internal radiation field:

* [`excitation_absorption`](@ref)`(p, λ′)` — photon conversion rate
  `[1/m]`; includes the quantum yield for fluorescence, is the Raman
  scattering coefficient for [`WaterRaman`](@ref). Do not apply the yield
  again.
* [`emission`](@ref)`(p, λ′, λ)` — photon redistribution `[1/nm]`,
  `∫ f dλ = 1`.
* [`radiance_source_kernel`](@ref)`(p, λ′, λ)` — `a^I f^I λ′/λ`, the
  kernel of the energy (radiance) source
  `J(λ) = (1/4π) ∫ b^I (λ′/λ) E₀(λ′) dλ′` for isotropic processes, with
  `E₀` the spectral scalar irradiance (total intensity over solid angle,
  direct and diffuse — not the downwelling plane irradiance).
* [`is_isotropic`](@ref)`(p)` — `WaterRaman` is not: it needs the angular
  excitation field.

The bundled water absorption does not add a separate Raman loss; whether
`b_R(λ′)` is added to the attenuation at the excitation wavelength is the
solver's convention and must match the reference it is compared with.

## The stable contract

* **Solver output:** [`column_optics`](@ref), [`OceanColumnOptics`](@ref),
  [`layer_optics`](@ref), [`OceanLayerOptics`](@ref),
  [`ScatteringComponent`](@ref), [`GreekCoefficients`](@ref) and the
  accessors above; [`BlackBottom`](@ref), [`LambertianBottom`](@ref),
  [`bottom_albedo`](@ref).
* **Constituents:** [`AbstractOceanConstituent`](@ref) and its role
  types; [`absorption`](@ref), [`scattering`](@ref),
  [`backscattering`](@ref), [`scattering_phase`](@ref),
  [`is_total_absorption`](@ref).
* **Phase functions:** [`AbstractOceanPhaseFunction`](@ref),
  [`greek_expansion`](@ref), [`phase_function_moments`](@ref),
  [`is_polarizable`](@ref), [`is_wavelength_dependent`](@ref).
* **Refractive index:** [`AbstractRefractiveIndexModel`](@ref),
  [`refractive_index`](@ref).
* **Inelastic processes:** [`AbstractOceanInelasticProcess`](@ref) and the
  kernels above.

Exported names do not collide with vSmartMOM's exports, so both packages
can be loaded together.
