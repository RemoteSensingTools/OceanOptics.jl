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
| `absorbers` | one [`AbsorptionComponent`](@ref) per absorber, or per absorption part (pigment group) |
| `scatterers` | one [`ScatteringComponent`](@ref) per scattering constituent |
| `inelastic` | one [`InelasticComponent`](@ref) per inelastic process, linked to the absorption that excites it |
| [`absorption_optical_depth`](@ref), [`optical_depth`](@ref), [`scattering_optical_depth`](@ref), [`single_scattering_albedo`](@ref) | derived totals |
| [`greek_expansion`](@ref)`(o, iλ)` | mixed expansion at a wavelength (scattering-weighted) |
| [`backscatter_fraction`](@ref)`(o)` | scattering-weighted backscatter fraction of the phase functions |

An [`AbsorptionComponent`](@ref) carries a `label` (as for scattering
components), a `part` (`:total`, or a pigment group such as
`:chlorophyll`; see [Pigments and fluorescence excitation](@ref)) and its
absorption optical depth `τ = a Δz` per wavelength. The elastic solve
needs only their sum; the split lets the solver report the light absorbed
by each constituent and pigment group.

A [`ScatteringComponent`](@ref) carries a `label` (constituent type, or
the `source` of tabulated data), its scattering optical depth `τ` per
wavelength, its distinct expansions `greek` ([`GreekCoefficients`](@ref);
one unless the phase function depends on wavelength, in which case
`phase_index[iλ]` selects the entry at each wavelength), the backscatter
fraction of its phase function per entry, a `phase_id` key per entry,
and a `polarization` label:

| `polarization` | Meaning |
|---|---|
| `:scalar` | phase function only (`polarized = false`) |
| `:full` | the constituent's actual phase matrix |
| `:depolarizing` | its phase function with the completely depolarizing matrix `diag(p, 0, 0, 0)`, selected explicitly |

Keeping components separate lets a solver build each component's phase
operators once per band and truncate them individually (as it does for
aerosol types), and keeps every phase matrix's provenance visible.
Layers are ordered from the surface down;
[`interface_depths`](@ref)`(col)` gives the layer boundaries.

**Expansion keys.** `phase_id[j]` is the content key of `greek[j]`
([`phase_id`](@ref)`(g)`: a 64-bit FNV-1a hash of the element type and of
the lengths and values of the six coefficient vectors), and
[`phase_id`](@ref)`(c, iλ)` is the key in use at the `iλ`-th wavelength.
Equal keys mean equal coefficients, `ℓ_max` and element type — in any
component, layer, column, session or Julia version — and expansions that
differ in any of these have different keys (up to hash collisions,
probability 2⁻⁶⁴ per pair). A wavelength-dependent component has one key
per distinct expansion. A solver can therefore build its phase operators
once per `(phase_id, Fourier mode)` and store the keys in fixtures. Within
one `column_optics` call, equal phase functions also share the same
`GreekCoefficients` objects.

**Scalar expansions.** With `polarized = false` the expansions carry `β`
only and `α, γ, δ, ϵ, ζ` are empty. A consumer that sizes its arrays by
`length(g.α)` takes [`zero_padded`](@ref)`(g)`: all six vectors of length
`ℓ_max + 1`, the five polarization vectors zero (a polarized expansion is
returned unchanged). Its key differs from that of the scalar expansion.

**Derivatives.** The element type `FT` of constituents, layers and
results is any `Real`, so a layer can be built from ForwardDiff dual
numbers the same way a `Float32` layer is built, and `layer_optics` /
`column_optics` return optical depths, Greek coefficients, backscatter
fractions and refractive indices that carry the derivatives:

```julia
using ForwardDiff
function τ_abs(Chl)                      # layer absorption optical depth
    FT    = typeof(Chl)
    layer = OceanLayer(0.0, 5.0, AbstractOceanConstituent{FT}[
                PureWater{FT}(), CDOM{FT}(), Phytoplankton{FT}(Chl = Chl)])
    return absorption_optical_depth(layer_optics(layer, [443.0, 551.0]; ℓ_max = 16))
end
ForwardDiff.derivative(τ_abs, 0.5)       # ∂τ_a/∂Chl per wavelength
```

All constituents of one layer share the element type (`PureWater{FT}()`
above carries zero derivatives). A constructor without an explicit type,
`Phytoplankton(Chl = x)`, takes `promote_type(Float64, …)` of its real
parameters, so it follows a dual `x`; `Phytoplankton{Float64}(Chl = x)`
with a dual `x` is a `MethodError`.

| differentiable with respect to | through |
|---|---|
| `Chl` (and the `GordonMorel1983` parameters) | `absorption`, `absorption_by_part`, `scattering`, `backscattering`, `layer_optics`, `column_optics` |
| CDOM `a_ref`, `slope`, `λ_ref` | the same |
| NAP `a_ref`, `slope`, `b_ref`, `spectral_exponent`, `c_ref`, `ω_nap`, `B` | the same |
| water `temperature`, `salinity`, `depolarization` | the same; `refractive_index(QuanFry1995(…), λ)` for `n` |
| Fournier-Forand `backscatter_fraction` (or `n`, `μ_j`), Henyey-Greenstein `g`, Rayleigh `depolarization` | `phase_function_value`, `phase_function_moments`, `greek_expansion`, `backscatter_fraction`, and the expansions stored by `layer_optics` |
| layer depths | `layer_optics`, `column_optics` |
| tabulated coefficients of `TabulatedAbsorber` / `TabulatedScatterer`, entries of `IOP` | their accessors and `layer_optics` |
| wavelength | scalar functions; `layer_optics(FT, layer, λ)` with a floating-point layer and a dual `FT` |
| scattering cosine `μ` | `phase_function_value` of Fournier-Forand and Henyey-Greenstein |
| irradiance | [`absorbed_power`](@ref), [`absorbed_photons`](@ref), [`absorbed_par`](@ref) |

The Fournier-Forand backscatter fraction is carried through the
expansion: the constructor refines the root of the inversion `B ↦ (n, μ_j)`
in the number type of `B` by four chord steps with `B′(μ_j)` (each step
makes one more derivative order exact, so nested duals and Hessians up to
fourth order are exact), and the Legendre quadrature is evaluated in
dual numbers. Range checks of the constructors compare values only, so a
dual at a bound (salinity 35 PSU in `QuanFry1995`, `ω_nap = 1`) is
accepted. Element types must be floating-point types or duals of them;
integer and rational types are rejected with an `ArgumentError`. Its refinement criterion looks at the values of the
moments, not at their derivatives; the derivatives of `β_ℓ`, `ℓ ≤ 16`,
agree with central differences to better than `1e-6` relative
(`test/test_ad.jl`).

Not differentiable: inelastic processes ([`WaterRaman`](@ref), the
fluorescence types and their emission spectra) are floating point only,
so a layer of dual element type has no inelastic processes
(`WaterRaman{FT}` with a dual `FT` is a `TypeError`,
`chlorophyll_fluorescence` of a dual `Phytoplankton` a `MethodError`);
their kernels are differentiable with respect to wavelength only.
`PetzoldPhase{FT}()` with a dual `FT` throws an `ArgumentError` (a
measured table without parameters), and `TabulatedMuellerPhase` accepts
floating-point tables only. Derivatives across a table node (wavelength)
are those of the linear segment on one side.

Chlorophyll sensitivities are not a scale factor of a stored optical
depth: the particulate absorption and scattering are non-linear in [Chl]
(`a_p = A(λ) [Chl]^{E(λ)}`, `b_p ∝ [Chl]^{0.62}`). For components that
are linear in their parameter (CDOM in `a_ref`) the derivative equals
the stored `τ` divided by the parameter.

**Element type.** `column_optics(FT, column, λ; …)` and
`layer_optics(FT, layer, λ; …)` return containers of element type `FT`
with an inferred return type; without `FT` the element type of the
column or layer is used.

**Deep water.** Layers have finite thickness. An optically deep ocean is
represented by a bottom layer thick enough for the bottom not to matter
(checked by convergence in depth), with `BlackBottom()`.

## Quantities and units

| Quantity | Units | Notes |
|---|---|---|
| Wavelength `λ` | nm | As tabulated by the source data; no air/vacuum conversion is applied |
| `a`, `b`, `b_b` | m⁻¹ | Volumetric coefficients, not spectral densities |
| Depth, thickness | m | Positive downward from the surface |
| Component `τ` | — | `a Δz` (absorption), `b Δz` (scattering) |
| Scalar irradiance `E₀` | W m⁻² nm⁻¹ | Supplied by the solver; total intensity over all directions |
| Absorbed photons, absorbed PAR | µmol photons m⁻² s⁻¹ (nm⁻¹) | Photon energy `hc/λ` |
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

Values are defined at their tabulated wavelengths. A query matches a
tabulated wavelength when the two agree to a relative `1e-6`, which
covers `Float32` wavelengths and grids computed from wavenumbers as
`1e7/ν`; any other wavelength throws. `TabulatedAbsorber` and
`TabulatedScatterer` accept `lookup = :linear` to interpolate linearly
between neighbouring tabulated wavelengths instead (never beyond the
table) — an explicit choice for spectral grids that do not coincide with
the table. **Total absorption:** published bulk absorption that
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
* `bottom = BlackBottom()` (deep water), `LambertianBottom(albedo)` (one
  value or one per wavelength) or
  [`TabulatedLambertianBottom`](@ref)`(λ, albedo; source)` — an albedo
  spectrum on its own grid (e.g. measured sand or seagrass reflectance),
  interpolated onto `λ`. For these the column carries a `BlackBottom` or
  a `LambertianBottom` with [`bottom_albedo`](@ref)`(col)` per wavelength.
  Any other `AbstractOceanBottom` subtype passes through unchanged for a
  solver that knows it, and `bottom_albedo` throws an `ArgumentError` for
  it.

`FixedRefractiveIndex(n)` requires `n ≥ 1` (`n = 1`: index-matched
interface). The column stores one real index per wavelength; reducing it
to one value per band is the solver's choice.

The solver applies the floor at [`bottom_depth`](@ref)`(col)`, below the
last layer: a Lambertian floor reflects intensity isotropically and
depolarizes (only the azimuthal mode `m = 0`, only `I`). OceanOptics only
describes it.

## Validity and provenance

Every interpolated table throws a `DomainError` outside its data range —
nothing is extrapolated silently. Each model documents its source and
range on the [Materials](@ref) page; bundled tables in `data/` carry
provenance headers. The default water absorption, [`IOCCG2018`](@ref),
covers 350–2500 nm; [`SmithBaker1981`](@ref) reaches down to 200 nm but
is several times too high below ~420 nm.

## Pigments and fluorescence excitation

Phytoplankton can attribute its absorption to pigment groups with a
partition ([`ProspectPigments`](@ref): chlorophyll and carotenoids with
the PROSPECT-D leaf spectra, by analogy with land; or
[`Bidigare1990Pigments`](@ref): chlorophylls a, b, c, photosynthetic and
photoprotective carotenoids and phycoerythrins). The total stays that of
the phytoplankton model; each group becomes an absorption component, and
absorption that the spectra cannot attribute (where all group spectra
vanish) goes to `:unassigned`. Outside the tabulated pigment spectra (below
400 nm) the shares of the nearest tabulated wavelength apply by default;
`outside = UnassignedOutsideTable()` leaves that absorption unassigned.

A fluorescence process is excited through [`ExcitationPathway`](@ref)s,
each an absorber (or absorption part) with an energy-transfer efficiency
`η ∈ [0, 1]`:

`excitation_absorption(p, λ′) = φ · Σ_k η_k · a_k(λ′)`.

Chlorophyll has `η = 1`; carotenoids acting as antenna pass on `η ≤ 1`;
photoprotective absorption and quenching have `η = 0`.
[`chlorophyll_fluorescence`](@ref) builds the pathways from the
partition: `carotenoid_transfer = 0` is chlorophyll-only excitation (as
in Fluspect), `1` a lossless carotenoid antenna, values between the
share passed on; `transfer` sets groups individually. Its default
emission is the chloroplast-level spectrum of Fluspect
([`chloroplast_emission`](@ref)); the quantum yield is always explicit.

In [`layer_optics`](@ref), each process becomes an
[`InelasticComponent`](@ref) whose `excitation` indexes the layer's
`absorbers` and whose `efficiency` holds the `η`s, so that
`quantum_yield · Σ_k efficiency[k] · absorbers[excitation[k]].τ / Δz`
equals `excitation_absorption(process, λ)` on the layer grid. An absorber
that excites fluorescence must be a constituent of the layer: its light
absorption has to attenuate the field it is excited by.

## Light absorbed per layer

Given the solver's layer-mean spectral scalar irradiance
`Ē₀(λ) = (1/Δz) ∫ E₀(z, λ) dz` `[W m⁻² nm⁻¹]` for a layer (total intensity
over all directions, direct and diffuse, integrated through the layer —
not interface values), every component absorbs `A_c(λ) = τ_c(λ) Ē₀(λ)`:

* [`absorbed_power`](@ref)`(o, Ē₀)` — `A_c(λ)` `[W m⁻² nm⁻¹]` per
  absorption component;
* [`absorbed_photons`](@ref)`(o, Ē₀)` — the same in µmol photons
  m⁻² s⁻¹ nm⁻¹;
* [`absorbed_par`](@ref)`(o, Ē₀; range = (400, 700))` — absorbed PAR per
  component `[µmol photons m⁻² s⁻¹]`, e.g. by chlorophyll, per layer;
* [`excitation_optical_depth`](@ref)`(o, j)` — `τᴵ = aᴵ Δz` of the
  layer's `j`-th inelastic process on the layer grid, from the stored
  components;
* [`absorber_index`](@ref)`(o, label, part)` — the row of a component,
  e.g. `absorber_index(o, "Phytoplankton", :chlorophyll)`;
* [`inelastic_emission`](@ref)`(o, j, Ē₀, λ_em)` — spectral power
  `[W m⁻² nm⁻¹]` emitted in all directions through the layer's `j`-th
  inelastic process. For fluorescence, emitted photons equal the quantum
  yield times the photons delivered by the excitation pathways.

Over a column, in the elastic case, `E_d(0⁻) − E_u(0⁻) = Σ_layers Σ_c A_c +`
the floor's absorption `(1 − albedo) E_d(floor)`: a closure test for the
solver.

## Inelastic source kernels

Each inelastic process is carried by its layer's
[`InelasticComponent`](@ref); its kernels are evaluated by the solver
against its own internal radiation field:

* [`conversion_units`](@ref)`(p)` — [`PhotonConversion`](@ref) (fluorescence) or
  [`EnergyConversion`](@ref) (the default [`WaterRaman`](@ref)): what the
  conversion coefficient counts.
* [`excitation_absorption`](@ref)`(p, λ′)` — conversion coefficient `κ`
  `[1/m]`; includes the quantum yield for fluorescence, is the Raman
  scattering coefficient `b_R` for [`WaterRaman`](@ref). Do not apply the
  yield again.
* [`emission`](@ref)`(p, λ′, λ)` — redistribution `f` `[1/nm]` of the
  converted quantity, `∫ f dλ = 1`.
* [`radiance_source_kernel`](@ref)`(p, λ′, λ)` — the energy kernel `K_E`
  (`κ f λ′/λ` for a photon conversion, `κ f` for an energy conversion) of
  the radiance source `J(λ) = (1/4π) ∫ K_E E₀(λ′) dλ′` for isotropic
  processes, with `E₀` the spectral scalar irradiance (total intensity over
  solid angle, direct and diffuse — not the downwelling plane irradiance).
  [`inelastic_coefficient`](@ref) is the photon kernel `K_E λ/λ′`.
* [`is_isotropic`](@ref)`(p)` — `WaterRaman` is not: its source is
  `(1/4π) ∫ K_E [E₀ + β₂ P₂(μ) E₂] dλ′` in the scalar `m = 0` case, with
  the phase function and phase matrix from [`raman_phase`](@ref)`(p)`.

The bundled water absorption does not add a separate Raman loss; whether
`b_R(λ′)` is added to the attenuation at the excitation wavelength is the
solver's convention and must match the reference it is compared with.

## Precision and devices

**Build models in the precision you like, request the element type of the
result.** `layer_optics(FT, layer, λ; …)` and `column_optics(FT, column,
λ; …)` return arrays of element type `FT`; without `FT` the element type
is that of the layer or column. Data tables and physical constants are
`Float64` on the host. A function that reads a table or evaluates a
formula with such constants (pure-water absorption and scattering,
phytoplankton, Fournier–Forand, Quan–Fry, water Raman, the emission
shapes) computes in `Float64` and rounds its result once, so a `Float32`
result is the correctly rounded `Float64` one for the same parameters;
closed forms of the model's own parameters (CDOM, non-algal particles,
Henyey–Greenstein, Rayleigh) compute in the model's precision, within a
few `eps(Float32)` of the `Float64` model.

A scalar function returns the floating-point type that the model and the
argument imply, `promote_type(FT, typeof(float(λ)))`:

```julia
absorption(PureWater{Float32}(), 550f0)   # Float32
absorption(PureWater{Float32}(), 550.0)   # Float64
absorption(PureWater{Float64}(), ForwardDiff.Dual(550.0, 1.0))   # Dual
```

The safest configuration for a single-precision solver is therefore
`Float64` models and `column_optics(Float32, column, λ; …)`: every stored
number is then one rounding away from the double-precision value.

**Results move to a device with `Adapt.adapt`.** The numeric arrays of
`GreekCoefficients`, `AbsorptionComponent`, `ScatteringComponent`,
`OceanLayerOptics`, `LambertianBottom` and `OceanColumnOptics` are type
parameters of these structs:

```julia
using Adapt, CUDA
col = column_optics(Float32, column, λ; refractive_index = QuanFry1995(temperature = 20, salinity = 35), ℓ_max = 32)
dev = adapt(CuArray, col)      # τ, λ, Greek coefficients, B, phase_index, n, albedo on the GPU
adapt(Array, dev)              # back on the host; adapt(Array, col) copies nothing
```

Labels, absorption parts, `phase_id` keys, the lists of components and
layers and the inelastic components (process, excitation indices and
efficiencies) stay on the host: the structs are containers, and a kernel
receives their arrays (`c.τ`, `g.β`, …). `adapt` keeps the element type and
refuses an adaptor that would change it (`CuArray{Float32}` for a `Float64`
column): request the precision from `column_optics`. Expansions that
layers share are moved once and stay shared, and the adapted structs are
not validated again. On
device arrays, construction and validation, `absorption_optical_depth`,
`scattering_optical_depth`, `optical_depth`, `single_scattering_albedo`,
`backscatter_fraction(layer)` and `bottom_albedo` run as broadcasts and
reductions without scalar indexing and agree with the host to `4 eps`.
Use `adapt`, not `cu`: `cu` converts to `Float32` whatever was requested.

Array methods of the constituents (`absorption(c, λ)`, `scattering(c,
λ)`, `backscattering(c, λ)`, `refractive_index(m, λ)` with a vector `λ`)
are one `map`. For the models without data tables (`CDOM`,
`NonAlgalParticles`, `QuanFry1995`, `FixedRefractiveIndex`) they accept a
device array and return one. `QuanFry1995{Float32}` evaluates its
polynomial in `Float64` for each element and rounds once, so this one
needs a device with `Float64` support (CUDA has it, Metal does not); on
other devices evaluate it on the host and `adapt` the result.

A [`SpectralTable`](@ref) is adaptable as well: `adapt(CuArray, tbl)`
moves its two vectors, and `linterp(tbl, λ)` with a device vector `λ`
checks the range once on the host (a `DomainError` outside the table) and
interpolates in one broadcast on the device. The models keep their
`Float64` host tables.

Not available on a device: the models themselves when they read data
tables (`PureWater`, `Phytoplankton`, pigment partitions, tabulated
constituents, emission tables), the Gauss–Legendre quadrature behind
`greek_expansion`, and the functions that take a wavelength index
(`greek_expansion(o, iλ)`, `phase_id(c, iλ)`) or work per component on the
host (`absorbed_power`, `absorbed_photons`, `absorbed_par`,
`excitation_optical_depth`, `inelastic_emission`, `absorber_index`).
Assemble on the host, apply those functions there, then `adapt`.

## Migrating from 0.1

Each breaking change of 0.2.0 with the code that changes. Results that
changed without an API change are listed in the CHANGELOG under "Changed
results".

```julia
# Raman: field/keyword renamed, convention is a type parameter
WaterRaman{Float64}(depolarization = 0.14)            # 0.1
WaterRaman{Float64}(linear_depolarization = 0.17)     # 0.2 (default 0.17, energy conversion)
WaterRaman{Float64}(conversion = PhotonConversion())  # 0.2, the 0.1 photon convention

# A user-defined inelastic process declares its convention
OceanOptics.conversion_units(::MyProcess) = PhotonConversion()   # new in 0.2, required

# Gordon1992 is removed; particulate scattering is a field of Phytoplankton
Phytoplankton{Float64, Gordon1992}(Chl = 0.5)                          # 0.1
Phytoplankton{Float64}(Chl = 0.5, scattering_model = GordonMorel1983()) # 0.2

# Element type is positional, not a keyword
layer_optics(layer, λ; FT = Float32)    # 0.1
layer_optics(Float32, layer, λ)         # 0.2

# Layers take `inelastic`
OceanLayer(0.0, 5.0, cs; fluorophores = [sif])   # 0.1
OceanLayer(0.0, 5.0, cs; inelastic = [sif])      # 0.2

# Quantum yields have no default
chlorophyll_fluorescence(phyto)                          # 0.1
chlorophyll_fluorescence(phyto; quantum_yield = 0.02)    # 0.2

# Pigment partitions: below 400 nm the edge shares apply; the 0.1 rule is
ProspectPigments(carotenoid_ratio = 0.3, outside = UnassignedOutsideTable())

# Photoprotective carotenoids no longer take carotenoid_transfer
chlorophyll_fluorescence(p; quantum_yield = 0.02, carotenoid_transfer = 0.7,
                         photoprotective_transfer = 0.7)   # the 0.1 result
```

Further changes a caller may meet: `layer_optics` throws when a layer lists
`PureWater` twice (merge them into one); `PureWater` and `RayleighWaterPhase`
reject a depolarization outside `[0, 1)`; `pigment_class` returns
`:photoprotective_carotenoid` for the photoprotective carotenoids of
`Bidigare1990Pigments`; integer and rational element types are rejected
(`CDOM{Float64}()` instead of `CDOM{Int}()`); positional inner constructors
of `Phytoplankton`, `NonAlgalParticles`, `ProspectPigments`,
`Bidigare1990Pigments` and `WaterRaman` have the new fields (prefer the
keyword constructors, which are unchanged). The CHANGELOG lists every
change with its reason.

## The stable contract

* **Solver output:** [`column_optics`](@ref), [`OceanColumnOptics`](@ref),
  [`layer_optics`](@ref), [`OceanLayerOptics`](@ref),
  [`AbsorptionComponent`](@ref), [`ScatteringComponent`](@ref),
  [`InelasticComponent`](@ref), [`GreekCoefficients`](@ref) and the
  accessors above; [`phase_id`](@ref), [`zero_padded`](@ref);
  [`BlackBottom`](@ref), [`LambertianBottom`](@ref),
  [`TabulatedLambertianBottom`](@ref), [`bottom_albedo`](@ref).
* **Absorbed light:** [`absorbed_power`](@ref), [`absorbed_photons`](@ref),
  [`absorbed_par`](@ref), [`inelastic_emission`](@ref),
  [`excitation_optical_depth`](@ref), [`absorber_index`](@ref).
* **Constituents:** [`AbstractOceanConstituent`](@ref) and its role
  types; [`absorption`](@ref), [`absorption_parts`](@ref),
  [`absorption_by_part`](@ref), [`scattering`](@ref),
  [`backscattering`](@ref), [`scattering_phase`](@ref),
  [`is_total_absorption`](@ref); pigment partitions
  ([`AbstractPigmentPartition`](@ref)).
* **Phase functions:** [`AbstractOceanPhaseFunction`](@ref),
  [`greek_expansion`](@ref), [`phase_function_moments`](@ref),
  [`is_polarizable`](@ref), [`is_wavelength_dependent`](@ref).
* **Refractive index:** [`AbstractRefractiveIndexModel`](@ref),
  [`refractive_index`](@ref).
* **Inelastic processes:** [`AbstractOceanInelasticProcess`](@ref),
  [`excitation_pathways`](@ref) and the kernels above.

Exported names do not collide with vSmartMOM's exports, so both packages
can be loaded together.
