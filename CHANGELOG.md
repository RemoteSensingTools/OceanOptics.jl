# Changelog

All notable changes to OceanOptics.jl are recorded here. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and
this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.0] — 2026-11-13

Release notes. Water Raman is an energy conversion by default and its
phase function follows Mobley et al. (1993); every physical parameter is a
field or keyword with a documented default (see the parameter reference);
model parameters can be differentiated with ForwardDiff through
`layer_optics` and `column_optics`; light absorbed per constituent,
chlorophyll fluorescence with chloroplast-level emission, and the ocean
floor were added for the coupled solver. Each entry below says what a user
must change; the "Migrating from 0.1" section of the Consumer API page shows
the code.

Light absorbed per constituent, chlorophyll fluorescence with
land-analogous chloroplast emission, and the ocean floor, for the
fluorescence and absorbed-PAR requirements of the coupled solver.

### Changed results

- **Photoprotective carotenoids no longer excite chlorophyll fluorescence by
  default.** `pigment_class` gives them their own class and
  `chlorophyll_fluorescence` the keyword `photoprotective_transfer = 0`;
  `carotenoid_transfer` now applies to photosynthetic (and Prospect's lumped)
  carotenoids only. With `Bidigare1990Pigments(chl_a = 0.5, chl_c = 0.1,
  psc = 0.3, ppc = 0.15)` and `carotenoid_transfer = 0.7` the excitation
  coefficient falls by the photoprotective share, 18 % at 440 nm and 12 % at
  360 nm (physics review DREV_PHYS). To keep the 0.1 result pass
  `photoprotective_transfer` equal to `carotenoid_transfer`.

- **Bricaud absorption of `Phytoplankton{Bricaud1998}`: `a = A(λ)·Chl^E(λ)`**
  (formerly `A(λ)·Chl^(1−E(λ))`). The bundled tables are HydroLight's
  `AEmidUVabs.txt`, whose header defines the particulate absorption
  `ap(wave) = A(wave)*[chl]^E(wave)` (Bricaud et al. 1998,
  `a_p(440) = 0.052 Chl^0.635`); the exponent `1 − E` belongs to the
  specific-absorption form of Bricaud et al. (1995), which has other
  coefficients. Results are unchanged at Chl = 1 mg m⁻³ and change
  everywhere else:

  | Chl [mg m⁻³] | a(440) before | a(440) now | a(676) before | a(676) now |
  |---|---|---|---|---|
  | 0.03 | 0.01446 | 0.00561 | 0.01089 | 0.00114 |
  | 0.1 | 0.02244 | 0.01205 | 0.01348 | 0.00306 |
  | 0.5 | 0.04038 | 0.03348 | 0.01795 | 0.01148 |
  | 1 | 0.05200 | 0.05200 | 0.02030 | 0.02030 |
  | 10 | 0.12050 | 0.22439 | 0.03056 | 0.13483 |

  `Chl = 0` gives zero absorption (formerly `Inf` at 698–700 nm, where
  `E > 1`). Everything derived from this absorption changes by the same
  factors: the pigment parts (`absorption_by_part`, the absorption
  components of `layer_optics`), the excitation of
  `chlorophyll_fluorescence`, `absorbed_power`, `absorbed_photons` and
  `absorbed_par`, and the documentation figures at Chl ≠ 1. The quantity is
  the particulate absorption `a_p` (phytoplankton plus detritus), not
  `a_φ`: the docstrings, data-file headers and documentation say so, and
  warn that a separate `NonAlgalParticles` constituent adds detrital
  absorption a second time.
- **`NonAlgalParticles{FT, FellMineral}` defaults are those of Fell (1997),
  p. 120:** single-scattering albedo `ω_nap = 0.95` (formerly 0.96) and
  backscatter fraction `B = 0.036` (formerly 0.03). For the default
  `c_ref = 0.4 m⁻¹`: `a(440)` 0.016 → 0.020 m⁻¹ (+25 %), `b(440)` 0.384 →
  0.380 m⁻¹, `b_b(440)` 0.01152 → 0.01368 m⁻¹ (+19 %). Explicit `ω_nap`
  and `B` give the same results as before.
- The particulate scattering law `b_p = 0.30·Chl^0.62·(550/λ)` is unchanged
  in value; it is now attributed to Gordon & Morel (1983) (as cited by
  Fell 1997, Eq. 5.5) instead of Morel & Maritorena (2001), and described
  as the total particulate scattering of Case-1 water.
- `TabulatedAbsorber` and `TabulatedScatterer` with linear lookup return
  the type that the query wavelength and the table promote to (formerly
  always the table's precision): a `Float64` query of a `Float32` table
  gives a `Float64`, for example 0.0175 at 425 nm between 0.02 (400 nm)
  and 0.01 (500 nm) instead of `0.0175f0`; values of
  `Float64` tables are unchanged. A `ForwardDiff.Dual` wavelength now
  gives the slope of the segment (formerly a `MethodError`).

- **Water Raman source: `b_R` is an energy conversion coefficient.**
  `radiance_source_kernel(::WaterRaman, λ′, λ)` is now `b_R(λ′) f_R(λ′, λ)`
  (formerly `b_R f_R λ′/λ`), the definition of Mobley et al. (1993, Problem 7
  and Eqs. A3–A4): `∫ K_E dλ = b_R(λ′)`. The Raman energy kernel and
  `inelastic_emission` of a Raman component rise by the inverse of the band
  mean of `λ′/λ`: ×1.135 at 350 nm excitation, ×1.165 at 417 nm (the Mobley
  Problem 7 source `b_R E₀/(4π)` was 14.2 % low), ×1.199 at 488 nm, ×1.256
  at 600 nm, ×1.31 at 700 nm. The photon kernel
  `inelastic_coefficient(::WaterRaman, λ′, λ)` rises by the same factors
  (it is now `K_E λ/λ′`, formerly `b_R f_R`). The former behaviour is
  `WaterRaman(conversion = PhotonConversion())`. Fluorescence kernels are
  unchanged.
- **Default Raman depolarization 0.17** (formerly 0.14, a field that no
  function used). It now sets the Raman phase function `raman_phase`.

- **Pigment partitions outside their tables hold the edge shares.**
  `pigment_fractions` of `ProspectPigments` and `Bidigare1990Pigments` below
  400 nm returns the shares at 400 nm (formerly all `:unassigned`); above the
  red end of the tables the shares are all `:unassigned` under both rules, so
  results there are unchanged. With
  Chl = 0.5 mg m⁻³, `excitation_range = (350, 700)` and unit transfer
  efficiencies, `excitation_absorption` of `chlorophyll_fluorescence` at
  399.999 nm is 5.32·10⁻⁴ m⁻¹ (formerly 0) and its 350–700 nm integral rises
  by 1/0.781 = 1.281 to the unpartitioned value; the pigment rows of
  `absorbed_power`, `absorbed_photons` and `layer_optics` below 400 nm change
  accordingly. Results at 400–750 nm, and with the default
  `excitation_range = (400, 700)`, are unchanged. The former rule is
  `outside = UnassignedOutsideTable()`.

- **`PureWater(depolarization = δ)` enters the Zhang–Hu–He scattering
  magnitude** (formerly fixed at 0.039 inside the routine; only the phase
  matrix used the field). Unchanged for the default `δ = 0.039` in
  `Float64`; for `δ = 0.09`, `b(500 nm, 20 °C, S = 35)` is 2.785·10⁻³ m⁻¹
  (formerly 2.547·10⁻³, +9.3 %); for `δ = 0`, 2.385·10⁻³ (−6.4 %). In
  `Float32` the default value is now the stored `0.039f0`, which changes
  the last bit of some results (2 of 55 wavelengths at 5 °C, S = 10).

### Breaking

- **Integer and rational element types are rejected** with an
  `ArgumentError` (`CDOM{Int}()`, `layer_optics(Rational{Int}, …)`; formerly an
  `InexactError` or a silent computation in rationals). Use a floating-point
  type. `IOP` of integer or rational arrays converts them to `Float64`, as the
  scalar `IOP` does.
- **`Phytoplankton{FT}` rejects a scattering model of a non-floating-point
  element type that `FT` cannot hold** (a dual `GordonMorel1983` in a
  `Phytoplankton{Float32}`); construct it without `FT`, which then follows the
  model's element type.
- **`pigment_class(::Bidigare1990Pigments, :photoprotective_carotenoids)`** is
  `:photoprotective_carotenoid` (formerly `:carotenoid`).

- **`Gordon1992` is removed.** The tag had no absorption (it threw "not yet
  implemented") and its scattering was the generic Gordon & Morel law, not
  Gordon's. Use `Phytoplankton{FT, Bricaud1998}`.
- **`layer_optics` rejects a layer that lists `PureWater` more than once**
  (`ArgumentError`); formerly the water was counted twice. Other
  constituents may repeat.
- **`PureWater` validates `depolarization ∈ [0, 1)`.**
- **`Phytoplankton` has a fifth type parameter and a field
  `scattering_model`** (default `GordonMorel1983()`), and
  `NonAlgalParticles` a field `spectral_exponent` (default 1); positional
  constructors changed accordingly.

- **`ProspectPigments` and `Bidigare1990Pigments` have a second type
  parameter and a field `outside`** (`HoldEdgeShares()` by default,
  `UnassignedOutsideTable()` for the former behaviour); the positional inner
  constructors changed accordingly.

- **`WaterRaman` has a second type parameter and a renamed field.**
  `WaterRaman{FT, U<:AbstractConversionUnits}` carries the convention of
  `b_R` in the new field `conversion` (default `EnergyConversion()`);
  `WaterRaman{Float64}` is no longer a concrete type (`WaterRaman{Float64}()`
  still constructs). The field and keyword `depolarization` is now
  `linear_depolarization`: the depolarization ratio for linearly polarized
  excitation (Mobley et al. 1993, Eq. 9), limited to `[0, 3/4]`.
- **Every inelastic process must define `conversion_units`.** A subtype of
  `AbstractOceanInelasticProcess` without a `conversion_units` method throws
  a `MethodError` in `inelastic_coefficient`, `radiance_source_kernel` and
  `inelastic_emission`.

- `layer_optics` and `column_optics` have no `FT` keyword: the element
  type is the optional first positional argument,
  `layer_optics(Float32, layer, λ; …)`. Passing `FT = …` is a
  `MethodError`.
- `InelasticComponent{FT, P}` is host data: `efficiency` is a
  `Vector{FT}` next to `excitation::Vector{Int}`, and `Adapt.adapt`
  returns the component unchanged (the layer's `inelastic` vector is the
  same object after `adapt`). There is no array type parameter.
- The table lookup is a type: `TabulatedAbsorber{FT, L}`,
  `TabulatedScatterer{FT, P, L}` and `WavelengthDependentPhase{FT, P, L}`
  with `lookup::L`, one of the new exported `ExactLookup`, `LinearLookup`
  and `NearestLookup` (`<: AbstractTableLookup`). The keyword `lookup`
  still accepts `:exact`, `:linear` and `:nearest`, and also the lookup
  objects; the field no longer holds a `Symbol`, and the inner
  constructors are `TabulatedAbsorber{FT, L}(λ, a, total, source, lookup)`,
  `TabulatedScatterer{FT, P, L}(…)` and
  `WavelengthDependentPhase{FT, P, L}(λ, phases, lookup)`.
- `Adapt.adapt` keeps the element type of every result struct and of
  `SpectralTable`, and refuses an adaptor that changes it
  (`adapt(CuArray{Float32}, col64)`, `adapt(Array{Float32}, g64)`) with an
  `ArgumentError` that says to build the optics in the wanted precision;
  formerly `GreekCoefficients` and `SpectralTable` converted and layers
  and columns failed with a `TypeError`.
- `OceanColumnOptics(layers, n, bottom)` rejects a refractive index below 1
  (formerly any positive value), with the message of
  `FixedRefractiveIndex`.
- `SpectralTable(λ, y)` rejects non-finite values `y` (formerly only
  non-finite wavelengths).
- The element-type parameter of the result structs, of `OceanLayer`,
  `OceanColumn` and the bottoms is bounded by `Real`:
  `LambertianBottom{Any, Vector{Any}}` no longer constructs.
- The column name of `data/bricaud_1998_A.csv` is `A_per_m` (formerly
  `A_m2_per_mg`).

- `OceanLayerOptics` carries `absorbers` (one `AbsorptionComponent` per
  absorber or absorption part: `label`, `part`, `τ`), `scatterers` (the
  former `components`) and `inelastic` (one `InelasticComponent` per
  process, with the indices of the absorption components that excite it
  and their energy-transfer efficiencies). The `τ_abs` field is gone;
  `absorption_optical_depth(o)` sums the absorption components.
- `OceanLayer` and the column helpers take `inelastic` (formerly
  `fluorophores`), since water Raman is not a fluorophore.
- `IsotropicFluorescence(excitation; emission, quantum_yield,
  excitation_range)` is driven by `ExcitationPathway`s (absorber, part,
  efficiency η): `a^I = φ Σ_k η_k a_k`. Quantum yields have no default
  (`chlorophyll_fluorescence`, `CDOMFluorescence`, `IsotropicFluorescence`);
  previously 0.003 and 0.01 were implied.
- `chlorophyll_fluorescence` defaults to the chloroplast-level emission of
  Fluspect instead of a 685-nm Gaussian (still available through
  `emission = GaussianEmission(...)`).
- `layer_optics` throws if a fluorescence process is excited by an
  absorber that is not a constituent of the layer.
- `ScatteringComponent.phase_id` is a `Vector{UInt64}` with one content
  key per entry of `greek` (formerly one `UInt64` derived from the
  identity of the phase-function object). The key is the 64-bit FNV-1a
  hash of the element type and coefficients, so it is equal for equal
  expansions, differs with `ℓ_max`, element type, polarization treatment
  and wavelength-dependent entry, and is reproducible between sessions.
  The inner constructor `ScatteringComponent{FT}(label, τ, greek, B,
  phase_index, polarization)` computes the keys; the `phase_id` argument
  and keyword are gone.
- `layer_optics` and `column_optics` default to the element type of the
  layer or column (formerly `Float64`); this changes the result only for
  `Float32` layers and columns.
- `InelasticComponent{FT, P}` has the process type as a second parameter
  (`process::P`). `InelasticComponent{FT}(label, process, excitation,
  efficiency)` still constructs it.
- `TabulatedAbsorber` and `TabulatedScatterer` have a `lookup` field;
  their inner constructors take it as the last argument.
- Exact-wavelength tables (`TabulatedAbsorber`, `TabulatedScatterer`,
  `TabulatedGreekPhase`, `WavelengthDependentPhase(lookup = :exact)`)
  match wavelengths to a relative tolerance of `1e-6` (formerly 8 ulp of
  the less precise operand), and reject tables whose neighbouring
  wavelengths are closer than twice that.
- `FixedRefractiveIndex(n)` rejects `n < 1` (formerly any positive `n`
  was accepted); `n = 1` is allowed.
- `PureWater` with `IOCCG2018` absorption throws a `DomainError` naming
  the composite range 350–2500 nm (formerly the range of the sub-table).
- `inelastic_emission` returns a vector whose element type is the promoted
  type of the layer, the irradiance and the emission wavelengths
  (formerly `Float64` for a `Float32` layer).
- `SpectralTable{T, V}` holds `λ`, `y` and the range limits `lo`, `hi` and
  is interpolated by the package itself; the `itp` field and the
  Interpolations.jl dependency are gone. Interpolated values are
  bit-identical to 0.1.0 on all shipped tables (checked at every tabulated
  wavelength and 40 000 others per table; for example pure-water
  absorption at 551.3 nm, 0.05732191899999997 1/m, unchanged), and the
  segment that a tabulated wavelength is assigned to is the same, so
  derivatives are unchanged as well. `SpectralTable(λ, y)` accepts any real vectors
  and rejects repeated wavelengths ("strictly increasing"; formerly
  "sorted ascending" and `Vector{T}` arguments only).
- The result structs carry their array types as type parameters, appended
  after the existing ones: `GreekCoefficients{FT, V}`,
  `AbsorptionComponent{FT, V}`, `ScatteringComponent{FT, V, I}`,
  `OceanLayerOptics{FT, V, I}`,
  `LambertianBottom{FT, V}` and `OceanColumnOptics{FT, B, V, I}` (`V` the
  array type of the floating-point fields, `I` that of `phase_index`).
  `OceanLayerOptics{FT}`, `GreekCoefficients{FT}` and
  `OceanColumnOptics{FT, BlackBottom}` in a signature dispatch as before,
  and the `{FT}` constructors still work. No longer concrete, and to be
  replaced where a concrete type is needed (struct fields, `Vector`
  element types, `===` on types): `OceanLayerOptics{FT}` (now
  `OceanLayerOptics{FT, Vector{FT}, Vector{Int}}` on the host),
  `GreekCoefficients{FT}`, `LambertianBottom{FT}` and therefore
  `OceanColumnOptics{FT, LambertianBottom{FT}}`, which matches no column:
  write `OceanColumnOptics{FT, <:LambertianBottom}`.
- The inner constructors take the array types
  (`GreekCoefficients{FT, V}(α, β, γ, δ, ϵ, ζ)` and likewise for the other
  result structs) and store their arguments without converting them. The
  `{FT}` constructors keep a dense vector of element type `FT` as it is
  (formerly always a `Vector{FT}` copy made by `convert`), convert one of
  another element type where it lives, and collect anything else; the
  keyword and positional convenience constructors still copy.
- `phase_function_value` of `RayleighWaterPhase{Float32}`,
  `PetzoldPhase{Float32}` and a `Float32` `TabulatedGreekPhase` returns
  `Float64` for a `Float64` argument (formerly `Float32`: the argument was
  converted to the model's type); a `Float32` argument still gives
  `Float32`.
- Scalar functions of `Float32` models return `Float32` for `Float32`
  arguments (formerly `Float64`): `absorption`, `scattering` and
  `backscattering` of `PureWater` and `Phytoplankton`,
  `phase_function_value` of `FournierForandPhase`, `emission` and
  `raman_peak_wavelength` of `WaterRaman`, `redistribution` of
  `GaussianEmission` and `WavenumberShiftEmission`; so do
  `absorbed_photons` and `absorbed_par` of a `Float32` layer, and the
  vector methods. The result
  type is the model's element type promoted with the argument's, so a
  `Float64` or dual-number argument still gives a `Float64` or dual result.
  These functions compute in `Float64` and round once: the `Float32` value
  is the correctly rounded `Float64` value for the same parameters.
  `Float64` results are unchanged.

### Added

- Keyword `photoprotective_transfer` (default 0) of
  `chlorophyll_fluorescence`.
- A "Migrating from 0.1" section on the Consumer API page.

- **Derivatives with respect to model parameters by ForwardDiff.** The
  element type of the elastic path is any `Real` (formerly
  `AbstractFloat`): `AbstractOceanConstituent{FT}`,
  `AbstractOceanPhaseFunction{FT}`, `CDOM`, `Phytoplankton`,
  `GordonMorel1983`, `NonAlgalParticles`, `PureWater`,
  `FournierForandPhase`, `HenyeyGreensteinPhase`, `RayleighWaterPhase`,
  `TabulatedGreekPhase`, `TabulatedAbsorber`, `TabulatedScatterer`,
  `FixedRefractiveIndex`, `QuanFry1995`, `IOP`, and
  `layer_optics(FT, …)` / `column_optics(FT, …)`. A layer built at a dual
  element type gives optical depths, single-scattering albedo, Greek
  coefficients, backscatter fractions and refractive index that carry
  the derivative with respect to chlorophyll, CDOM `a_ref` and `slope`,
  the non-algal particle parameters, water temperature, salinity and
  depolarization, layer thickness and the Fournier-Forand backscatter
  fraction (the inversion `B ↦ (n, μ_j)` is linearized at its root and
  the Legendre quadrature runs in dual numbers). ForwardDiff is not a
  dependency. Floating-point results are unchanged bit for bit. Not
  covered: inelastic processes and their emission spectra,
  `PetzoldPhase` (an `ArgumentError` for a non-floating-point element
  type) and `TabulatedMuellerPhase`. See "Forward-mode automatic
  differentiation" in the worked example and "Derivatives" in the
  consumer interface; `test/test_ad.jl` checks every parameter against
  central differences (rtol 1e-6) or closed forms.
- Constructors without an explicit element type (`CDOM(…)`,
  `Phytoplankton(…)`, `NonAlgalParticles(…)`, `PureWater(…)`,
  `GordonMorel1983(…)`, `FournierForandPhase(…)`,
  `HenyeyGreensteinPhase(…)`, `RayleighWaterPhase(…)`) take
  `promote_type(Float64, …)` of their real parameters: `Float64` for
  plain numbers of any precision as before, the dual type for a dual
  parameter. A `BigFloat` parameter now gives a `BigFloat` model
  (formerly rounded to `Float64`).
- `IOP(a, b, bb)` accepts dual numbers; three equal-typed arrays of any
  `Real` element type are accepted (formerly floating point only).
- Parameter exposure: `GordonMorel1983(b_ref = 0.30, chl_exponent = 0.62,
  λ_ref = 550, spectral_exponent = 1)` as `Phytoplankton(scattering_model = …)`;
  `NonAlgalParticles(spectral_exponent = 1)`; keyword `roundoff_rtol` of
  `absorbed_power`, `absorbed_photons`, `absorbed_par` and
  `inelastic_emission`; the pure-water depolarization 0.039, the PAR limits
  (400, 700) nm and the O–H stretch band (3400, 85 cm⁻¹) are each defined
  once. Defaults reproduce the former `Float64` results bit for bit.
- `HoldEdgeShares` and `UnassignedOutsideTable` (`outside` keyword of the
  pigment partitions).
- Documentation page "Parameters" (`docs/src/reference/parameters.md`): every
  physical parameter with default, unit, source and how to set it.

- `AbstractConversionUnits`, `PhotonConversion`, `EnergyConversion` and
  `conversion_units(process)`: whether the conversion coefficient of an
  inelastic process counts photons or energy; the energy and photon kernels
  dispatch on it.
- `raman_phase(::WaterRaman)` (the Rayleigh-type phase function and phase
  matrix of Raman emission, `β₂ = (1−ρ)/(2(1+2ρ))`) and
  `natural_depolarization(::WaterRaman)` (`ρ_n = 2ρ/(1+ρ)`).

- `layer_optics(FT, layer, λ; …)` and `column_optics(FT, column, λ; …)`:
  positional element type, with an inferred return type.
- `phase_id(g::GreekCoefficients)` and `phase_id(c::ScatteringComponent, iλ)`:
  the content key of an expansion.
- `zero_padded(g)`: a scalar expansion with `α, γ, δ, ϵ, ζ` as zero
  vectors of length `ℓ_max + 1`, for consumers that size by `length(α)`.
- `TabulatedAbsorber(...; lookup = :linear)` and
  `TabulatedScatterer(...; lookup = :linear)`: opt-in linear interpolation
  between tabulated wavelengths (never extrapolation).
- `absorbed_power`, `absorbed_photons`, `absorbed_par` and
  `inelastic_emission` treat negative irradiance entries down to
  `−√eps · maximum(abs, E₀)` as zero; more negative entries still throw.
- `bottom_albedo` throws an `ArgumentError` naming the bottom type for
  bottoms other than `BlackBottom` and `LambertianBottom` (formerly a
  `MethodError`).
- Device-adaptable results: `Adapt.adapt(to, x)` for `GreekCoefficients`,
  `AbsorptionComponent`, `ScatteringComponent`, `InelasticComponent`,
  `OceanLayerOptics`, `LambertianBottom` and `OceanColumnOptics` moves the
  numeric arrays (for example `adapt(CuArray, column)`) and keeps the
  element type; labels, parts, `phase_id` keys, the component and layer
  lists and the inelastic process models stay on the host. Construction
  and validation of these structs, `absorption_optical_depth`,
  `scattering_optical_depth`, `optical_depth`, `single_scattering_albedo`,
  `backscatter_fraction(layer)` and `bottom_albedo` use reductions and
  broadcasts only, so they run on device arrays with scalar indexing
  disallowed. Adapt.jl is a new dependency.
- `linterp(tbl, λ::AbstractVector)`: one range check and one broadcast.
  `Adapt.adapt(to, tbl)` moves a `SpectralTable` to a device, where the
  vector method interpolates device wavelengths; an out-of-range
  wavelength throws the `DomainError` on the host.
- Array methods of `absorption`, `scattering`, `backscattering` and
  `refractive_index` are a `map` over the wavelengths (same host results);
  for `CDOM`, `NonAlgalParticles`, `QuanFry1995` and `FixedRefractiveIndex`
  they accept and return device arrays. `QuanFry1995{Float32}` evaluates
  its polynomial in `Float64` per element and rounds once, so its array
  method needs a device with `Float64` support.
- `adapt` builds the adapted result structs without validating them again
  and keeps expansions that components and layers share shared (one device
  array per distinct host expansion), so `GreekCoefficients` and
  `LambertianBottom` can be captured by a device broadcast (the components,
  layers and columns carry labels and lists and remain host containers).
- `hash(::GreekCoefficients)`, consistent with `==`, so expansions work as
  `Dict` keys and in `unique`.
- Tests: Aqua, ExplicitImports and JET (`test_package`, and `test_opt` /
  `test_call` on the scalar functions and expansions), `@inferred` checks
  including `layer_optics`, `column_optics`, `excitation_optical_depth`
  and `inelastic_emission`; `Float32` tests of the result types and of
  `Float32` against `Float64` with the stated difference; device tests
  (`test/test_device.jl`) in the default suite on JLArrays 0.5 with scalar
  indexing disallowed, and the same file on CUDA from `test/gpu`.
- Absorption parts: `absorption_parts`, `absorption_by_part`,
  `absorption_part`; `Phytoplankton(...; pigments)` with the pigment
  partitions `ProspectPigments` (chlorophyll and carotenoids with the
  PROSPECT-D leaf spectra used by Fluspect; land analogy) and
  `Bidigare1990Pigments` (chlorophylls a, b, c, photosynthetic and
  photoprotective carotenoids, phycoerythrins; Bidigare et al. 1990). The
  total absorption is unchanged; unattributable absorption is
  `:unassigned`.
- `chlorophyll_fluorescence(phyto; quantum_yield, carotenoid_transfer,
  phycobilin_transfer, transfer)`: chlorophyll-only excitation
  (`carotenoid_transfer = 0`), carotenoids as lossless antenna (`1`), or
  any share in between, per pigment group if needed.
- `TabulatedEmission` and `chloroplast_emission([FT]; psi_share)`:
  Fluspect's combined chloroplast emission, or a PSI/PSII mixture.
- `absorbed_power`, `absorbed_photons`, `absorbed_par` and
  `inelastic_emission`: light absorbed per layer and component, and the
  emission it drives, from the solver's layer-mean scalar irradiance;
  `excitation_optical_depth(o, j)` and `absorber_index(o, label, part)`.
- `TabulatedLambertianBottom(λ, albedo; source)`: a floor albedo spectrum
  on its own grid, interpolated onto the column grid by `column_optics`.
- Data: `fluspect_chloroplast_emission.csv`, `prospect_d_pigment_absorption.csv`
  (both from SCOPE's Optipar2017_ProspectD) and `bidigare_1990_pigments.csv`
  (Ocean Optics Web Book, CC BY).

### Fixed

- Range checks accept a dual number at their bounds. ForwardDiff 1.x
  orders duals by value and then partials, so `Dual(35, 1) ≤ 35` is false and
  `∂n/∂S` of `QuanFry1995` at 35 PSU, `∂n/∂T` at 30 °C, `∂/∂λ` at a table end
  and `∂b/∂ω_nap` at `ω_nap = 1` threw; the checks now compare values.
- `∂b_w/∂S` of `ZhangHu2009` scattering at `S = 0` is the one-sided
  derivative (1.764·10⁻⁵ m⁻¹ PSU⁻¹ at 500 nm, 20 °C) instead of `NaN`.
- Second and higher derivatives with respect to the Fournier-Forand
  backscatter fraction: `d²μ_j/dB²` was 0 and `d²β₂/dB²` `NaN` without a
  warning; the inversion is refined by chord steps in the dual type
  (derivatives exact to fourth order) and `δ^ν` is evaluated as
  `exp(ν log δ)` for duals, whose second derivative overflowed at the forward
  quadrature floor. Floating-point results are unchanged.
- `Phytoplankton(; scattering_model = GordonMorel1983(b_ref = dual))` without
  an element type builds a dual constituent (formerly a `Float64` one that
  failed when evaluated).

- `excitation_optical_depth` and `inelastic_emission` have inferred return
  types; `pigment_fractions(::Bidigare1990Pigments, λ)` is inferred on
  Julia 1.10; `greek_expansion(::PetzoldPhase{Float32}, ℓ)` runs without
  runtime dispatch (the phase function returned `Float64` or `Float32`
  depending on the angle).
- Static-analysis findings: an unbound type parameter in
  `_pigment_shares`, a missing docstring for `quantum_yield`, implicit
  imports from DocStringExtensions, and four JET reports on Julia 1.12.
- Documentation no longer claims that model parameters can be
  differentiated with ForwardDiff or that optical depths are linear in
  chlorophyll.
- Docstrings that disagreed with the code (no function changed):
  `absorption` does not accept `Unitful` quantities (wavelengths are
  numbers in nm, or vectors of them); `raman_peak_wavelength` for 440 nm
  excitation and the default shift is 517.4 nm, not 491 nm; the bundled
  `Kou1993` table has steps of 0.3 to 18 nm, not about 1 nm; `ZhangHu2009`
  checks no range of validity; `IOP` is a standalone value type that `layer_optics`
  does not use; `is_isotropic(::WaterRaman)` refers to the process's own
  depolarization ratio.
- Every exported name has a docstring with a signature line (or type
  definition) and an imperative first sentence, at most 92 characters per
  line; `absorption`, `raman_peak_wavelength` and `fwhm_to_sigma` carry
  doctests.

### Packaging

- `CI.yml` runs the tests only; `Documentation.yml` is the one docs build:
  on pull requests it builds without deploying, on pushes to `main` and tags
  it builds and deploys. JET 0.12 stays excluded from `[compat]`: it does not
  load on Julia 1.12 with the current LoweredCodeUtils
  (`record_termination_points!` undefined); 0.9 (Julia 1.10) and 0.13 (1.12)
  are used.

- Runtime dependencies are Adapt, DocStringExtensions and
  FastGaussQuadrature; Interpolations is dropped. JLArrays 0.5 is a test
  dependency. CUDA is a dependency of `test/gpu` only.
- The documentation build is strict: a missing docstring of an exported
  name, an unresolved cross-reference or a failing doctest is an error
  (`warnonly` removed from `docs/make.jl`). `Documentation.yml` runs it.

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
