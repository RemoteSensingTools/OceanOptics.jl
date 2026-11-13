# vSmartMOM integration

How [vSmartMOM.jl](https://github.com/RemoteSensingTools/vSmartMOM.jl)
consumes OceanOptics. The adapter below was run against the vSmartMOM
release candidate `d441061` (October 2026) in a joint environment; it is
documentation, not part of either package.

## Compatibility

* OceanOptics' dependency bounds include the versions vSmartMOM resolves
  to (`FastGaussQuadrature` 0.5 through its CanopyOptics dependency,
  `Adapt` 4); the OceanOptics test suite passes with them.
* No exported names collide, so `using vSmartMOM, OceanOptics` works
  without qualification. (vSmartMOM's own `greek_coefficients`,
  `phase_function` and `scattering_matrix` refer to its analytic phase
  functions; OceanOptics uses [`greek_expansion`](@ref),
  [`scattering_phase`](@ref) and [`evaluate_scattering_matrix`](@ref).)
* OceanOptics v0.1.0 is in the General registry's queue; v0.2.0 (this
  documentation) adds absorption and inelastic components and breaks the
  v0.1 layer layout. Declare the version the adapter is written against in
  vSmartMOM's `[weakdeps]`/`[compat]`.

## Mapping

| OceanOptics | vSmartMOM |
|---|---|
| [`GreekCoefficients`](@ref) `(α, β, γ, δ, ϵ, ζ)` | `Scattering.GreekCoefs(α, β, γ, δ, ϵ, ζ)` — same fields, order and signs |
| [`ScatteringComponent`](@ref) | one `CoreScatteringOpticalProperties` per component, `ϖ = 1`, truncated like an aerosol type |
| `absorption_optical_depth(layer)` (sum of `layer.absorbers`) | `CoreAbsorptionOpticalProperties(τ)` |
| scattering components within a layer | combined with vSmartMOM's `+` (weights `τ·ϖ`, identical to OceanOptics' scattering weights) |
| layers of an [`OceanColumnOptics`](@ref) | vertical composition by adding on the water quadrature grid, below the refracting interface — never `*`, which concatenates spectral points |
| `col.refractive_index` | interface index per spectral point (relative to air) |
| `bottom_albedo(col)` at `bottom_depth(col)` | Lambertian reflection operator (`m = 0`, intensity only) below the last layer; zero for `BlackBottom` |
| `layer.absorbers` (individual `τ`) | absorbed power per component from the layer-mean scalar irradiance ([`absorbed_power`](@ref), [`absorbed_par`](@ref)) |
| `layer.inelastic` | distributed sources (V3/V4): kernels of each `process`, excitation from the linked absorption components |
| `λ` `[nm]` | evaluate at `λ = 10⁷ ./ ν`; decreasing grids are accepted |

Rayleigh Greek values agree with `get_greek_rayleigh` exactly, and the
expansion of angular matrices was cross-checked against vSmartMOM's
`reconstruct_phase` on the Siewert (2000) aerosol
(`test/data/siewert_IIA_scattering_matrix.csv`). vSmartMOM's Stokes `Q`
(and `U`, `V`) have the opposite sign to the standard convention of
benchmark tables such as Chowdhary et al. (2020), as for its VLIDORT
comparisons: convert the truth, not the optical inputs.

## Adapter sketch (tested)

```julia
using OceanOptics, vSmartMOM
const S, RT = vSmartMOM.Scattering, vSmartMOM.CoreRT

to_vsmartmom(g::GreekCoefficients) = S.GreekCoefs(g.α, g.β, g.γ, g.δ, g.ϵ, g.ζ)

function layer_core_optics(layer, μ, m, pol, l_trunc)
    nλ   = n_wavelengths(layer)
    core = RT.CoreAbsorptionOpticalProperties(absorption_optical_depth(layer))
    scat = nothing
    for c in layer.scatterers
        is_wavelength_dependent(c) &&               # per-band expansion assumed here
            error("$(c.label): one expansion per band expected")
        g = greek_expansion(c, 1)
        if any(!iszero, g.γ) && any(!iszero, g.ϵ)   # see "δ-BGE" below
            tr = S.truncate_phase(S.δBGE(l_trunc),
                                  S.AerosolOptics(greek_coefs = to_vsmartmom(g),
                                                  ω̃ = 1.0, k = 1.0, fᵗ = 0.0))
            greek, fᵗ = tr.greek_coefs, tr.fᵗ
        else
            greek, fᵗ = to_vsmartmom(resize_moments(g, l_trunc)), 0.0
        end
        Z⁺⁺, Z⁻⁺ = S.compute_Z_moments(pol, μ, greek, m)
        comp = RT.CoreScatteringOpticalProperties((1 - fᵗ) .* c.τ, ones(nλ), Z⁺⁺, Z⁻⁺)
        scat = scat === nothing ? comp : scat + comp
    end
    return scat === nothing ? core : scat + core
end

ν   = collect(range(1e7 / 560, 1e7 / 540; length = 64))       # cm⁻¹
col = column_optics(column, 1e7 ./ ν; refractive_index = QuanFry1995(temperature = 20, salinity = 35),
                    bottom = BlackBottom(), ℓ_max = 400, polarized = true,
                    scalar_phase_approximation = :depolarizing)
cores = [layer_core_optics(l, μ_water, m, S.Stokes_IQUV(), l_trunc) for l in col]
```

`μ_water` stands for the water-side quadrature, including the refracted
nodes and the total-internal-reflection interval, which the interface
implementation defines.

## Precision and devices

The models are evaluated on the host; the solver chooses the element type
and the device of the results:

```julia
using Adapt
col = column_optics(FT, column, λ; refractive_index = index, ℓ_max = l_trunc)   # FT = Float32 or Float64
dev = adapt(architecture_array_type, col)      # e.g. CuArray; `Array` returns `col` without copying
```

Build the constituents in `Float64` and request `FT` here: tables and
constants are `Float64`, and each stored value is rounded once. `adapt`
moves every numeric array (`λ`, the `τ` of each component, the Greek
coefficients, backscatter fractions, `phase_index`, refractive index,
bottom albedo) and keeps the element type; labels, `phase_id` keys and
the component and layer lists stay host objects, so the adapter loops
over components on the host and hands their arrays to the kernels.
`absorption_optical_depth(layer)`, `scattering_optical_depth(layer)`,
`single_scattering_albedo(layer)` and `bottom_albedo(col)` work on the
adapted column and return device arrays.

Out of scope on a device: the quadrature that produces the Greek
coefficients (`greek_expansion`), the constituent, phase-function and
inelastic model structs (they hold or read host tables), and the absorbed
light and emission diagnostics. Compute these on the host, before `adapt`
(see "Precision and devices" in the [Consumer API](consumer_api.md)).

## Points for the solver side

* **δ-BGE with vanishing elements.** vSmartMOM's `truncate_phase(δBGE…)`
  weights its `γ` and `ϵ` fits by `1/f₁₂²` and `1/f₃₄²` and fails when
  those elements vanish identically — Rayleigh water (`ϵ ≡ 0`) and any
  `:depolarizing` component. Rayleigh water needs no truncation (`ℓ ≤ 2`);
  depolarizing components need a β-only truncation path, or δ-BGE must
  skip identically zero elements.
* **Per-band vs per-point phase matrices.** Components whose phase
  function is wavelength independent carry a single expansion; build
  their `Z` once per band, and reuse it across layers for equal
  `phase_id`. A wavelength-dependent component stores its distinct
  expansions plus `phase_index`, so `Z` is needed only per distinct entry
  (e.g. per band with `WavelengthDependentPhase(...; lookup = :nearest)`).
  `phase_id` holds one content key per entry of `greek`
  ([`phase_id`](@ref)`(c, iλ)` at a wavelength): equal keys mean equal
  coefficients, `ℓ_max` and element type, so `(phase_id, m)` is a valid
  cache key for `Z` before truncation.
* **Scalar expansions.** `GreekCoefs`-based code that sizes by
  `length(α)` needs [`zero_padded`](@ref)`(g)` for expansions computed with
  `polarized = false`.
* **Forward peaks.** Particle expansions (Fournier-Forand, hydrosol
  matrices) are not converged at practical `ℓ_max` — e.g. the Chowdhary
  AOS-IV hydrosols still have `β₁₀₀₀/2001 ≈ 3×10⁻³` — so truncation is
  required, with its own convergence study.
* **Absorbed light per layer.** Chlorophyll light-absorption profiles and
  fluorescence need, per layer and spectral point, the layer-mean scalar
  irradiance `(1/Δz) ∫ E₀ dz` (direct and diffuse, all directions),
  integrated within each layer — sublayers or an analytic layer integral,
  not interface values alone. With it, [`absorbed_power`](@ref),
  [`absorbed_par`](@ref) and [`inelastic_emission`](@ref) give the
  diagnostics and the reference fluorescence budget; energy closure
  (net irradiance below the surface = absorbed in the layers + absorbed by
  the floor) is the test.
* **Inelastic sources**: [`radiance_source_kernel`](@ref) is the energy
  kernel of every process, and [`conversion_units`](@ref) says whether the
  photon-energy factor λ′/λ enters it ([`PhotonConversion`](@ref) for
  fluorescence, [`EnergyConversion`](@ref) for the default
  [`WaterRaman`](@ref)); the kernel already contains that factor, so a
  solver must not apply λ′/λ again (see the [Consumer API](@ref) page). For
  fluorescence, [`excitation_optical_depth`](@ref)`(layer, j)` gives the
  excitation optical depth on the layer grid from the linked absorption
  components; [`excitation_absorption`](@ref) gives the same coefficient
  at any wavelength.
