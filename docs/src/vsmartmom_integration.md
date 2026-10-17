# vSmartMOM integration

How [vSmartMOM.jl](https://github.com/RemoteSensingTools/vSmartMOM.jl)
consumes OceanOptics. The adapter below was run against the vSmartMOM
release candidate `d441061` (October 2026) in a joint environment; it is
documentation, not part of either package.

## Compatibility

* OceanOptics' dependency bounds include the versions vSmartMOM resolves
  to (`Interpolations` 0.16, `FastGaussQuadrature` 0.5 through its
  CanopyOptics dependency); the OceanOptics test suite passes with them.
* No exported names collide, so `using vSmartMOM, OceanOptics` works
  without qualification. (vSmartMOM's own `greek_coefficients`,
  `phase_function` and `scattering_matrix` refer to its analytic phase
  functions; OceanOptics uses [`greek_expansion`](@ref),
  [`scattering_phase`](@ref) and [`evaluate_scattering_matrix`](@ref).)
* OceanOptics is not registered yet: register and tag it before a
  vSmartMOM release declares it in `[weakdeps]`/`[compat]`.

## Mapping

| OceanOptics | vSmartMOM |
|---|---|
| [`GreekCoefficients`](@ref) `(α, β, γ, δ, ϵ, ζ)` | `Scattering.GreekCoefs(α, β, γ, δ, ϵ, ζ)` — same fields, order and signs |
| [`ScatteringComponent`](@ref) | one `CoreScatteringOpticalProperties` per component, `ϖ = 1`, truncated like an aerosol type |
| `absorption_optical_depth(layer)` | `CoreAbsorptionOpticalProperties(τ)` |
| components within a layer | combined with vSmartMOM's `+` (weights `τ·ϖ`, identical to OceanOptics' scattering weights) |
| layers of an [`OceanColumnOptics`](@ref) | vertical composition by adding on the water quadrature grid, below the refracting interface — never `*`, which concatenates spectral points |
| `col.refractive_index` | interface index per spectral point (relative to air) |
| `bottom_albedo(col)` | lower boundary of the water grid |
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
    for c in layer.components
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
* **Forward peaks.** Particle expansions (Fournier-Forand, hydrosol
  matrices) are not converged at practical `ℓ_max` — e.g. the Chowdhary
  AOS-IV hydrosols still have `β₁₀₀₀/2001 ≈ 3×10⁻³` — so truncation is
  required, with its own convergence study.
* **Inelastic sources** follow the photon contract on the
  [Consumer API](@ref) page ([`radiance_source_kernel`](@ref)).
