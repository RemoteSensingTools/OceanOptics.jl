# vSmartMOM integration audit

This page records the result of a convention-by-convention audit of
OceanOptics.jl's optical-property outputs against the expectations of
[vSmartMOM.jl](https://github.com/RemoteSensingTools/vSmartMOM.jl) on
its `unified-vsmartmom` branch. It is the reference a future
solver-side adapter author should read before writing the
`OceanOptics → CoreScatteringOpticalProperties` bridge.

**Bottom line: all optical-property conventions match exactly.** Two
translation details are noted below — a unit conversion (nm ↔ μm) and
a shape adaptation (2D matrix → per-band 1D vector).

## Convention-by-convention

### 1. `CoreScatteringOpticalProperties` fields — **match**

From `src/CoreRT/types.jl:798-822`:

```julia
Base.@kwdef struct CoreScatteringOpticalProperties{FT,FT2,FT3} <: AbstractOpticalProperties
    "Absorption optical depth (scalar or wavelength dependent)"
    τ::FT
    "Single scattering albedo"
    ϖ::FT2
    "Z scattering matrix (forward)"
    Z⁺⁺::FT3
    "Z scattering matrix (backward)"
    Z⁻⁺::FT3
end
```

- `τ` is **layer-integrated** `(a + b) · Δz`, matching our
  [`OceanLayerOptics`](@ref)`.τ`.
- `ϖ = b / (a + b) = σ_s / σ_t`, matching our `.ϖ`.
- `Z⁺⁺`, `Z⁻⁺` are built *on the solver side* from the Greek
  coefficients and the solver's quadrature; OceanOptics.jl does not
  supply them.

### 2. `GreekCoefs` mapping to the `B` matrix — **match**

From `src/Scattering/types.jl:197-224` and
`src/Scattering/mie_helper_functions.jl:579`:

```julia
construct_B_matrix(mod::Stokes_IQUV, ...) =
    SMatrix{4,4}([β[l] γ[l]  0     0 ;
                  γ[l] α[l]  0     0 ;
                  0    0    ζ[l]  ϵ[l];
                  0    0   -ϵ[l]  δ[l]])
```

This is the Sanghavi (2014) Eq. 16 layout. OceanOptics.jl's
[`phase_matrix_moments`](@ref) returns exactly the same `(α, β, γ, δ,
ϵ, ζ)` named tuple with the same mapping, and the
[`OceanLayerOptics`](@ref) `.α, .β, .γ, .δ, .ϵ, .ζ` fields document the
same positions.

### 3. Vector-index convention — **match**

vSmartMOM stores each Greek coefficient as `Array{FT, 1}` indexed
`[ℓ+1]`, so `β[1]` is `ℓ = 0`. OceanOptics.jl uses the same: `β_ℓ` at
Julia index `ℓ + 1` (row 1 of the `(ℓ_max+1, nλ)` matrix is the
`ℓ = 0` row). Extracting a single wavelength with `opt.β[:, iλ]`
yields exactly vSmartMOM's expected 1D layout.

### 4. Rayleigh closed-form Greek values — **match**

OceanOptics.jl's [`phase_matrix_moments`](@ref) for
[`RayleighWaterPhase`](@ref) is ported verbatim from vSmartMOM's
`get_greek_rayleigh` at
`src/Scattering/mie_helper_functions.jl:454-468`:

```julia
dpl_p = (1 - ρ)  / (1 + ρ/2)
dpl_r = (1 - 2ρ) / (1 - ρ)

β[0] = 1            β[2] = ½·dpl_p
α[2] = 3·dpl_p      γ[2] = √(3/2)·dpl_p
δ[1] = (3/2)·dpl_p·dpl_r
ϵ[ℓ] = ζ[ℓ] = 0
```

### 5. Pipeline composition operators — **match**

vSmartMOM composes `CoreScatteringOpticalProperties` with
`Base.:+` for per-wavelength / per-material spectral mixing and
`Base.:*` for vertical stacking across layers
(`src/CoreRT/LayerOpticalProperties/compEffectiveLayerProperties.jl:20-45`).
OceanOptics.jl's [`IOP`](@ref) uses the same `+` semantics for
constituent mixing, and its layers already arrive as independent
`OceanLayerOptics` values that stack naturally.

## Translation details (not mismatches)

### 6. Wavelength units — OceanOptics uses **nm**, vSmartMOM uses **μm**

This is the ocean-optics vs atmospheric-RT unit-system split; both are
standard in their own communities. The adapter must divide wavelength
by 1000 when handing values into vSmartMOM:

```julia
# At the integration boundary:
core = CoreScatteringOpticalProperties(
    τ   = opt.τ[iλ],
    ϖ   = opt.ϖ[iλ],
    Z⁺⁺ = …,
    Z⁻⁺ = …,
)
# …or convert the OceanOptics wavelength grid up front:
λ_μm = opt.λ ./ 1000
```

OceanOptics.jl deliberately keeps nm internally because every ocean-
optics data table (Pope-Fry, Smith-Baker, Bricaud, Mason-Cone-Fry,
Pegau, Petzold) is tabulated in nm.

### 7. Greek-coefficient shape: 2D `(ℓ+1, nλ)` vs per-band 1D vector

vSmartMOM's `GreekCoefs` stores `Array{FT, 1}` — one coefficient
vector per band (wavelength variation handled by having one
`GreekCoefs` per band). OceanOptics.jl stores `(ℓ+1, nλ)` matrices on
`OceanLayerOptics` for compactness. The adapter builds per-wavelength
`GreekCoefs` by column extraction:

```julia
gc = GreekCoefs(
    α = opt.α[:, iλ],
    β = opt.β[:, iλ],
    γ = opt.γ[:, iλ],
    δ = opt.δ[:, iλ],
    ϵ = opt.ϵ[:, iλ],
    ζ = opt.ζ[:, iλ],
)
```

## Still open in vSmartMOM (not OceanOptics scope)

### 8. `OceanSurface <: AbstractSurfaceType` does not yet exist

`src/CoreRT/types.jl` currently defines `CoxMunkSurface` (wind-
roughened BRDF) but no ocean-scattering-layer surface type. The DESIGN
memo (§10, Phase 3) places this work on the vSmartMOM side, not here:
the `OceanSurface` would consume `OceanLayerOptics` bundles and hand
them to the existing CoreRT pipeline alongside the atmospheric layers
via `Base.:*`. Implementing it is solver-side work outside this
repository's scope.

## Adapter skeleton

For reference, the concrete adapter a solver-side PR would provide
(pseudocode; lives in vSmartMOM, not here):

```julia
function core_scattering_from_ocean(opt::OceanLayerOptics,
                                    iλ::Int,
                                    quad_pts,
                                    pol_type)
    # 1. per-wavelength Greek coefficients
    gc = GreekCoefs(α = opt.α[:, iλ], β = opt.β[:, iλ],
                    γ = opt.γ[:, iλ], δ = opt.δ[:, iλ],
                    ϵ = opt.ϵ[:, iλ], ζ = opt.ζ[:, iλ])
    # 2. forward / backward Z matrices from Greek + quadrature
    Z⁺⁺, Z⁻⁺ = compute_Z_matrices(gc, quad_pts, pol_type)
    # 3. wavelength unit convention is vSmartMOM's: convert if needed
    return CoreScatteringOpticalProperties(
        τ   = opt.τ[iλ],
        ϖ   = opt.ϖ[iλ],
        Z⁺⁺ = Z⁺⁺,
        Z⁻⁺ = Z⁻⁺,
    )
end
```

The inelastic adapter is analogous: iterate
`layer.fluorophores`, evaluate
[`excitation_absorption`](@ref), [`emission`](@ref), and
[`inelastic_coefficient`](@ref) against the pass-1 irradiance, and
build the `OceanRS <: AbstractRamanType` the solver expects.

## When this audit was last verified

- vSmartMOM branch `unified-vsmartmom`, as of April 2026.
- OceanOptics.jl commit `c7359da` (Phase 4 polarized-Rayleigh +
  Documenter).

If vSmartMOM changes its `CoreScatteringOpticalProperties` signature,
its `GreekCoefs` layout, or the Rayleigh closed-form coefficients in
`get_greek_rayleigh`, this page must be re-verified.
