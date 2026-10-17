# vSmartMOM integration audit

This page records the result of a convention-by-convention audit of
OceanOptics.jl's optical-property outputs against the expectations of
[vSmartMOM.jl](https://github.com/RemoteSensingTools/vSmartMOM.jl) on
its `unified-vsmartmom` branch. It is the reference a future
solver-side adapter author should read before writing the
`OceanOptics → CoreScatteringOpticalProperties` bridge.

**Bottom line:** the per-layer optical quantities (`τ`, `ϖ`) and the
Greek-coefficient layout, indexing and Rayleigh values match. Three
translation details are noted below — spectral units (OceanOptics nm,
vSmartMOM core RT cm⁻¹), a shape adaptation (2D matrix → per-wavelength
1D vector), and vertical composition, which vSmartMOM performs by
adding, not by an operator on optical-property objects.

!!! note "Re-checked 2026-10-16"
    Sections 1, 4, 5 and 6 were re-checked against the vSmartMOM release
    candidate `d441061` (October 2026). The original April 2026 audit
    described `Base.:*` as vertical stacking and vSmartMOM's spectral unit
    as μm (true only for its Mie inputs); both are corrected below.

## Convention-by-convention

### 1. `CoreScatteringOpticalProperties` fields — **match**

`CoreScatteringOpticalProperties` (vSmartMOM `src/CoreRT/types.jl`)
carries `τ`, `ϖ`, the forward/backward phase operators `Z⁺⁺`, `Z⁻⁺`, and
two optional direct-solar columns `Z₀⁺`, `Z₀⁻` used only by the
external-solar path:

- `τ` is **layer-integrated** `(a + b) · Δz`, matching our
  [`OceanLayerOptics`](@ref)`.τ`.
- `ϖ = b / (a + b) = σ_s / σ_t`, matching our `.ϖ`.
- `Z⁺⁺`, `Z⁻⁺` are built *on the solver side* from the Greek
  coefficients and the solver's quadrature; OceanOptics.jl does not
  supply them.

### 2. `GreekCoefs` mapping to the `B` matrix — **match**

From vSmartMOM `src/Scattering/types.jl` (`GreekCoefs`) and
`src/Scattering/mie_helper_functions.jl` (`construct_B_matrix`):

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
`get_greek_rayleigh` (`src/Scattering/mie_helper_functions.jl`), including
the Hovenier-positive `γ₂`:

```julia
dpl_p = (1 - ρ)  / (1 + ρ/2)
dpl_r = (1 - 2ρ) / (1 - ρ)

β[0] = 1            β[2] = ½·dpl_p
α[2] = 3·dpl_p      γ[2] = √(3/2)·dpl_p
δ[1] = (3/2)·dpl_p·dpl_r
ϵ[ℓ] = ζ[ℓ] = 0
```

### 5. Composition operators — **different roles**

vSmartMOM's `Base.:+` on `CoreScatteringOpticalProperties` mixes
components *within one layer*: optical depths add, and `ϖ` and the phase
operators are weighted by each component's `τ·ϖ`. That is the same
scattering-weighted mixing [`layer_optics`](@ref) applies to constituents
(and [`IOP`](@ref) `+` applies to coefficients).

vSmartMOM's `Base.:*` does **not** stack layers vertically: it
concatenates spectral points (`[x.τ; y.τ]`, phase operators along the
third dimension). Vertical composition happens in the RT kernel's adding
(interaction) step. An ocean column therefore enters the solver as a
sequence of layers for adding — separately from the atmosphere, below a
refracting air–water interface — never through `*`.

## Translation details (not mismatches)

### 6. Spectral units — OceanOptics uses **nm**, vSmartMOM core RT uses **cm⁻¹**

vSmartMOM's RT spectral grids are wavenumbers in cm⁻¹ (its Mie code
takes wavelengths in μm). The adapter evaluates OceanOptics on the
wavelengths `λ = 10⁷ / ν` nm of the solver's wavenumber grid:

```julia
λ_nm = 1e7 ./ ν_cm⁻¹          # decreasing when ν is increasing
opt  = layer_optics(layer, λ_nm; ℓ_max)
```

[`OceanLayerOptics`](@ref) accepts strictly decreasing wavelength grids,
so the bundle can be built directly in the solver's order. Per-nm
spectral densities (`J_source`, inelastic redistribution) additionally
need the Jacobian `|dλ/dν| = 10⁷/ν²` and the solver's bin weights.

OceanOptics.jl keeps nm internally because every ocean-optics data table
(Pope-Fry, Smith-Baker, Bricaud, Mason-Cone-Fry, Pegau, Petzold) is
tabulated in nm.

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
an ocean boundary would solve the water column (from `OceanLayerOptics`
bundles) and its bottom internally, couple it to the atmosphere through
the four reflection/transmission blocks of a refracting air–water
interface, and present the result to the atmospheric adding step.
Implementing it is solver-side work outside this repository's scope.

## Adapter skeleton

Illustrative pseudocode for the per-layer part of a solver-side adapter
(lives in vSmartMOM, not here). Function names on the vSmartMOM side are
placeholders; check the current API — e.g. its phase operators are built
per Fourier moment by `Scattering.compute_Z_moments`. A polarized solver
must check [`is_polarized`](@ref)`(opt)` before reading `α…ζ`.

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

- vSmartMOM branch `unified-vsmartmom`, as of April 2026; sections 1,
  4, 5 and 6 re-checked against release candidate `d441061`, October 2026.
- OceanOptics.jl commit `c7359da` (Phase 4 polarized-Rayleigh +
  Documenter).

If vSmartMOM changes its `CoreScatteringOpticalProperties` signature,
its `GreekCoefs` layout, or the Rayleigh closed-form coefficients in
`get_greek_rayleigh`, this page must be re-verified.
