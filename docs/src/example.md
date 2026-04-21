# End-to-end example

A realistic Case-1 water workflow: build the biogeochemistry, compute
inherent optical properties over the visible, inspect the inelastic
source kernels, and differentiate the elastic-path output w.r.t. the
chlorophyll concentration.

All code blocks below are executed at documentation-build time by
Documenter's `@example` machinery, so the numerical outputs in the
rendered page are the *actual* package outputs.

## Biogeochemistry

Open-ocean oligotrophic water at `Chl = 0.5 mg/m³`, surface temperature
`T = 20 °C`, practical salinity `S = 35`.

```@example ex
using OceanOptics

water = PureWater{Float64}(
    absorption_model = CombinedPopeFryMason(),
    scattering_model = ZhangHu2009(),
    temperature      = 20.0,
    salinity         = 35.0,
)
cdom  = CDOM{Float64}(a_ref = 0.03, λ_ref = 440.0, slope = 0.014)
phyto = Phytoplankton{Float64, Bricaud1998}(Chl = 0.5)
nothing # hide
```

Each constituent exposes the same `absorption` / `scattering` /
`backscattering` traits. Spot-checks at the chlorophyll blue peak:

```@example ex
λ = 440.0
a_water = absorption(water, λ)
b_water = scattering(water, λ)
a_cdom  = absorption(cdom, λ)
a_phyto = absorption(phyto, λ)
b_phyto = scattering(phyto, λ)
(a_water, b_water, a_cdom, a_phyto, b_phyto)
```

## Building a layer and computing the optical bundle

Stack the three constituents in a 5-m surface layer and compute the
elastic optical bundle at 10-nm resolution across the visible:

```@example ex
layer = OceanLayer(0.0, 5.0,
                   AbstractOceanConstituent{Float64}[water, cdom, phyto])
λ_grid = collect(400.0:10.0:700.0)
opt    = layer_optics(layer, λ_grid; ℓ_max = 32)
opt
```

Inspect the single-scattering albedo and the backscatter fraction
across the spectrum:

```@example ex
using Printf
println(" λ [nm]   τ         ϖ         B")
for i in eachindex(opt.λ)
    @printf("%6.0f  %.5f   %.4f   %.4f\n",
            opt.λ[i], opt.τ[i], opt.ϖ[i], opt.B[i])
end
```

The Legendre-moment matrix `opt.β` has size `(ℓ_max+1, n_λ)` with
`β_0 ≡ 1` by construction:

```@example ex
size(opt.β), opt.β[1, 1], opt.β[2, 1], opt.β[3, 1]
```

## Multi-layer column via the Fell (1997) grid

The Fell (1997) §5.2.1 grid (1 m in 0–10 m, 2 m in 10–20 m, 5 m
below) is available as a one-liner:

```@example ex
column = fell_column(AbstractOceanConstituent{Float64}[water, cdom, phyto];
                     depth = 60.0)
(n_layers(column), thickness(column), bottom_depth(column))
```

```@example ex
layer_thicknesses(column)
```

## Inelastic kernels

Chlorophyll-a fluorescence is built from the phytoplankton absorber
plus a Gaussian emission band at 685 nm. The package returns kernels
a downstream two-pass RT solver convolves against the pass-1
irradiance; we do not solve RT here.

```@example ex
sif = chlorophyll_fluorescence(phyto)
sif
```

```@example ex
excitation_absorption(sif, 440.0)           # φ · a_φ(440)
```

The emission density peaks at 685 nm with FWHM ≈ 25 nm:

```@example ex
using Printf
println(" λ [nm]    f(440, λ) [1/nm]")
for λ in (640.0, 670.0, 680.0, 685.0, 690.0, 700.0, 730.0)
    @printf("%6.1f    %.5f\n", λ, emission(sif, 440.0, λ))
end
```

Water Raman at 440-nm excitation peaks around 517 nm (the `3400 cm⁻¹`
Stokes shift of the O–H stretch):

```@example ex
ram = WaterRaman{Float64}()
(raman_peak_wavelength(ram, 440.0),
 excitation_absorption(ram, 440.0),
 emission(ram, 440.0, 517.0))
```

## Forward-mode automatic differentiation

Every forward function along the elastic path is AD-transparent via
`ForwardDiff.Dual` inputs on the *wavelength* and most model
parameters. Spectral sensitivity of the total absorption coefficient
at 440 nm:

```@example ex
using ForwardDiff
∂a_∂λ = ForwardDiff.derivative(
    λ -> absorption(water, λ) + absorption(cdom, λ) + absorption(phyto, λ),
    440.0,
)
```

Struct parameters themselves are typed `FT<:AbstractFloat`, so
differentiating with respect to e.g. the CDOM slope goes through the
analytic formula rather than the struct:

```@example ex
a_ref, λ_ref = 0.03, 440.0
g(S) = sum(a_ref * exp(-S * (λ - λ_ref)) for λ in (500.0, 600.0, 700.0))
ForwardDiff.derivative(g, 0.014)
```

See `test/test_ad.jl` for a wider set of worked AD examples, including
CDOM, Pure Water, Morel / Zhang-Hu, Fournier-Forand and
Henyey-Greenstein sensitivities.
