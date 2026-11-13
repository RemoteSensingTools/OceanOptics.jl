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
    absorption_model = IOCCG2018(),
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

Stack the three constituents in a 5-m surface layer and compute its
optical properties at 10-nm resolution across the visible. The result
holds one absorption component per absorber and one scattering component
per scattering constituent (here water and phytoplankton; CDOM only
absorbs):

```@example ex
layer = OceanLayer(0.0, 5.0,
                   AbstractOceanConstituent{Float64}[water, cdom, phyto])
λ_grid = collect(400.0:10.0:700.0)
opt    = layer_optics(layer, λ_grid; ℓ_max = 32)
[(c.label, c.part) for c in opt.absorbers], [(c.label, c.polarization) for c in opt.scatterers]
```

Mixed quantities are derived from the components:

```@example ex
using Printf
τ, ϖ, B = optical_depth(opt), single_scattering_albedo(opt), backscatter_fraction(opt)
println(" λ [nm]   τ         ϖ         B")
for i in 1:5:length(λ_grid)
    @printf("%6.0f  %.5f   %.4f   %.4f\n", λ_grid[i], τ[i], ϖ[i], B[i])
end
```

Each component carries its own expansion; the mixture at a wavelength
is their scattering-weighted average, with `β₀ = 1`:

```@example ex
mixed = greek_expansion(opt, 1)
(n_moments(mixed), mixed.β[1], mixed.β[2] / 3)   # terms, β₀, asymmetry parameter at 400 nm
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

For a coupled solver, [`column_optics`](@ref) evaluates every layer on
one grid and adds the interface refractive index and the bottom; both
must be chosen explicitly:

```@example ex
col = column_optics(column, [450.0, 550.0, 650.0];
                    refractive_index = QuanFry1995(temperature = 20.0, salinity = 35.0),
                    bottom = BlackBottom(), ℓ_max = 32)
(n_layers(col), col.refractive_index)
```

## Chlorophyll absorption and fluorescence

To follow the light absorbed by chlorophyll, give the phytoplankton a
pigment partition. [`ProspectPigments`](@ref) splits its absorption into
chlorophyll and carotenoids with the leaf spectra of PROSPECT-D (the land
analogy); [`Bidigare1990Pigments`](@ref) uses ocean pigment groups. The
total absorption is unchanged:

```@example ex
pig = Phytoplankton(Chl = 0.5, pigments = ProspectPigments(carotenoid_ratio = 0.3))
a = absorption_by_part(pig, 490.0)
(absorption_parts(pig), a, sum(a) ≈ absorption(pig, 490.0))
```

Fluorescence is excited by chlorophyll absorption and by the share `η`
of carotenoid absorption passed to chlorophyll (`0`: chlorophyll only;
`1`: lossless antenna). The default emission is the chloroplast-level
spectrum of the leaf model Fluspect; the quantum yield is always given:

```@example ex
chl_only = chlorophyll_fluorescence(pig; quantum_yield = 0.02, carotenoid_transfer = 0.0)
antenna  = chlorophyll_fluorescence(pig; quantum_yield = 0.02, carotenoid_transfer = 0.5)
excitation_absorption(chl_only, 490.0), excitation_absorption(antenna, 490.0)
```

The emission density (photons per nm) has its maxima at 685 and 739 nm:

```@example ex
using Printf
println(" λ [nm]    f(λ) [1/nm]")
for λ in (660.0, 685.0, 700.0, 720.0, 739.0, 760.0, 800.0)
    @printf("%6.1f    %.5f\n", λ, emission(antenna, 440.0, λ))
end
```

[`radiance_source_kernel`](@ref) adds the photon-energy ratio `λ'/λ`
needed for an energy (radiance) source:

```@example ex
radiance_source_kernel(antenna, 440.0, 685.0) / inelastic_coefficient(antenna, 440.0, 685.0)
```

In a layer, the process is linked to the absorption components that
excite it, and the light each component absorbs follows from the
solver's layer-mean scalar irradiance. Here a uniform 1 W m⁻² nm⁻¹ stands
in for the solver's field:

```@example ex
lay = layer_optics(OceanLayer(0.0, 5.0, AbstractOceanConstituent{Float64}[water, pig];
                              inelastic = AbstractOceanInelasticProcess{Float64}[antenna]),
                   collect(400.0:5.0:750.0))
E₀  = ones(n_wavelengths(lay))
par = absorbed_par(lay, E₀)                  # µmol photons m⁻² s⁻¹ per absorption component
[(c.label, c.part, round(p; sigdigits = 4)) for (c, p) in zip(lay.absorbers, par)]
```

```@example ex
λ_em = collect(640.0:1.0:850.0)
P    = inelastic_emission(lay, 1, E₀, λ_em)  # W m⁻² nm⁻¹ emitted by the layer
λ_em[argmax(P)], maximum(P)
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
