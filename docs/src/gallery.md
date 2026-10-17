# Gallery

Visual sanity checks of OceanOptics.jl outputs against the canonical
published figures they should reproduce. All charts below are generated
at documentation-build time from the actual package code by the
`@example` block shown above each one. They are interactive: hover for
values, click a legend entry to hide a series, and open *Data table* for
the plotted numbers.

Each figure cites the published work it is intended to match. If you
spot a qualitative disagreement with one of those references, it is a
bug.

```@setup gallery
using OceanOptics
using Main.DocCharts
```

## 1. Pure-water absorption, `a_w(λ)`

Comparison target: Pope & Fry (1997) Fig. 4 and the consolidated
pure-water absorption curve on the Ocean Optics Web Book
["Absorption by pure water"](https://www.oceanopticsbook.info/view/absorption/absorption-by-pure-sea-water)
page. Expect a pronounced minimum around 420 nm (≈ 0.005 m⁻¹), a
steep rise into the NIR, and the characteristic water-vapor / O–H
combination bands producing shoulders near 605 and 740 nm. Each table
is drawn only over its own range (Pope & Fry: 380–727.5 nm); outside it a
table has no data and OceanOptics throws. See also the
[Comparison with HydroLight IOPs](@ref).

```@example gallery
λ   = collect(350.0:2.0:800.0)
λpf = filter(x -> 380 ≤ x ≤ 727.5, λ)

w_popefry  = PureWater{Float64}(absorption_model = PopeFry1997())
w_smithb   = PureWater{Float64}(absorption_model = SmithBaker1981())
w_ioccg    = PureWater{Float64}(absorption_model = IOCCG2018())

a_pf = [absorption(w_popefry,  λi) for λi in λpf]
a_sb = [absorption(w_smithb,   λi) for λi in λ]
a_io = [absorption(w_ioccg,    λi) for λi in λ]

linechart([Series("Pope & Fry 1997", λpf, a_pf),
           Series("Smith & Baker 1981", λ, a_sb; dash = true),
           Series("IOCCG 2018 + Kou 1993 (default)", λ, a_io; dot = true)];
          title = "Pure-water absorption", xlabel = "wavelength λ [nm]",
          ylabel = "a_w(λ) [m⁻¹]", xname = "λ [nm]", ylog = true)
```

## 2. Pure-water scattering — Morel (1974) vs Zhang-Hu (2009)

Comparison target: Zhang, Hu & He (2009) Fig. 3 (seawater at S =
38.4 ‰, T = 20 °C, 350–700 nm), who report agreement with Morel's
measurements to about 1 %. Morel's (1974) power-law fit for seawater runs
somewhat above the Zhang-Hu model (7–16 % at 18 °C, S = 35; see the
[Comparison with HydroLight IOPs](@ref)), and freshwater (S = 0)
scatters roughly 30 % less than seawater.

```@example gallery
λ = collect(350.0:5.0:700.0)

w_morel   = PureWater{Float64}(scattering_model = Morel1974())
w_zhh_38  = PureWater{Float64}(scattering_model = ZhangHu2009(),
                               salinity = 38.4, temperature = 20.0)
w_zhh_00  = PureWater{Float64}(scattering_model = ZhangHu2009(),
                               salinity = 0.0,  temperature = 20.0)

b_mo = [scattering(w_morel,  λi) for λi in λ]
b_zh = [scattering(w_zhh_38, λi) for λi in λ]
b_zp = [scattering(w_zhh_00, λi) for λi in λ]

linechart([Series("Morel 1974", λ, b_mo),
           Series("Zhang-Hu 2009, S = 38.4", λ, b_zh),
           Series("Zhang-Hu 2009, S = 0 (pure)", λ, b_zp; dash = true)];
          title = "Pure (sea)water scattering", xlabel = "wavelength λ [nm]",
          ylabel = "b_w(λ) [m⁻¹]", xname = "λ [nm]", ylog = true)
```

## 3. Case-1 water IOPs at varying chlorophyll

Comparison target: Morel & Maritorena (2001) Fig. 8 (total absorption
and scattering of Case-1 water versus chlorophyll) and IOCCG Report 5
§2.4. At 440 nm, higher Chl produces both higher total absorption
(more `a_φ`) and higher particulate scattering `b_p`, with the
classical greening of `a(λ)` away from the pure-water blue minimum.

```@example gallery
λ   = collect(400.0:2.0:700.0)
Chls = (0.03, 0.1, 0.3, 1.0, 3.0)

water = PureWater{Float64}()           # IOCCG2018 + Zhang-Hu
cdom  = CDOM{Float64}(a_ref = 0.03)

function total_a(phyto)
    return [absorption(water, λi) + absorption(cdom, λi) + absorption(phyto, λi)
            for λi in λ]
end
function total_b(phyto)
    return [scattering(water, λi) + scattering(phyto, λi) for λi in λ]
end

phytos = [Phytoplankton{Float64, Bricaud1998}(Chl = Chl) for Chl in Chls]
linechart([Series("Chl = $Chl mg/m³", λ, total_a(p)) for (Chl, p) in zip(Chls, phytos)];
          title = "Total absorption a = a_w + a_g + a_φ", xlabel = "wavelength λ [nm]",
          ylabel = "a(λ) [m⁻¹]", xname = "λ [nm]", ylog = true, ordinal = true)
```

```@example gallery
linechart([Series("Chl = $Chl mg/m³", λ, total_b(p)) for (Chl, p) in zip(Chls, phytos)];
          title = "Total scattering b = b_w + b_p", xlabel = "wavelength λ [nm]",
          ylabel = "b(λ) [m⁻¹]", xname = "λ [nm]", ylog = true, ordinal = true)
```

## 4. Layer optics `τ, ϖ, B` for a 5-m Case-1 layer

Comparison target: typical HydroLight / Mobley "clear ocean" reference
output, e.g. Mobley (1994) Table 4.2 and the Ocean Optics Web Book
Case-1 water IOP reference pages. At `Chl = 0.5 mg/m³` the
single-scattering albedo `ϖ(λ)` should peak near 500 nm (minimum
water absorption combined with non-negligible particulate scattering)
and fall toward the red as water absorption takes over.

```@example gallery
λ     = collect(400.0:5.0:700.0)
water = PureWater{Float64}()
cdom  = CDOM{Float64}(a_ref = 0.03)
phyto = Phytoplankton{Float64, Bricaud1998}(Chl = 0.5)
layer = OceanLayer(0.0, 5.0,
                   AbstractOceanConstituent{Float64}[water, cdom, phyto])
opt   = layer_optics(layer, λ; ℓ_max = 32)

linechart([Series("τ", opt.λ, optical_depth(opt))]; title = "Layer optical depth",
          xlabel = "λ [nm]", ylabel = "τ(λ)", xname = "λ [nm]", height = 300)
```

```@example gallery
linechart([Series("ϖ", opt.λ, single_scattering_albedo(opt))]; title = "Single-scattering albedo",
          xlabel = "λ [nm]", ylabel = "ϖ(λ)", xname = "λ [nm]", height = 300)
```

```@example gallery
linechart([Series("B", opt.λ, backscatter_fraction(opt))]; title = "Backscatter fraction",
          xlabel = "λ [nm]", ylabel = "B(λ)", xname = "λ [nm]", height = 300)
```

## 5. Phase-function comparison — Fournier-Forand vs Henyey-Greenstein vs Rayleigh-water

Comparison target: Mobley, Sundman & Boss (2002) Fig. 1 (Fournier-
Forand at several backscatter fractions vs. Petzold) and Hansen &
Travis (1974) Fig. 7 for the Rayleigh-with-depolarization shape.
Fournier-Forand should be strongly forward-peaked and fall off faster
than Henyey-Greenstein at intermediate angles; pure Rayleigh is
symmetric forward / back with a p(90°) minimum.

```@example gallery
μ   = collect(-1.0:0.01:1.0)

pf_ff   = FournierForandPhase{Float64}(backscatter_fraction = 0.018)
pf_hg95 = HenyeyGreensteinPhase{Float64}(g = 0.95)
pf_hg50 = HenyeyGreensteinPhase{Float64}(g = 0.50)
pf_ray  = RayleighWaterPhase{Float64}(depolarization = 0.039)

ff   = [phase_function_value(pf_ff,   μi) for μi in μ]
hg95 = [phase_function_value(pf_hg95, μi) for μi in μ]
hg50 = [phase_function_value(pf_hg50, μi) for μi in μ]
ray  = [phase_function_value(pf_ray,  μi) for μi in μ]

θ_deg = rad2deg.(acos.(μ))
linechart([Series("Fournier-Forand, B = 0.018", θ_deg, ff),
           Series("Henyey-Greenstein, g = 0.95", θ_deg, hg95),
           Series("Henyey-Greenstein, g = 0.50", θ_deg, hg50),
           Series("Rayleigh-water, ρ = 0.039", θ_deg, ray)];
          title = "Phase functions for ocean scattering",
          xlabel = "scattering angle θ [deg]", ylabel = "p(θ)", xname = "θ [deg]",
          ylog = true, digits = 4)
```

## 6. Chlorophyll fluorescence emission spectrum

Comparison target: Fell (1997) Fig. 5.5 and Gordon (1979) — Gaussian
emission band centred at 685 nm with σ = 10.6 nm (FWHM ≈ 25 nm). The
band is excitation-independent in the 400–700 nm PAR window, so the
shape does not shift with `λ'`.

```@example gallery
phyto = Phytoplankton{Float64, Bricaud1998}(Chl = 0.5)
sif   = chlorophyll_fluorescence(phyto)

λ_emit = collect(640.0:1.0:740.0)
f_exc  = [(λ′, [emission(sif, λ′, λ) for λ in λ_emit]) for λ′ in (400.0, 500.0, 600.0)]

# The three curves coincide (excitation-independent band); dashes keep
# them distinguishable where they overlap.
styles = ((dash = false, dot = false), (dash = true, dot = false), (dash = false, dot = true))
linechart([[Series("λ' = $(Int(λ′)) nm", λ_emit, y; st...)
            for ((λ′, y), st) in zip(f_exc, styles)];
           Series("peak = 685 nm", [685.0, 685.0], [0.0, maximum(last(f_exc[1]))];
                  reference = true, dot = true)];
          title = "Chlorophyll fluorescence emission (Fell 1997, Fig. 5.5)",
          xlabel = "emission wavelength λ [nm]", ylabel = "f(λ', λ) [nm⁻¹]",
          xname = "λ [nm]")
```

## 7. Water Raman emission vs excitation wavelength

Comparison target: Haltrin & Kattawar (1991) and Bartlett et al.
(1998) — a Stokes-shifted Gaussian emission in wavenumber (fixed
`Δν ≈ 3400 cm⁻¹`) that maps onto a variable-width band in wavelength.
Excitation at 440 nm produces an emission peak at ≈ 517 nm; at 550 nm
the peak is at ≈ 673 nm.

```@example gallery
ram = WaterRaman{Float64}()

λ_emit = collect(430.0:1.0:720.0)
λ_excs = (400.0, 440.0, 500.0, 550.0)

linechart([Series("λ' = $(Int(round(λ′))) nm (peak $(round(Int, raman_peak_wavelength(ram, λ′))) nm)",
                  λ_emit, [emission(ram, λ′, λ) for λ in λ_emit]) for λ′ in λ_excs];
          title = "Water Raman emission — Stokes shift = 3400 cm⁻¹",
          xlabel = "emission wavelength λ [nm]", ylabel = "f_R(λ', λ) [nm⁻¹]",
          xname = "λ [nm]", ordinal = true)
```
