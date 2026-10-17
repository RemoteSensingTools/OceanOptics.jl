# Comparison with HydroLight IOPs

Loisel et al. (2023) published the inherent optical properties used as
input for 29 880 HydroLight v5.0 simulations: 3320 combinations of
phytoplankton, non-algal particle and CDOM properties plus pure seawater,
at 350–750 nm in 5 nm steps
([ESSD article](https://doi.org/10.5194/essd-15-3711-2023),
[Dryad dataset](https://doi.org/10.6076/D1630T), CC0). This page compares
OceanOptics.jl with those inputs.

The file `data/reference/loisel2023_hydrolight_iops.csv` holds the extract
used here and in the test suite: the pure-seawater coefficients (the same
in every scenario; derived as total minus non-water) and three scenarios
at the 5th, 50th and 95th percentile of total `a(440)`. Its header records
the source file and checksum. Database values are single precision, with
totals rounded to about five significant figures.

The paper's pure-seawater assumptions:

| Quantity | HydroLight input (Loisel et al. 2023, §2) |
|---|---|
| `a_w(λ)` | IOCCG Protocol Series (2018) composite: Morel et al. (2007) 340–415 nm, Pope & Fry (1997) 420–725 nm, Kou et al. (1993) 730–750 nm |
| `b_w(λ)` | Zhang et al. (2009) at 18 °C and salinity 35 |
| particles | Fournier-Forand phase functions, `b_b/b = 0.01` (phytoplankton), `0.018` (non-algal) |

```@setup hl
using OceanOptics, Printf
using Main.DocCharts

ref = OceanOptics._read_csv_columns(
    OceanOptics.data_path("reference/loisel2023_hydrolight_iops.csv"))
λ   = ref["lambda_nm"]
water(; a = CombinedPopeFryMason(), b = ZhangHu2009(), T = 18.0, S = 35.0) =
    PureWater{Float64}(absorption_model = a, scattering_model = b,
                       temperature = T, salinity = S)
```

## Pure-seawater scattering

OceanOptics' port of the Zhang, Hu & He (2009) model, evaluated at the
same temperature and salinity, reproduces the HydroLight `b_w` to the
precision of the stored data over the whole range. Morel's (1974) power
law is 7–16 % higher.

```@example hl
zhh   = scattering(water(), λ) ./ ref["bw"]
morel = scattering(water(b = Morel1974()), λ) ./ ref["bw"]
@printf("Zhang-Hu 2009 / HydroLight: %.6f – %.6f\n", extrema(zhh)...)
linechart([Series("Zhang-Hu 2009, 18 °C, S = 35", λ, zhh),
           Series("Morel 1974", λ, morel; dash = true),
           Series("HydroLight", λ, ones(length(λ)); reference = true, dot = true)];
          xlabel = "wavelength λ [nm]", ylabel = "b_w(OceanOptics) / b_w(HydroLight)",
          xname = "λ [nm]", ymin = 0.95, ymax = 1.2)
```

## Pure-water absorption

Between 420 and 725 nm the HydroLight values are the Pope & Fry (1997)
table, which OceanOptics reproduces exactly. Outside that range the two
use different sources: the default `CombinedPopeFryMason` model falls back
to Smith & Baker (1981) below 380 nm and above 700 nm, because the
Mason, Cone & Fry (2016) UV table is not bundled yet. Smith & Baker's UV
values are several times higher than those of Morel et al. (2007).
`PopeFry1997` alone is defined for 380–727.5 nm and is shown only there:
outside its table it would return the end values.

```@example hl
in_range(lo, hi) = lo .≤ λ .≤ hi
λpf  = λ[in_range(380, 727.5)]
a_hl = ref["aw"]
a_co = absorption(water(T = 20.0), λ)                       # default model, no T correction
a_pf = absorption(water(a = PopeFry1997(), T = 20.0), λpf)
a_sb = absorption(water(a = SmithBaker1981(), T = 20.0), λ)
linechart([Series("HydroLight (IOCCG 2018)", λ, a_hl; points = true, reference = true),
           Series("CombinedPopeFryMason (default)", λ, a_co),
           Series("PopeFry1997", λpf, a_pf; dash = true),
           Series("SmithBaker1981", λ, a_sb; dot = true)];
          xlabel = "wavelength λ [nm]", ylabel = "a_w(λ) [m⁻¹]", xname = "λ [nm]",
          ylog = true)
```

The ratio makes the differences explicit:

```@example hl
linechart([Series("CombinedPopeFryMason (default)", λ, a_co ./ a_hl),
           Series("PopeFry1997", λpf, a_pf ./ a_hl[in_range(380, 727.5)]; dash = true),
           Series("SmithBaker1981", λ, a_sb ./ a_hl; dot = true),
           Series("HydroLight", λ, ones(length(λ)); reference = true, dot = true)];
          xlabel = "wavelength λ [nm]", ylabel = "a_w(OceanOptics) / a_w(HydroLight)",
          xname = "λ [nm]", ylog = true)
```

```@example hl
ranges = [(350, 375, "Morel et al. 2007", "Smith & Baker 1981"),
          (380, 415, "Morel et al. 2007", "Pope & Fry 1997"),
          (420, 700, "Pope & Fry 1997",   "Pope & Fry 1997"),
          (705, 725, "Pope & Fry 1997",   "Smith & Baker 1981"),
          (730, 750, "Kou et al. 1993",   "Smith & Baker 1981")]
@printf("%-9s  %-18s  %-20s  %s\n", "λ [nm]", "HydroLight", "OceanOptics default", "ratio")
for (lo, hi, hl, oo) in ranges
    r = (a_co ./ a_hl)[in_range(lo, hi)]
    @printf("%3d–%3d    %-18s  %-20s  %.3f – %.3f\n", lo, hi, hl, oo, extrema(r)...)
end
```

OceanOptics also applies the Pegau et al. (1997) temperature correction
between 600 and 800 nm, referenced to 20 °C; HydroLight's tabulated values
carry no temperature correction. At 18 °C the correction lowers
`a_w` by up to about 1 %:

```@example hl
δ = absorption(water(), λ) ./ absorption(water(T = 20.0), λ) .- 1
@printf("largest change at 18 °C: %.2f %% at %.0f nm\n", 100minimum(δ), λ[argmin(δ)])
```

## Bulk IOPs through the tabulated route

The three scenarios can be built from tabulated coefficients: the bulk
absorption as a [`TabulatedAbsorber`](@ref) with `total = true`, and pure
water, phytoplankton and non-algal particles as
[`TabulatedScatterer`](@ref)s with their phase functions. The resulting
single-scattering albedo and backscatter fraction match the HydroLight
totals to the precision of the stored values.

```@example hl
cases  = ("clear", "median", "turbid")
series = Series[]
for (k, case) in enumerate(cases)
    c(v) = ref["$(case)_$v"]
    layer = OceanLayer(0.0, 1.0, AbstractOceanConstituent{Float64}[
        TabulatedAbsorber(λ, c("a"); total = true),
        TabulatedScatterer(λ, ref["bw"], RayleighWaterPhase{Float64}()),
        TabulatedScatterer(λ, c("bph"), FournierForandPhase{Float64}(backscatter_fraction = 0.01)),
        TabulatedScatterer(λ, c("bd"), FournierForandPhase{Float64}(backscatter_fraction = 0.018))])
    opt  = layer_optics(layer, λ; ℓ_max = 16)
    ϖ_hl = c("b") ./ (c("a") .+ c("b"))
    @printf("%-7s max |Δϖ| = %.1e   max |ΔB/B| = %.1e\n", case,
            maximum(abs, opt.ϖ .- ϖ_hl), maximum(abs, opt.B ./ (c("bb") ./ c("b")) .- 1))
    push!(series, Series("$case: OceanOptics", λ, opt.ϖ; slot = k - 1))
    push!(series, Series("$case: HydroLight", λ, ϖ_hl; points = true, slot = k - 1))
end
linechart(series; xlabel = "wavelength λ [nm]", ylabel = "single-scattering albedo ϖ",
          xname = "λ [nm]", ymin = 0, ymax = 1)
```

## What this establishes

* The `ZhangHu2009` seawater-scattering port agrees with HydroLight's
  input at 18 °C and salinity 35 within the data's precision
  (`test/test_hydrolight_reference.jl` checks 1 × 10⁻⁴).
* The bundled Pope & Fry table matches over 420–725 nm.
* The default pure-water absorption differs strongly below 420 nm and
  moderately above 700 nm. Choosing and bundling primary UV and near-IR
  tables (Mason, Cone & Fry 2016 or Morel et al. 2007; Kou et al. 1993),
  and reconciling the temperature-correction reference, are open work.
  Until then, pass HydroLight-matched water absorption as a
  `TabulatedAbsorber` when comparing with these simulations.
* Particle backscattering follows directly from the Fournier-Forand
  backscatter fractions; the phase-function shapes themselves are not part
  of the database and are not compared here.
