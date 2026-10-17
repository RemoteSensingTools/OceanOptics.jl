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
water(; a = IOCCG2018(), b = ZhangHu2009(), T = 18.0, S = 35.0) =
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

The default absorption model, [`IOCCG2018`](@ref), is the composite
HydroLight used here — its 350–750 nm table is taken from this database
— continued above 750 nm with Kou et al. (1993), the composite's own
source at 730–750 nm. The individual primary tables bundled with
OceanOptics compare as follows (each drawn over its own range; at 22 °C,
where no temperature correction applies):

```@example hl
in_range(lo, hi) = lo .≤ λ .≤ hi
a_hl  = ref["aw"]
models = [("PopeFry1997", PopeFry1997(), in_range(380, 727.5)),
          ("SmithBaker1981", SmithBaker1981(), in_range(350, 750)),
          ("Kou1993", Kou1993(), in_range(667, 750))]
a(m, sel) = absorption(water(a = m, T = 22.0), λ[sel])
linechart([Series("HydroLight (IOCCG 2018)", λ, a_hl; points = true, reference = true);
           [Series(name, λ[sel], a(m, sel); dash = name == "SmithBaker1981",
                   dot = name == "Kou1993") for (name, m, sel) in models]];
          xlabel = "wavelength λ [nm]", ylabel = "a_w(λ) [m⁻¹]", xname = "λ [nm]", ylog = true)
```

```@example hl
linechart([[Series(name, λ[sel], a(m, sel) ./ a_hl[sel]; dash = name == "SmithBaker1981",
                   dot = name == "Kou1993") for (name, m, sel) in models];
           Series("HydroLight", λ, ones(length(λ)); reference = true, dot = true)];
          xlabel = "wavelength λ [nm]", ylabel = "a_w(model) / a_w(HydroLight)",
          xname = "λ [nm]", ylog = true)
```

```@example hl
for (name, m, sel) in models
    for (lo, hi) in ((350, 415), (420, 725), (730, 750))
        s = sel .& in_range(lo, hi)
        any(s) || continue
        r = a(m, s) ./ a_hl[s]
        @printf("%-15s %3d–%3d nm: ratio %.3f – %.3f\n", name, max(lo, Int(minimum(λ[s]))),
                hi, extrema(r)...)
    end
end
```

Pope & Fry is the composite between 420 and 725 nm; Smith & Baker's
"clearest natural water" values are several times higher in the UV,
which is why they are no longer used in the default (the previous
default fell back to them below 380 nm); Kou et al. (1993) supplies
the composite's 730–750 nm values to within 1 %.

OceanOptics applies the Röttgers et al. (2014) temperature correction
from 550 nm upward, relative to the tables' 22 °C reference; HydroLight's
tabulated values carry no correction. At the database's 18 °C the
correction changes `a_w` by up to a few percent near the 740 nm O–H band:

```@example hl
δ = absorption(water(), λ) ./ absorption(water(T = 22.0), λ) .- 1
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
    ϖ, B = single_scattering_albedo(opt), backscatter_fraction(opt)
    @printf("%-7s max |Δϖ| = %.1e   max |ΔB/B| = %.1e\n", case,
            maximum(abs, ϖ .- ϖ_hl), maximum(abs, B ./ (c("bb") ./ c("b")) .- 1))
    push!(series, Series("$case: OceanOptics", λ, ϖ; slot = k - 1))
    push!(series, Series("$case: HydroLight", λ, ϖ_hl; points = true, slot = k - 1))
end
linechart(series; xlabel = "wavelength λ [nm]", ylabel = "single-scattering albedo ϖ",
          xname = "λ [nm]", ymin = 0, ymax = 1)
```

## What this establishes

* The `ZhangHu2009` seawater-scattering port agrees with HydroLight's
  input at 18 °C and salinity 35 within the data's precision
  (`test/test_hydrolight_reference.jl` checks 1 × 10⁻⁴).
* The bundled Pope & Fry table matches over 420–725 nm, and Kou et al.
  (1993) over 730–750 nm within 1 %; the default `IOCCG2018` model is
  the HydroLight composite by construction.
* Bulk single-scattering albedo and backscatter fraction are reproduced
  through the tabulated route, with Fournier-Forand particles on the
  Mobley et al. (2002) family that HydroLight uses.
* Particle phase-function shapes are not part of the database and are not
  compared here.
