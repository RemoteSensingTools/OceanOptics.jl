# =============================================================================
# Pure (sea)water material and its IOPs
# =============================================================================
#
# A `PureWater` instance carries *two* model-choice tags because absorption
# and scattering come from independent literatures that have evolved at
# different rates:
#
#   PureWater{FT, A <: AbstractPureWaterAbsorptionModel,
#                 S <: AbstractPureWaterScatteringModel}
#
# This lets you, for example, mix Pope & Fry (1997) visible-light absorption
# with Zhang, Hu & He (2009) salinity-aware scattering — which is exactly
# the combination that HydroLight 5 uses by default. Or swap either
# independently as the literature advances.
#
# Temperature/salinity corrections are carried as optional fields on the
# absorption-model tag type. For the scattering model, salinity is a
# *required* input (Zhang et al. 2009 cannot produce a pure-water result
# without knowing S).
#
# Shipped model tags
# ------------------
#
# Absorption:
#   SmithBaker1981        — 200–800 nm, low-precision, good NIR
#   PopeFry1997           — 380–700 nm, visible gold standard
#   MasonConeFry2016      — 250–550 nm, modern UV measurements
#   CombinedPopeFryMason  — piecewise: Mason 250–380, PopeFry 380–700, SmithBaker elsewhere
#                            (recommended default; equivalent to HydroLight's "aw" table)
#
# Scattering:
#   Morel1974            — simple power law, salinity-independent, legacy
#   ZhangHu2009          — modern default, depends on T, S, λ
# =============================================================================

# -----------------------------------------------------------------------------
# Absorption-model tags and metadata
# -----------------------------------------------------------------------------

"""
$(TYPEDEF)

Pure-water absorption dataset.  Each subtype represents a published
spectral absorption table; concrete `absorption(::PureWater{FT, M, …}, λ)`
methods dispatch on `M` to pick the right one.
"""
abstract type AbstractPureWaterAbsorptionModel end

"Smith & Baker (1981). 200–800 nm, 10 nm resolution. Superseded in visible."
struct SmithBaker1981 <: AbstractPureWaterAbsorptionModel end

"Pope & Fry (1997). 380–700 nm, 2.5 nm resolution. Visible gold standard."
struct PopeFry1997     <: AbstractPureWaterAbsorptionModel end

"Mason, Cone & Fry (2016). 250–550 nm. Modern UV definitive."
struct MasonConeFry2016 <: AbstractPureWaterAbsorptionModel end

"""
Combined piecewise model: Mason-Cone-Fry below 380 nm, Pope-Fry in
380–700 nm, Smith-Baker elsewhere (mainly NIR). Matches the default
HydroLight "aw" table.
"""
struct CombinedPopeFryMason <: AbstractPureWaterAbsorptionModel end

# -----------------------------------------------------------------------------
# Scattering-model tags
# -----------------------------------------------------------------------------

"""
$(TYPEDEF)

Pure-water scattering model. `ZhangHu2009` is the modern default; the
older `Morel1974` is kept for reproducing legacy codes.
"""
abstract type AbstractPureWaterScatteringModel end

"Morel (1974). Salinity-independent: b_w(λ) = 0.00193 (550/λ)^4.32."
struct Morel1974  <: AbstractPureWaterScatteringModel end

"Zhang, Hu & He (2009). Includes salinity and temperature dependence. 200–1100 nm."
struct ZhangHu2009 <: AbstractPureWaterScatteringModel end

# -----------------------------------------------------------------------------
# Concrete PureWater type
# -----------------------------------------------------------------------------

"""
$(TYPEDEF)

Pure (sea)water as an optical constituent.

# Fields
$(TYPEDFIELDS)

# Examples

```julia
# Default: piecewise absorption, Zhang-Hu scattering, S=35 PSU, T=20 °C
water = PureWater{Float64}()

# Legacy Pope-Fry + Morel combination
water = PureWater{Float64}(
    absorption_model = PopeFry1997(),
    scattering_model = Morel1974(),
)

# Cold, freshwater Arctic case
water = PureWater{Float64}(temperature = 2.0, salinity = 2.0)
```
"""
struct PureWater{FT,
                 A<:AbstractPureWaterAbsorptionModel,
                 S<:AbstractPureWaterScatteringModel} <: AbstractAbsorbingScatterer{FT}
    "Absorption-model choice (tag)"
    absorption_model::A
    "Scattering-model choice (tag)"
    scattering_model::S
    "Temperature [°C]"
    temperature::FT
    "Salinity [PSU] (Practical Salinity Unit)"
    salinity::FT
    "Depolarization ratio for the Rayleigh-type scattering phase function"
    depolarization::FT
end

function PureWater{FT}(;
        absorption_model = CombinedPopeFryMason(),
        scattering_model = ZhangHu2009(),
        temperature      = FT(20),
        salinity         = FT(35),
        depolarization   = FT(0.039),
) where {FT}
    PureWater{FT, typeof(absorption_model), typeof(scattering_model)}(
        absorption_model, scattering_model,
        FT(temperature), FT(salinity), FT(depolarization))
end

PureWater(; kwargs...) = PureWater{Float64}(; kwargs...)

# Phase function interface — water scatters as depolarized Rayleigh.
phase_function(w::PureWater{FT}) where {FT} =
    RayleighWaterPhase{FT}(w.depolarization)

# =============================================================================
# Absorption dispatch — one method per model tag
# =============================================================================
#
# The numerical tables are kept in a companion module `PureWaterData`.
# Each `absorption(::PureWater{FT, M, S}, λ)` method delegates to
# `_absorption_table(M)` which returns the (λ, a) vectors for that model.
# That lets us add a new model by adding one method and one table, without
# touching `PureWater` itself.
# =============================================================================

"""
    absorption(w::PureWater, λ) -> FT

Pure-water absorption coefficient [1/m] at wavelength `λ` [nm].

Linearly interpolates the spectral table of the configured
`absorption_model`. A small temperature correction (Pegau, Gray &
Zaneveld 1997) is applied if the model supports it; salinity affects
absorption only above 700 nm.
"""
function absorption(w::PureWater{FT, A, S}, λ::Real) where {FT, A, S}
    # Temperature correction (Pegau, Gray & Zaneveld 1997) is blended in via
    # the Pegau slope table, which is explicitly zero outside the calibrated
    # 600–800 nm window — so we always-apply it without branching on λ.
    # Branch-free form is required for AD transparency (DESIGN §12.5).
    a₀      = _interpolate_pure_water_absorption(A(), λ)
    ∂a_∂T   = _pegau_temperature_slope(λ)
    return a₀ + ∂a_∂T * (w.temperature - FT(20))
end

# Vectorized convenience
absorption(w::PureWater, λ::AbstractVector) = [absorption(w, λi) for λi in λ]

# -----------------------------------------------------------------------------
# Spectral reference tables
# -----------------------------------------------------------------------------
# Loaded once at package init from `data/` via `load_spectral_csv`. The
# provenance header on each CSV records the upstream citation, DOI, URL,
# date accessed, and any unit conversion applied at load time — see
# DESIGN §3.4 for the provenance-first rationale.

"Pure-water absorption, Pope & Fry (1997), 380–727.5 nm, 2.5 nm grid."
const _POPE_FRY_1997_TABLE    = load_spectral_csv(data_path("pope_fry_1997.csv"))

"Pure-water absorption, Smith & Baker (1981), 200–800 nm, 10 nm grid."
const _SMITH_BAKER_1981_TABLE = load_spectral_csv(data_path("smith_baker_1981.csv"))

# Mason, Cone & Fry (2016) UV extension — optional. The Table 2 transcription
# is NOT shipped with this release (DESIGN §8.1). `_MASON_CONE_FRY_AVAILABLE`
# is a compile-time snapshot of data-file presence; drop
# `data/mason_cone_fry_2016.csv` into place and restart Julia to activate.
const _MASON_CONE_FRY_2016_PATH  = data_path("mason_cone_fry_2016.csv")
const _MASON_CONE_FRY_AVAILABLE  = isfile(_MASON_CONE_FRY_2016_PATH)
const _MASON_CONE_FRY_2016_TABLE = Ref{Union{Nothing, SpectralTable{Float64}}}(nothing)

function _mason_cone_fry_2016_table()
    tbl = _MASON_CONE_FRY_2016_TABLE[]
    tbl === nothing || return tbl
    isfile(_MASON_CONE_FRY_2016_PATH) || error(
        "MasonConeFry2016 pure-water absorption requires ",
        "`data/mason_cone_fry_2016.csv`, which is not bundled with this ",
        "OceanOptics.jl release. Transcribe Table 2 of Mason, Cone & Fry ",
        "(2016), Appl. Opt. 55, 7163 (see DESIGN §8.1). Until then, use ",
        "PopeFry1997 or CombinedPopeFryMason (which falls back to Smith-Baker ",
        "below 380 nm when Mason data are absent).")
    tbl = load_spectral_csv(_MASON_CONE_FRY_2016_PATH)
    _MASON_CONE_FRY_2016_TABLE[] = tbl
    return tbl
end

_interpolate_pure_water_absorption(::SmithBaker1981,   λ) = linterp(_SMITH_BAKER_1981_TABLE, λ)
_interpolate_pure_water_absorption(::PopeFry1997,      λ) = linterp(_POPE_FRY_1997_TABLE,    λ)
_interpolate_pure_water_absorption(::MasonConeFry2016, λ) = linterp(_mason_cone_fry_2016_table(), λ)

function _interpolate_pure_water_absorption(::CombinedPopeFryMason, λ)
    # Piecewise HydroLight-style "aw" composite. The UV segment (λ < 380 nm)
    # prefers Mason-Cone-Fry (2016); Smith-Baker (1981), which spans 200-800
    # nm, serves as a standing fallback until the Mason table is bundled.
    if λ < 380
        return _MASON_CONE_FRY_AVAILABLE ?
                   linterp(_mason_cone_fry_2016_table(), λ) :
                   linterp(_SMITH_BAKER_1981_TABLE, λ)
    elseif λ ≤ 700
        return linterp(_POPE_FRY_1997_TABLE, λ)
    else
        return linterp(_SMITH_BAKER_1981_TABLE, λ)
    end
end

# Pegau, Gray & Zaneveld (1997) temperature slope `∂a/∂T` `[1/m/°C]`.
# The bundled CSV pads the 600–800 nm calibrated window with sentinel zeros
# at 0 nm and 10 000 nm so this interpolator returns exactly zero
# temperature sensitivity anywhere outside [600, 800] nm — letting
# `absorption` always-apply the correction without branching on λ.
const _PEGAU_1997_TABLE = load_spectral_csv(data_path("pegau_1997_temperature.csv"))
_pegau_temperature_slope(λ) = linterp(_PEGAU_1997_TABLE, λ)

# =============================================================================
# Scattering dispatch — one method per model tag
# =============================================================================

"""
    scattering(w::PureWater, λ) -> FT

Pure-(sea)water scattering coefficient [1/m] at wavelength `λ` [nm].
Dispatches on the `scattering_model` tag.
"""
scattering(w::PureWater{FT, A, Morel1974},   λ::Real) where {FT, A} = _morel_1974_scattering(λ)
scattering(w::PureWater{FT, A, ZhangHu2009}, λ::Real) where {FT, A} =
    _zhang_hu_2009_scattering(λ, w.temperature, w.salinity)

# Morel (1974) pure-water scattering, salinity-independent power law
#   b_w(λ) = 0.00193 · (550/λ)^4.32    [1/m],  with λ in [nm]
# Left un-annotated for AD transparency — integer/float literals promote
# against any `Real` (including `ForwardDiff.Dual`) input.
_morel_1974_scattering(λ) = 0.00193 * (550 / λ)^4.32

# Zhang, Hu & He (2009) seawater scattering — faithful port
# ----------------------------------------------------------
#
# Verbatim Julia port of the Zhang, Hu & He (2009) `betasw_ZHH2009.m`
# reference MATLAB implementation:
#
#   Zhang, Hu & He (2009), "Scattering by pure seawater: Effect of
#   salinity," Opt. Express 17, 5698–5710. doi:10.1364/OE.17.005698.
#
# Auxiliary formulations (see paper Table 1 and §2):
#
#   * n_sw(λ, T, S)   Quan & Fry (1995) Appl. Opt. 34, 3477 — seawater
#                     refractive index, scaled by Ciddor (1996) air
#                     refractive index to obtain absolute (in-vacuum) n.
#   * ρ_sw(T, S)      UNESCO (1981) equation of state.
#   * β_T(T, S)       Millero (1980) Deep-Sea Res. secant-bulk modulus.
#   * ∂ ln a₀/∂ S    Millero & Leung (1976) Am. J. Sci., polynomial fit.
#   * ρ ∂n²/∂ρ        Proutière, Megnassan & Hucteau (1992) J. Phys. Chem.
#                     — the "PMH" density derivative of the refractive index.
#   * f(δ) = (6+6δ)/(6-7δ)    Cabannes depolarization factor.
#
# The total scattering coefficient closes by the standard angular
# integration of the Cabannes phase function:
#
#   β(90°) = β_density(90°) + β_concentration(90°)
#   b_sw    = (8π/3) · β(90°) · (2+δ) / (1+δ)
#
# Agreement with Morel's (1966, 1968) measurements reported in the paper
# is within 1 % over 350–700 nm at S = 38.4 ‰, T = 20 °C (Fig. 3).
#
# Reference MATLAB source:
#   https://raw.githubusercontent.com/ooici/ion-functions/master/
#           ion_functions/data/matlab_scripts/flort/betasw_ZHH2009.m
# =============================================================================

const _ZHH_Na  = 6.0221417930e23   # Avogadro's constant [1/mol]
const _ZHH_Kbz = 1.3806503e-23     # Boltzmann constant  [J/K]
const _ZHH_M0  = 18e-3             # Molecular weight of water [kg/mol]

"Refractive index of seawater (Quan-Fry 1995 + Ciddor 1996) and `∂n_sw/∂S`."
function _zhh_refractive_index(λ_nm, Tc, S)
    λ_um = λ_nm / 1_000                      # µm, for Ciddor 1996
    n_air = 1 + (5_792_105 / (238.0185 - 1/λ_um^2) +
                 167_917   / (57.362   - 1/λ_um^2)) * 1e-8
    # Quan-Fry 1995 seawater refractive index (w.r.t. air)
    n0 = 1.31405;   n1 = 1.779e-4;  n2 = -1.05e-6;  n3 = 1.6e-8;    n4 = -2.02e-6
    n5 = 15.868;    n6 = 0.01155;   n7 = -0.00423;  n8 = -4382.0;   n9 = 1.1455e6
    nsw_rel_air = n0 + (n1 + n2*Tc + n3*Tc^2)*S + n4*Tc^2 +
                  (n5 + n6*S + n7*Tc)/λ_nm + n8/λ_nm^2 + n9/λ_nm^3
    nsw         = nsw_rel_air * n_air         # absolute (in-vacuum) index
    dnds        = (n1 + n2*Tc + n3*Tc^2 + n6/λ_nm) * n_air
    return nsw, dnds
end

"Isothermal compressibility of seawater `β_T` `[1/Pa]`, Millero (1980)."
function _zhh_isothermal_compressibility(Tc, S)
    # Pure-water secant bulk modulus
    kw = 19652.21 + 148.4206*Tc - 2.327105*Tc^2 + 1.360477e-2*Tc^3 - 5.155288e-5*Tc^4
    # Salinity corrections
    a0 = 54.6746 - 0.603459*Tc + 1.09987e-2*Tc^2 - 6.167e-5*Tc^3
    b0 = 7.944e-2 + 1.6483e-2*Tc - 5.3009e-4*Tc^2
    Ks = kw + a0*S + b0*S^1.5
    # Ks in bar; β_T = 1/Ks in 1/bar → multiply by 1e-5 to get 1/Pa
    return 1 / Ks * 1e-5
end

"Seawater density `ρ_sw` `[kg/m³]`, UNESCO (1981) EOS."
function _zhh_density(Tc, S)
    a0=8.24493e-1;  a1=-4.0899e-3; a2=7.6438e-5;  a3=-8.2467e-7; a4=5.3875e-9
    a5=-5.72466e-3; a6=1.0227e-4;  a7=-1.6546e-6; a8=4.8314e-4
    b0=999.842594;  b1=6.793952e-2; b2=-9.09529e-3
    b3=1.001685e-4; b4=-1.120083e-6; b5=6.536332e-9
    ρ_w = b0 + b1*Tc + b2*Tc^2 + b3*Tc^3 + b4*Tc^4 + b5*Tc^5
    return ρ_w + ((a0 + a1*Tc + a2*Tc^2 + a3*Tc^3 + a4*Tc^4)*S +
                  (a5 + a6*Tc + a7*Tc^2)*S^1.5 +
                  a8*S^2)
end

"Partial derivative of ln(water activity) w.r.t. salinity, Millero-Leung (1976)."
function _zhh_dlna_dS(Tc, S)
    return (-5.58651e-4 + 2.40452e-7*Tc - 3.12165e-9*Tc^2 + 2.40808e-11*Tc^3) +
           1.5 * (1.79613e-5 - 9.9422e-8*Tc + 2.08919e-9*Tc^2 - 1.39872e-11*Tc^3) * S^0.5 +
           2 * (-2.31065e-6 - 1.37674e-9*Tc - 1.93316e-11*Tc^2) * S
end

"Proutière-Megnassan-Hucteau density derivative `ρ ∂n²/∂ρ = (n²-1)·[1 + (2/3)(n²+2)·((n²-1)/(3n))²]`."
function _zhh_pmh(n)
    n2 = n^2
    return (n2 - 1) * (1 + (2/3) * (n2 + 2) * (n/3 - 1/(3n))^2)
end

"""
    _zhang_hu_2009_scattering(λ_nm, Tc, S; δ = 0.039) -> FT

Total seawater scattering coefficient `b_sw` `[1/m]` at wavelength
`λ_nm` `[nm]`, temperature `Tc` `[°C]` and practical salinity `S` `[PSU]`,
with depolarization ratio `δ` (Farinato & Rowell 1976 value 0.039).

Valid range: `λ ∈ [220, 1100] nm`, `T ∈ [0, 30] °C`, `S ∈ [0, 40] PSU`.
"""
function _zhang_hu_2009_scattering(λ_nm, Tc, S)
    δ  = 0.039                              # Farinato & Rowell (1976)
    Tk = Tc + 273.15

    nsw, dnds = _zhh_refractive_index(λ_nm, Tc, S)
    β_T       = _zhh_isothermal_compressibility(Tc, S)
    ρ_sw      = _zhh_density(Tc, S)
    dlna_dS   = _zhh_dlna_dS(Tc, S)
    DFRI      = _zhh_pmh(nsw)

    λm        = λ_nm * 1e-9                 # nm → m
    cabannes  = (6 + 6δ) / (6 - 7δ)

    # Density-fluctuation contribution to β(90°)
    β_df = (π^2 / 2) * λm^(-4) * _ZHH_Kbz * Tk * β_T * DFRI^2 * cabannes

    # Concentration-fluctuation contribution to β(90°)
    flu_con = S * _ZHH_M0 * dnds^2 / ρ_sw / (-dlna_dS) / _ZHH_Na
    β_cf    = 2 * π^2 * λm^(-4) * nsw^2 * flu_con * cabannes

    β_90    = β_df + β_cf
    # Angular integration of the Cabannes phase function → b
    return (8π / 3) * β_90 * (2 + δ) / (1 + δ)
end

scattering(w::PureWater, λ::AbstractVector) = [scattering(w, λi) for λi in λ]

# =============================================================================
# Backscatter fraction — Rayleigh-like, exactly 1/2 for any depolarization
# =============================================================================

# Backscatter fraction of a depolarized-Rayleigh molecular scatterer is
# exactly 1/2 for any depolarization (the phase function is symmetric
# forward/back), so `b_b = b/2`.
backscattering(w::PureWater{FT}, λ::Real) where {FT} = FT(0.5) * scattering(w, λ)
backscattering(w::PureWater, λ::AbstractVector) = [backscattering(w, λi) for λi in λ]
