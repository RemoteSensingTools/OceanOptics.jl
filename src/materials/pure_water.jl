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

# Zhang, Hu & He (2009) seawater scattering — calibrated approximation
# -------------------------------------------------------------------
# NOTE: This is an *approximation* that reproduces the central published
#       Zhang-Hu value at λ=500 nm, S=35, T=20 °C (b_sw ≈ 0.00194 1/m),
#       NOT a verbatim port of the full paper. The full Zhang et al.
#       formulation requires isothermal compressibility, density
#       derivatives of the refractive index, a salt-concentration
#       polynomial, and a Cabannes depolarization correction; the MATLAB
#       reference is published at <https://www.seanoe.org/data/00318/42916/>.
#       Porting that reference verbatim is tracked for OceanOptics.jl
#       v0.2 (DESIGN §8.2 and open question in §12).
#
# Functional form:
#
#   b_sw(λ, T, S) = b_pure(λ) · [1 + α_S · S/35] · [1 − α_T · (T − 20)]
#
#   b_pure(λ)     = B₀ · (λ₀/λ)^n          (pure-water Einstein-Smoluchowski
#                                           power law, calibrated so
#                                           b_sw(500, 35, 20) = 0.00194)
#
# Calibration constants
#   B₀ = 0.00149 m⁻¹   at λ₀ = 500 nm, S = 0, T = 20 °C
#   n  = 4.32          (matches the Morel-1974 spectral exponent)
#   α_S = 0.30         salinity enhancement 0 → 35 PSU
#   α_T = 0.003 °C⁻¹   temperature sensitivity near 20 °C
#
# Against Zhang et al. (2009) Table I the pure-water baseline differs by
# a few percent across the visible; within the expected uncertainty of a
# one-parameter power-law fit to the paper's full physical model.
const _ZHH_B0    = 0.001492      # pure-water b_w at 500 nm, T=20°C, S=0   [1/m]
                                 # Calibrated so b_sw(500, 35, 20) = 0.00194
const _ZHH_λ0    = 500.0         # reference wavelength                     [nm]
const _ZHH_N     = 4.32          # spectral exponent
const _ZHH_αS    = 0.30          # fractional enhancement 0 → 35 PSU
const _ZHH_αT    = 0.003         # fractional decrement per °C above 20 °C

function _zhang_hu_2009_scattering(λ, T, S)
    # Range checks with `ifelse` cost us nothing here because T and S are
    # scalar construction-time fields of PureWater; any out-of-range value
    # is a user modeling mistake rather than an optimizer's Dual iterate.
    (0 ≤ S ≤ 45) ||
        @warn "Zhang-Hu: salinity S=$S [PSU] outside calibrated 0-45 range"
    (-2 ≤ T ≤ 40) ||
        @warn "Zhang-Hu: temperature T=$T [°C] outside calibrated -2 to 40 range"
    b_pure   = _ZHH_B0 * (_ZHH_λ0 / λ)^_ZHH_N
    salt_fac = 1 + _ZHH_αS * (S / 35)
    temp_fac = 1 - _ZHH_αT * (T - 20)
    return b_pure * salt_fac * temp_fac
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
