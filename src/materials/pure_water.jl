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
# Either can be swapped independently.
#
# Absorption (each table throws outside its own range):
#   IOCCG2018       — 350–2500 nm, default: the IOCCG (2018) recommended
#                     composite (Morel et al. 2007, Pope & Fry 1997, Kou
#                     et al. 1993) as used for HydroLight by Loisel et al.
#                     (2023) over 350–750 nm, continued with Kou et al. 1993
#   PopeFry1997     — 380–727.5 nm
#   Kou1993         — 667–2500 nm
#   SmithBaker1981  — 200–800 nm, "clearest natural waters"; several times
#                     higher than modern pure-water values below ~420 nm
#
# The temperature correction a(T) = a(T_ref) + Ψ_T (T − T_ref) uses the
# Röttgers, McKee & Utschig (2014) coefficients from 550 nm upward (below,
# they are smaller than their uncertainty), referenced to each table's
# measurement temperature (22 °C for Pope & Fry and Kou); Smith & Baker give
# none and are not corrected.
#
# Scattering:
#   Morel1974       — power law, salinity-independent, legacy
#   ZhangHu2009     — default; depends on T, S, λ
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

"""
$(TYPEDEF)

Absorption of Smith & Baker (1981), 200–800 nm at 10 nm: upper bounds for the clearest
natural waters.

Below about 420 nm it is several times higher than modern pure-water measurements. No
temperature correction is applied.
"""
struct SmithBaker1981 <: AbstractPureWaterAbsorptionModel end

"""
$(TYPEDEF)

Absorption of Pope & Fry (1997), 380–727.5 nm at 2.5 nm, measured at 22 °C.
"""
struct PopeFry1997 <: AbstractPureWaterAbsorptionModel end

"""
$(TYPEDEF)

Absorption of Kou, Labrie & Chylek (1993), 667–2500 nm, measured at 22 °C.

The bundled table has 548 points with steps of 0.3 to 18 nm.
"""
struct Kou1993 <: AbstractPureWaterAbsorptionModel end

"""
$(TYPEDEF)

Pure-water absorption recommended by the IOCCG Protocol Series (2018); the default model.

It combines Morel et al. (2007) for 340–415 nm, Pope & Fry (1997) for 420–725 nm and Kou
et al. (1993) for 730–750 nm, at 5 nm over 350–750 nm as used for the HydroLight database
of Loisel et al. (2023), and is continued with Kou et al. (1993) up to 2500 nm. The
composite's 730–750 nm values differ from the bundled Kou table by up to 0.8 % (0.06 %
step at the 750 nm junction).
"""
struct IOCCG2018 <: AbstractPureWaterAbsorptionModel end

# -----------------------------------------------------------------------------
# Scattering-model tags
# -----------------------------------------------------------------------------

"""
$(TYPEDEF)

Pure-water scattering model: [`ZhangHu2009`](@ref) (the default) or [`Morel1974`](@ref).
"""
abstract type AbstractPureWaterScatteringModel end

"""
$(TYPEDEF)

Power law fitted to the seawater scattering of Morel (1974),

``b_w(λ) = 0.00193\\,(550/λ)^{4.32}`` `[1/m]`, `λ` in nm.

The formula represents seawater and does not depend on the `salinity`, `temperature`
or `depolarization` of the [`PureWater`](@ref): with `salinity = 0` it still returns
the seawater value, about 30 % above fresh water (compare [`ZhangHu2009`](@ref),
1.958·10⁻³ m⁻¹ at 500 nm for S = 0 against 2.913·10⁻³ here). It is a smooth fit, not
Morel's table: against the pure-sea-water values that Mobley et al. (1993, Problem 7)
quote from Smith & Baker (1981), 0.0063 m⁻¹ at 417 nm and 0.0032 m⁻¹ at 486 nm, it
gives 0.00638 (+1.3 %) and 0.00329 (+2.9 %); the printed values carry two digits. It
is 7–16 % above [`ZhangHu2009`](@ref) at S = 35 over 350–750 nm.
"""
struct Morel1974  <: AbstractPureWaterScatteringModel end

"""
$(TYPEDEF)

Scattering of seawater after Zhang, Hu & He (2009), a function of wavelength,
temperature and salinity.

The Rayleigh magnitude uses the `depolarization` of the [`PureWater`](@ref) (default
0.039) in the Cabannes factor `(6 + 6δ)/(6 − 7δ)` and in the angular integral.

No range of validity is checked: the routine returns values at any wavelength,
temperature and salinity. The refractive index inside it is the Quan & Fry (1995) fit
(400–700 nm, 0–30 °C, 0–35 PSU); outside those ranges the result is an extrapolation.
"""
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
# Default: IOCCG2018 absorption, Zhang-Hu scattering, S = 35 PSU, T = 20 °C
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
    "Depolarization ratio `δ` (unpolarized light, 90°) of the phase matrix and of the `ZhangHu2009` scattering magnitude"
    depolarization::FT
end

function PureWater{FT}(;
        absorption_model = IOCCG2018(),
        scattering_model = ZhangHu2009(),
        temperature      = 20,
        salinity         = 35,
        depolarization   = _WATER_DEPOLARIZATION,
) where {FT}
    _check_parameter_type(FT)
    _within(FT(depolarization), 0, 1) && FT(depolarization) < one(FT) ||
        throw(ArgumentError("PureWater: depolarization must be in [0, 1), got $depolarization"))
    PureWater{FT, typeof(absorption_model), typeof(scattering_model)}(
        absorption_model, scattering_model,
        FT(temperature), FT(salinity), FT(depolarization))
end

PureWater(; kwargs...) = PureWater{_parameter_type(values(kwargs))}(; kwargs...)

# Phase function interface — water scatters as depolarized Rayleigh.
scattering_phase(w::PureWater{FT}) where {FT} =
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

Return the pure-water absorption coefficient `[1/m]` at wavelength `λ` `[nm]`.

The table of the water's `absorption_model` is interpolated linearly; a wavelength
outside the tabulated range throws. The temperature correction
`a(T) = a(T_ref) + Ψ_T (T − T_ref)` uses the coefficients `Ψ_T` of Röttgers et al. (2014)
for `λ ≥ 550 nm` (up to 0.016 m⁻¹ K⁻¹ near 740 nm), relative to the measurement
temperature `T_ref` of the table. Salinity effects on absorption are not modelled.
"""
function absorption(w::PureWater{FT, A, S}, λ::Real) where {FT, A, S}
    T     = _result_type(FT, λ)
    x     = _wide(λ)
    a₀    = _interpolate_pure_water_absorption(A(), x)
    T_ref = _reference_temperature(A())
    T_ref === nothing && return convert(T, a₀)
    return convert(T, a₀ + _temperature_coefficient(x) * (_wide(w.temperature) - T_ref))
end

# Vectorized convenience
absorption(w::PureWater, λ::AbstractVector) = map(Base.Fix1(absorption, w), λ)

# -----------------------------------------------------------------------------
# Spectral reference tables
# -----------------------------------------------------------------------------
# Loaded once at package init from `data/` via `load_spectral_csv`. The
# provenance header on each CSV records the upstream citation, DOI, URL,
# date accessed, and any unit conversion applied at load time.

"Pure-water absorption, Pope & Fry (1997), 380–727.5 nm, 2.5 nm grid."
const _POPE_FRY_1997_TABLE    = load_spectral_csv(data_path("pope_fry_1997.csv"))

"Pure-water absorption, Smith & Baker (1981), 200–800 nm, 10 nm grid."
const _SMITH_BAKER_1981_TABLE = load_spectral_csv(data_path("smith_baker_1981.csv"))

"Liquid-water absorption, Kou, Labrie & Chylek (1993), 667–2500 nm."
const _KOU_1993_TABLE         = load_spectral_csv(data_path("kou_1993.csv"))

"IOCCG (2018) recommended pure-water absorption, 350–750 nm, 5 nm grid."
const _IOCCG_2018_TABLE       = load_spectral_csv(data_path("ioccg_2018_pure_water.csv"))

_interpolate_pure_water_absorption(::SmithBaker1981, λ) = linterp(_SMITH_BAKER_1981_TABLE, λ)
_interpolate_pure_water_absorption(::PopeFry1997,    λ) = linterp(_POPE_FRY_1997_TABLE,    λ)
_interpolate_pure_water_absorption(::Kou1993,        λ) = linterp(_KOU_1993_TABLE,         λ)
function _interpolate_pure_water_absorption(::IOCCG2018, λ)
    lo, hi = first(_IOCCG_2018_TABLE.λ), last(_KOU_1993_TABLE.λ)
    _within(λ, lo, hi) ||
        throw(DomainError(λ, "IOCCG2018 pure-water absorption covers $(lo)–$(hi) nm (the " *
                             "IOCCG 2018 table up to 750 nm, Kou et al. 1993 above); choose " *
                             "another absorption model or restrict the wavelength grid"))
    return λ ≤ last(_IOCCG_2018_TABLE.λ) ? linterp(_IOCCG_2018_TABLE, λ) : linterp(_KOU_1993_TABLE, λ)
end

# Measurement temperature [°C] each table refers to; `nothing`: no correction.
_reference_temperature(::SmithBaker1981) = nothing
_reference_temperature(::PopeFry1997)    = 22.0
_reference_temperature(::Kou1993)        = 22.0
_reference_temperature(::IOCCG2018)      = _reference_temperature(Kou1993())  # Pope & Fry / Kou in 600–800 nm

# Röttgers, McKee & Utschig (2014) temperature coefficient Ψ_T [1/(m K)],
# 550–2500 nm; zero below the first wavelength of the table (550 nm), where it is
# below the measurement uncertainty (and crosses zero near 550 nm, so there is no step).
const _ROTTGERS_2014_PSI_T = load_spectral_csv(data_path("rottgers_2014_temperature.csv"))
_temperature_coefficient(λ) =
    λ < _ROTTGERS_2014_PSI_T.lo ? zero(λ) : linterp(_ROTTGERS_2014_PSI_T, λ)

# =============================================================================
# Scattering dispatch — one method per model tag
# =============================================================================

"""
    scattering(w::PureWater, λ) -> FT

Return the scattering coefficient `[1/m]` of pure water or seawater at wavelength `λ`
`[nm]`, from the model named by the water's `scattering_model`.
"""
scattering(w::PureWater{FT, A, Morel1974},   λ::Real) where {FT, A} =
    convert(_result_type(FT, λ), _morel_1974_scattering(_wide(λ)))
scattering(w::PureWater{FT, A, ZhangHu2009}, λ::Real) where {FT, A} =
    convert(_result_type(FT, λ),
            _zhang_hu_2009_scattering(_wide(λ), _wide(w.temperature), _wide(w.salinity),
                                      _wide(w.depolarization)))

# Morel (1974) seawater power law b_w(λ) = b_ref (λ_ref/λ)^n: b_ref [1/m], λ_ref [nm], n.
const _MOREL_1974_B_REF    = 0.00193
const _MOREL_1974_λ_REF    = 550
const _MOREL_1974_EXPONENT = 4.32
# Un-annotated so that any `Real` wavelength (including `ForwardDiff.Dual`) promotes.
_morel_1974_scattering(λ) = _MOREL_1974_B_REF * (_MOREL_1974_λ_REF / λ)^_MOREL_1974_EXPONENT

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

# CODATA-2006 values of Avogadro's [1/mol] and Boltzmann's [J/K] constants as written in
# betasw_ZHH2009.m, kept for bit-reproduction of that code (the SI-2019 values differ by
# 3e-8 and 9e-7 relative).
const _ZHH_2009_AVOGADRO = 6.0221417930e23
const _ZHH_2009_BOLTZMANN = 1.3806503e-23
const _ZHH_M0  = 18e-3             # Molecular weight of water [kg/mol]

"Refractive index of seawater (Quan-Fry 1995 + Ciddor 1996) and `∂n_sw/∂S`."
function _zhh_refractive_index(λ_nm, Tc, S)
    λ_um = λ_nm / 1_000                      # µm, for Ciddor 1996
    n_air = 1 + (5_792_105 / (238.0185 - 1/λ_um^2) +
                 167_917   / (57.362   - 1/λ_um^2)) * 1e-8
    # Quan-Fry 1995 seawater refractive index (w.r.t. air), refractive_index.jl
    nsw_rel_air = _quan_fry_index(λ_nm, Tc, S)
    nsw         = nsw_rel_air * n_air         # absolute (in-vacuum) index
    dnds        = (1.779e-4 - 1.05e-6*Tc + 1.6e-8*Tc^2 + 0.01155/λ_nm) * n_air   # ∂n/∂S
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

"""
    _zhh_pmh(n)

Return the density derivative of Proutière, Megnassan & Hucteau,
`ρ ∂n²/∂ρ = (n²-1)·[1 + (2/3)(n²+2)·((n²-1)/(3n))²]`.
"""
function _zhh_pmh(n)
    n2 = n^2
    return (n2 - 1) * (1 + (2/3) * (n2 + 2) * (n/3 - 1/(3n))^2)
end

"""
    _zhang_hu_2009_scattering(λ_nm, Tc, S, δ)

Return the total seawater scattering coefficient `b_sw` `[1/m]` at wavelength `λ_nm`
`[nm]`, temperature `Tc` `[°C]`, practical salinity `S` `[PSU]` and depolarization
ratio `δ` (0.039 in Zhang, Hu & He 2009, after Farinato & Rowell 1976).

No range of validity is checked.
"""
function _zhang_hu_2009_scattering(λ_nm, Tc, S, δ)
    Tk = Tc + 273.15

    nsw, dnds = _zhh_refractive_index(λ_nm, Tc, S)
    β_T       = _zhh_isothermal_compressibility(Tc, S)
    ρ_sw      = _zhh_density(Tc, S)
    # ∂ln a/∂S has a √S term whose derivative is infinite at S = 0, where it multiplies
    # S = 0 in `flu_con`; evaluating it at `zero(S)` there gives the one-sided derivative
    # of the product (a dual S = 0 would otherwise give 0·∞ = NaN). Same value.
    dlna_dS   = _zhh_dlna_dS(Tc, _within(S, 0, 0) ? zero(S) : S)
    DFRI      = _zhh_pmh(nsw)

    λm        = λ_nm * 1e-9                 # nm → m
    cabannes  = (6 + 6δ) / (6 - 7δ)

    # Density-fluctuation contribution to β(90°)
    β_df = (π^2 / 2) * λm^(-4) * _ZHH_2009_BOLTZMANN * Tk * β_T * DFRI^2 * cabannes

    # Concentration-fluctuation contribution to β(90°)
    flu_con = S * _ZHH_M0 * dnds^2 / ρ_sw / (-dlna_dS) / _ZHH_2009_AVOGADRO
    β_cf    = 2 * π^2 * λm^(-4) * nsw^2 * flu_con * cabannes

    β_90    = β_df + β_cf
    # Angular integration of the Cabannes phase function → b
    return (8π / 3) * β_90 * (2 + δ) / (1 + δ)
end

scattering(w::PureWater, λ::AbstractVector) = map(Base.Fix1(scattering, w), λ)

# =============================================================================
# Backscatter fraction — Rayleigh-like, exactly 1/2 for any depolarization
# =============================================================================

# Backscatter fraction of a depolarized-Rayleigh molecular scatterer is
# exactly 1/2 for any depolarization (the phase function is symmetric
# forward/back), so `b_b = b/2`.
backscattering(w::PureWater{FT}, λ::Real) where {FT} = FT(0.5) * scattering(w, λ)
backscattering(w::PureWater, λ::AbstractVector) = map(Base.Fix1(backscattering, w), λ)
