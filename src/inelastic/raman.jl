# =============================================================================
# Water Raman scattering as an inelastic process
# =============================================================================
#
# Water Raman scattering (the O-H stretch band of liquid water) shifts light
# from the excitation wavenumber ν' to ν' − Δν. The process is described by
#
#   * the Stokes shift Δν and a Gaussian band of standard deviation σ_ν in
#     wavenumber,
#   * the Raman scattering coefficient b_R(λ') = b_R(λ_ref)·(λ_ref/λ')^n,
#   * the convention of b_R (energy or photon conversion, a type parameter),
#   * a Rayleigh-type phase function set by the depolarization ratio ρ.
#
# In Fell's factorization b_R takes the role of a^I(λ') in b^I = a^I · f^I.
# =============================================================================

"""
$(TYPEDEF)

Water Raman scattering as an inelastic process: a Gaussian emission band in wavenumber
centred at `ν' − shift_cm_inv`, with the scattering coefficient

``b_R(λ') = b_R(λ_\\mathrm{ref}) (λ_\\mathrm{ref}/λ')^n`` `[1/m]`.

# Convention of `b_R`

The type parameter `U` (field `conversion`) states what `b_R` converts
(see [`conversion_units`](@ref)):

- [`EnergyConversion`](@ref) (default): `b_R(λ') E₀(λ')` is the power emitted over the
  whole Raman band per unit volume and per nm of excitation wavelength. The energy
  kernel is `K_E(λ', λ) = b_R(λ') f_R(λ', λ)` with `∫ K_E dλ = b_R(λ')`, and the photon
  kernel is `K_N = K_E λ/λ'`. This is the definition of Mobley et al. (1993, Appl. Opt.
  32, 7484), whose Problem 7 source is `b_R E₀/(4π)` for the 417 → 486 nm transfer.
- [`PhotonConversion`](@ref): `b_R` counts photons, `K_N = b_R f_R` with
  `∫ K_N dλ = b_R(λ')`, and `K_E = K_N λ'/λ`. Its energy source is 12 % (excitation at
  350 nm) to 24 % (700 nm) below that of the energy convention for the same `b_R`.

# Angular distribution

``p_R(\\cos Θ) = \\frac{1}{4π}\\left[1 + β_2 P_2(\\cos Θ)\\right]``,
``β_2 = \\frac{1 - ρ}{2 (1 + 2ρ)}``

(Mobley et al. 1993, Eq. 9), with `ρ = linear_depolarization` the depolarization ratio
for linearly polarized excitation; `β₂ = 0.3097` for `ρ = 0.17`. See
[`raman_phase`](@ref) for the phase function and phase matrix.

# Fields
$(TYPEDFIELDS)

# Parameters and their sources

| field | default | source |
|:------|:--------|:-------|
| `shift_cm_inv`, `σ_cm_inv` | 3400, 85 cm⁻¹ | single-Gaussian band attributed to Haltrin & Kattawar (1993, Appl. Opt. 32, 5356) |
| `cross_section_ref` at `cross_section_λ_ref` | 2.7·10⁻⁴ m⁻¹ at 488 nm | attributed to Bartlett et al. (1998, JGR 103, 17875); the HydroLight runs of Loisel et al. (2023) are reproduced with 2.6·10⁻⁴ m⁻¹ |
| `cross_section_n` | 5.5 | attributed to Bartlett et al. (1998), for `b_R` in energy units |
| `linear_depolarization` | 0.17 | Mobley et al. (1993), Eq. 9 and Problem 7 |
| `excitation_range` | (350, 700) nm | choice of this package; emission then starts at 397 nm |

# Examples
```jldoctest
julia> ram = WaterRaman();

julia> conversion_units(ram)
EnergyConversion()

julia> round(excitation_absorption(ram, 417.0); sigdigits = 4)
0.0006411

julia> WaterRaman(conversion = PhotonConversion()) isa WaterRaman{Float64, PhotonConversion}
true
```
"""
struct WaterRaman{FT<:AbstractFloat, U<:AbstractConversionUnits} <: AbstractOceanInelasticProcess{FT}
    "Stokes shift between excitation and emission wavenumbers `[1/cm]`"
    shift_cm_inv::FT
    "Gaussian standard deviation of the emission band in wavenumber `[1/cm]`"
    σ_cm_inv::FT
    "Raman scattering coefficient `b_R` at the reference wavelength `[1/m]`"
    cross_section_ref::FT
    "Reference wavelength of the power law `[nm]`"
    cross_section_λ_ref::FT
    "Exponent `n` of `b_R(λ') = b_ref (λ_ref / λ')^n`"
    cross_section_n::FT
    "Depolarization ratio `ρ` for linearly polarized excitation, observed at 90°"
    linear_depolarization::FT
    "Excitation-wavelength integration window `[nm]`"
    excitation_range::Tuple{FT, FT}
    "Quantity that `b_R` converts: `EnergyConversion()` or `PhotonConversion()`"
    conversion::U
end

"""
    WaterRaman{FT}(; shift_cm_inv=3400, σ_cm_inv=85,
                     cross_section_ref=2.7e-4, cross_section_λ_ref=488,
                     cross_section_n=5.5, linear_depolarization=0.17,
                     excitation_range=(350, 700), conversion=EnergyConversion())

Construct the water Raman process; `WaterRaman(; ...)` is `WaterRaman{Float64}(; ...)`.

# Arguments
- `shift_cm_inv`: Stokes shift `[1/cm]`.
- `σ_cm_inv`: standard deviation of the Gaussian band `[1/cm]`, positive.
- `cross_section_ref`: `b_R` `[1/m]` at `cross_section_λ_ref` `[nm]`, non-negative.
- `cross_section_n`: exponent of the power law.
- `linear_depolarization`: depolarization ratio `ρ` for linearly polarized excitation,
  in `[0, 3/4]`; the natural-light ratio is `ρ_n = 2ρ/(1+ρ)`.
- `excitation_range`: `(lo, hi)` `[nm]`.
- `conversion`: [`EnergyConversion`](@ref)`()` or [`PhotonConversion`](@ref)`()`.
"""
function WaterRaman{FT}(;
        shift_cm_inv                            = FT(_OH_STRETCH_SHIFT_CM_INV),
        σ_cm_inv                                = FT(_OH_STRETCH_SIGMA_CM_INV),
        cross_section_ref                       = FT(2.7e-4),
        cross_section_λ_ref                     = FT(488.0),
        cross_section_n                         = FT(5.5),
        linear_depolarization                   = FT(0.17),
        excitation_range::Tuple{<:Real, <:Real} = (FT(350), FT(700)),
        conversion::AbstractConversionUnits     = EnergyConversion(),
) where {FT<:AbstractFloat}
    FT(σ_cm_inv) > zero(FT) ||
        throw(ArgumentError("WaterRaman: σ_cm_inv must be positive, got $σ_cm_inv"))
    FT(cross_section_ref) ≥ zero(FT) ||
        throw(ArgumentError("WaterRaman: cross_section_ref must be non-negative, got $cross_section_ref"))
    FT(cross_section_λ_ref) > zero(FT) ||
        throw(ArgumentError("WaterRaman: cross_section_λ_ref must be positive, got $cross_section_λ_ref"))
    # ρ = 3/4 is the largest ratio a Raman band can have (fully depolarized line).
    zero(FT) ≤ FT(linear_depolarization) ≤ FT(0.75) ||
        throw(ArgumentError("WaterRaman: linear_depolarization (the ratio for linearly " *
                            "polarized excitation) must be in [0, 3/4], got " *
                            "$linear_depolarization; a natural-light ratio ρ_n converts as " *
                            "ρ = ρ_n/(2 − ρ_n)"))
    excitation_range[1] < excitation_range[2] ||
        throw(ArgumentError("WaterRaman: excitation_range must be (lo < hi)"))
    return WaterRaman{FT, typeof(conversion)}(
        FT(shift_cm_inv),
        FT(σ_cm_inv),
        FT(cross_section_ref),
        FT(cross_section_λ_ref),
        FT(cross_section_n),
        FT(linear_depolarization),
        (FT(excitation_range[1]), FT(excitation_range[2])),
        conversion)
end

WaterRaman(; kwargs...) = WaterRaman{Float64}(; kwargs...)

# -----------------------------------------------------------------------------
# Trait implementations
# -----------------------------------------------------------------------------

conversion_units(p::WaterRaman) = p.conversion

"""
    excitation_absorption(p::WaterRaman, λ_prime) -> FT

Return the Raman scattering coefficient `[1/m]` at excitation wavelength `λ'` `[nm]`,

```
b_R(λ') = b_ref · (λ_ref / λ')^n
```

an energy or a photon conversion coefficient according to
[`conversion_units`](@ref)`(p)`.
"""
excitation_absorption(p::WaterRaman, λ_prime::Real) =
    p.cross_section_ref * (p.cross_section_λ_ref / λ_prime)^p.cross_section_n

"""
    emission(p::WaterRaman, λ_prime, λ) -> FT

Return the Raman redistribution density `[1/nm]`: a Gaussian in wavenumber centered at
`ν' − shift`, times the Jacobian `|dν/dλ| = 10⁷/λ²` so that `∫ f(λ', λ) dλ = 1`.
"""
emission(p::WaterRaman{FT}, λ_prime::Real, λ::Real) where {FT} =
    convert(_result_type(FT, λ_prime, λ),
            _wavenumber_shift_density(_wide(p.shift_cm_inv), _wide(p.σ_cm_inv),
                                      _wide(λ_prime), _wide(λ)))

excitation_range(p::WaterRaman) = p.excitation_range

"""
    is_isotropic(p::WaterRaman) -> Bool

Return `false`: water Raman emission follows the Rayleigh-type phase function of
[`raman_phase`](@ref).
"""
is_isotropic(p::WaterRaman) = false

"""
    natural_depolarization(p::WaterRaman) -> FT

Return the depolarization ratio for unpolarized (natural) excitation,
`ρ_n = 2ρ/(1 + ρ)`, with `ρ = p.linear_depolarization` the ratio for linearly polarized
excitation.

`ρ_n = 0.2906` for the default `ρ = 0.17`. `ρ_n` is the ratio that
[`RayleighWaterPhase`](@ref) takes.
"""
natural_depolarization(p::WaterRaman) =
    2 * p.linear_depolarization / (1 + p.linear_depolarization)

"""
    raman_phase(p::WaterRaman) -> RayleighWaterPhase

Return the angular distribution of Raman emission as the Rayleigh-type phase function
with the natural-light depolarization ratio [`natural_depolarization`](@ref)`(p)`.

In the normalization of this package (`½ ∫ p dμ = 1`) it is

``p(\\cos Θ) = 1 + β_2 P_2(\\cos Θ)``,
``β_2 = \\frac{1 - ρ_n}{2 + ρ_n} = \\frac{1 - ρ}{2 (1 + 2ρ)}``

which is `4π` times Eq. 9 of Mobley et al. (1993),
`(3/16π) (1+3ρ)/(1+2ρ) [1 + (1−ρ)/(1+3ρ) cos²Θ]`. `phase_function_value`,
`phase_function_moments` and `greek_expansion(...; polarized = true)` of the returned
object give the phase function, its Legendre moments and the Greek coefficients of the
phase matrix (the Rayleigh matrix with depolarization `ρ_n`).

# Examples
```jldoctest
julia> round(phase_function_moments(raman_phase(WaterRaman()), 2)[3]; digits = 4)
0.3097
```
"""
raman_phase(p::WaterRaman{FT}) where {FT} =
    RayleighWaterPhase{FT}(depolarization = natural_depolarization(p))

function Base.show(io::IO, p::WaterRaman{FT}) where {FT}
    print(io, "WaterRaman{$FT}(Δν=", p.shift_cm_inv, " 1/cm, σ=",
              p.σ_cm_inv, " 1/cm, b_R(",
              Int(round(p.cross_section_λ_ref)), ")=",
              p.cross_section_ref, " 1/m, ", nameof(typeof(p.conversion)), ")")
end

"""
    raman_peak_wavelength(p::WaterRaman, λ_prime::Real) -> FT

Return the wavelength `[nm]` of the emission peak for excitation at `λ_prime` `[nm]`:

```
ν_emission = ν' - Δν,   λ_peak = 10⁷ / ν_emission
```

with wavenumbers in `1/cm`.

# Examples
```jldoctest
julia> round(raman_peak_wavelength(WaterRaman(), 440.0); digits = 1)
517.4
```
"""
raman_peak_wavelength(p::WaterRaman{FT}, λ_prime::Real) where {FT} =
    convert(_result_type(FT, λ_prime), 1e7 / (1e7 / _wide(λ_prime) - _wide(p.shift_cm_inv)))
