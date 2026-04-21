# =============================================================================
# Water Raman scattering as an inelastic process
# =============================================================================
#
# In the inelastic-source formalism, water Raman is the one "scattering-
# that-shifts-wavelength" process in pure water. It has:
#
#   * a fixed Stokes shift `Δν ≈ 3400 cm⁻¹` between excitation and
#     emission wavenumbers (the symmetric O-H stretch of liquid water),
#   * a Gaussian-like emission band of width ~200 cm⁻¹ FWHM (σ_ν ≈ 85
#     cm⁻¹) at a given excitation,
#   * a wavelength-dependent scattering cross-section
#     `b_w^R(λ') = b_R(λ_ref) · (λ_ref / λ')^n` with `n ≈ 5.5` per
#     Haltrin & Kattawar (1991, Appl. Opt. 30, 630; 1993, Appl. Opt. 32,
#     5356) and Bartlett et al. (1998, JGR 103, 17875),
#   * a Cabannes-depolarized phase function (anisotropic Rayleigh-like).
#
# In Fell's factorization the Raman cross-section plays the role of
# `a^I(λ')` in `b^I = a^I · f^I`, i.e. the "absorption-equivalent" that
# removes photons from the excitation field and deposits them in the
# emission band.
# =============================================================================

"""
$(TYPEDEF)

Water Raman scattering as an inelastic process. Emission spectrum is a
Gaussian in wavenumber centred at `ν' − shift_cm_inv`; cross-section
follows a Haltrin-Kattawar power law.

# Fields
$(TYPEDFIELDS)

# Construction

```julia
# Default: Haltrin-Kattawar/Bartlett parameterization
WaterRaman{Float64}()

# Customize (e.g. with Bartlett salt- and temperature-aware values)
WaterRaman{Float64}(
    shift_cm_inv       = 3400.0,
    σ_cm_inv           = 85.0,
    cross_section_ref  = 2.7e-4,    # 1/m at reference wavelength
    cross_section_λ_ref = 488.0,    # nm
    cross_section_n     = 5.5,
    depolarization      = 0.14,     # Cabannes ρ for liquid water
    excitation_range    = (350.0, 700.0),
)
```
"""
struct WaterRaman{FT<:AbstractFloat} <: AbstractOceanInelasticProcess{FT}
    "Stokes shift between excitation and emission wavenumbers `[1/cm]`"
    shift_cm_inv::FT
    "Gaussian standard deviation in wavenumber space `[1/cm]`"
    σ_cm_inv::FT
    "Raman cross-section `b_w^R` at the reference wavelength `[1/m]`"
    cross_section_ref::FT
    "Reference wavelength for the cross-section power law `[nm]`"
    cross_section_λ_ref::FT
    "Exponent `n` in `b_w^R(λ') = b_ref · (λ_ref / λ')^n`"
    cross_section_n::FT
    "Cabannes depolarization ratio ρ_Cab (≈ 0.14 for liquid water)"
    depolarization::FT
    "Excitation-wavelength integration window `[nm]`"
    excitation_range::Tuple{FT, FT}
end

"""
    WaterRaman{FT}(; shift_cm_inv=3400, σ_cm_inv=85,
                     cross_section_ref=2.7e-4, cross_section_λ_ref=488,
                     cross_section_n=5.5, depolarization=0.14,
                     excitation_range=(350, 700))
"""
function WaterRaman{FT}(;
        shift_cm_inv                            = FT(3400.0),
        σ_cm_inv                                = FT(85.0),
        cross_section_ref                       = FT(2.7e-4),
        cross_section_λ_ref                     = FT(488.0),
        cross_section_n                         = FT(5.5),
        depolarization                          = FT(0.14),
        excitation_range::Tuple{<:Real, <:Real} = (FT(350), FT(700)),
) where {FT<:AbstractFloat}
    FT(σ_cm_inv) > zero(FT) ||
        throw(ArgumentError("WaterRaman: σ_cm_inv must be positive, got $σ_cm_inv"))
    FT(cross_section_ref) ≥ zero(FT) ||
        throw(ArgumentError("WaterRaman: cross_section_ref must be non-negative, got $cross_section_ref"))
    zero(FT) ≤ FT(depolarization) < one(FT) ||
        throw(ArgumentError("WaterRaman: depolarization must be in [0, 1), got $depolarization"))
    excitation_range[1] < excitation_range[2] ||
        throw(ArgumentError("WaterRaman: excitation_range must be (lo < hi)"))
    return WaterRaman{FT}(
        FT(shift_cm_inv),
        FT(σ_cm_inv),
        FT(cross_section_ref),
        FT(cross_section_λ_ref),
        FT(cross_section_n),
        FT(depolarization),
        (FT(excitation_range[1]), FT(excitation_range[2])))
end

WaterRaman(; kwargs...) = WaterRaman{Float64}(; kwargs...)

# -----------------------------------------------------------------------------
# Trait implementations
# -----------------------------------------------------------------------------

"""
    excitation_absorption(p::WaterRaman, λ_prime) -> FT

Raman "absorption-role" cross-section at excitation wavelength `λ'`,
using the Haltrin-Kattawar power law

```
b_w^R(λ') = b_ref · (λ_ref / λ')^n
```

in units of `[1/m]`.
"""
excitation_absorption(p::WaterRaman, λ_prime::Real) =
    p.cross_section_ref * (p.cross_section_λ_ref / λ_prime)^p.cross_section_n

"""
    emission(p::WaterRaman, λ_prime, λ) -> FT

Gaussian redistribution in wavenumber space centered at `ν' − shift`,
evaluated against the Jacobian `|dν/dλ| = 10⁷/λ²` so that
`∫ f(λ', λ) dλ = 1`.
"""
emission(p::WaterRaman{FT}, λ_prime::Real, λ::Real) where {FT} =
    redistribution(WavenumberShiftEmission{FT}(p.shift_cm_inv, p.σ_cm_inv),
                   λ_prime, λ)

excitation_range(p::WaterRaman) = p.excitation_range

"""
    is_isotropic(p::WaterRaman) -> Bool

Water Raman is Cabannes-depolarized: the re-emitted angular
distribution is Rayleigh-like with depolarization `ρ_Cab ≈ 0.14`.
Strictly, that is not isotropic, so we report `false`. Solvers
targeting scalar RT can still project onto the `m = 0` mode with a
small (sub-percent) error.
"""
is_isotropic(p::WaterRaman) = false

function Base.show(io::IO, p::WaterRaman{FT}) where {FT}
    print(io, "WaterRaman{$FT}(Δν=", p.shift_cm_inv, " 1/cm, σ=",
              p.σ_cm_inv, " 1/cm, b_R(",
              Int(round(p.cross_section_λ_ref)), ")=",
              p.cross_section_ref, " 1/m)")
end

"""
    raman_peak_wavelength(p::WaterRaman, λ_prime::Real) -> FT

Emission-peak wavelength `[nm]` for excitation at `λ_prime`:

```
ν_emission = ν' - Δν,   λ_peak = 10⁷ / ν_emission
```

Diagnostic helper — for `λ' = 440 nm` and the default `Δν = 3400 cm⁻¹`
this returns ≈ 491 nm.
"""
raman_peak_wavelength(p::WaterRaman, λ_prime::Real) =
    1e7 / (1e7 / λ_prime - p.shift_cm_inv)
