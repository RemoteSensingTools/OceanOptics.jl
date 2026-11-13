# =============================================================================
# Non-Algal Particles (NAP)
# =============================================================================
#
# Non-algal particulate matter — the union of detrital organic material
# (cell fragments, faecal pellets) and mineral/inorganic particles
# (sediments, aeolian dust, resuspended bottom material) — is the second
# major particulate-absorber and -scatterer pool in natural waters after
# phytoplankton. Two parameterizations are shipped:
#
# BabinBricaud2003  — detritus + minerals combined, spectral exponential
#                     absorption à la CDOM with a distinct shallower slope
#                     (Babin et al. 2003, JGR 108, 3211):
#
#                     a_nap(λ) = a_nap(λ_ref) · exp(-S · (λ - λ_ref))
#
#                     Typical S ≈ 0.011 1/nm (compare CDOM S ≈ 0.014).
#
# FellMineral       — beam attenuation of inorganic suspended matter
#                     (Fell 1997 §5.1.4, Eqs. 5.10-5.11):  c(λ) = c_ref · (440/λ)
#                     with a fixed single-scattering albedo ω₀ = 0.95, so that
#                     a = (1 − ω₀)·c and b = ω₀·c, and a backscatter
#                     probability of 0.036.
#
# Both parameterizations scatter into a Fournier-Forand phase function
# with a larger backscatter fraction than living phytoplankton (default
# B = 0.02 for BabinBricaud2003, 0.036 for FellMineral): Mobley, Sundman &
# Boss (2002) recommend FF over Petzold. Fell (1997) used a Petzold San Diego
# Harbor phase function modified to the same backscatter probability.
# =============================================================================

"""
$(TYPEDEF)

Babin et al. (2003) detritus + minerals absorption model,
`a = a_ref exp(−S (λ − λ_ref))`, with the particulate scattering law
`b = b_ref (λ_ref/λ)^γ` (`γ = spectral_exponent`, default 1).
"""
struct BabinBricaud2003 <: AbstractNAPModel end

"""
$(TYPEDEF)

Fell (1997) §5.1.4 beam-attenuation model of inorganic suspended matter,
`c = c_ref · (λ_ref/λ)^γ` with `λ_ref = 440` nm and `γ = spectral_exponent = 1`
(Eq. 5.10), split into absorption and scattering by a fixed single-scattering albedo.

The defaults are the thesis values: single-scattering albedo `ω_nap = 0.95`
(Eq. 5.11) and backscatter fraction `B = 0.036` (p. 120). The phase function is
Fournier-Forand with that backscatter fraction; Fell used a Petzold San Diego
Harbor phase function modified to it.
"""
struct FellMineral <: AbstractNAPModel end

"""
$(TYPEDEF)

Non-algal particulate matter (detritus + minerals) as an absorbing
scatterer.

# Fields
$(TYPEDFIELDS)

# Examples

```julia
# Babin et al. (2003) default: a_nap(443) = 0.05 m⁻¹, S = 0.011 1/nm
nap = NonAlgalParticles{Float64}()

# Mineral-dominated case per Fell (1997)
nap = NonAlgalParticles{Float64, FellMineral}(
    c_ref  = 0.40,   # beam attenuation at 440 nm
    ω_nap  = 0.95,   # single-scattering albedo (the default)
)
```
"""
struct NonAlgalParticles{FT, M<:AbstractNAPModel,
                            P<:AbstractOceanPhaseFunction{FT}} <: AbstractAbsorbingScatterer{FT}
    "Parameterization-model tag (e.g. `BabinBricaud2003()` or `FellMineral()`)"
    model::M
    "Reference absorption `a_nap(λ_ref)` `[1/m]` (BabinBricaud2003)"
    a_ref::FT
    "Reference scattering `b_nap(λ_ref)` `[1/m]` (BabinBricaud2003)"
    b_ref::FT
    "Reference wavelength `[nm]`"
    λ_ref::FT
    "Exponential absorption slope `[1/nm]` (BabinBricaud2003)"
    slope::FT
    "Beam attenuation at 440 nm `[1/m]` (FellMineral)"
    c_ref::FT
    "Single-scattering albedo at 440 nm (FellMineral)"
    ω_nap::FT
    "Exponent `γ` of the `(λ_ref/λ)^γ` law of scattering (BabinBricaud2003) or attenuation (FellMineral)"
    spectral_exponent::FT
    "Particulate-scattering phase function"
    phase::P
end

# BabinBricaud2003 constructor --------------------------------------------------

function NonAlgalParticles{FT, BabinBricaud2003}(;
        a_ref  = 0.05,
        b_ref  = 0.45,
        λ_ref  = 443,
        slope  = 0.011,
        B      = 0.02,
        spectral_exponent = 1,
) where {FT<:Real}
    _check_parameter_type(FT)
    a_ref ≥ zero(a_ref) ||
        throw(ArgumentError("NAP: a_ref must be non-negative, got $a_ref"))
    b_ref ≥ zero(b_ref) ||
        throw(ArgumentError("NAP: b_ref must be non-negative, got $b_ref"))
    _within(slope, 0.005, 0.020) ||
        @warn "NAP: BabinBricaud2003 slope $slope 1/nm outside typical 0.005-0.020 range"
    phase = FournierForandPhase{FT}(backscatter_fraction = FT(B))
    return NonAlgalParticles{FT, BabinBricaud2003, typeof(phase)}(
        BabinBricaud2003(), FT(a_ref), FT(b_ref), FT(λ_ref), FT(slope),
        zero(FT), zero(FT), FT(spectral_exponent), phase)
end

# FellMineral constructor -------------------------------------------------------

function NonAlgalParticles{FT, FellMineral}(;
        c_ref  = 0.40,
        ω_nap  = 0.95,
        λ_ref  = 440,
        B      = 0.036,
        spectral_exponent = 1,
) where {FT<:Real}
    _check_parameter_type(FT)
    c_ref ≥ zero(c_ref) ||
        throw(ArgumentError("NAP: c_ref must be non-negative, got $c_ref"))
    _within(ω_nap, 0, 1) ||
        throw(ArgumentError("NAP: single-scattering albedo ω_nap must be in [0, 1], got $ω_nap"))
    phase = FournierForandPhase{FT}(backscatter_fraction = FT(B))
    return NonAlgalParticles{FT, FellMineral, typeof(phase)}(
        FellMineral(), zero(FT), zero(FT), FT(λ_ref), zero(FT),
        FT(c_ref), FT(ω_nap), FT(spectral_exponent), phase)
end

# Default model tag is BabinBricaud2003
NonAlgalParticles{FT}(; kwargs...) where {FT<:Real} =
    NonAlgalParticles{FT, BabinBricaud2003}(; kwargs...)
NonAlgalParticles(; kwargs...) = NonAlgalParticles{_parameter_type(values(kwargs))}(; kwargs...)

scattering_phase(nap::NonAlgalParticles) = nap.phase

# -----------------------------------------------------------------------------
# BabinBricaud2003 dispatch: exponential a, (440/λ) b
# -----------------------------------------------------------------------------

function absorption(nap::NonAlgalParticles{FT, BabinBricaud2003}, λ::Real) where {FT}
    return nap.a_ref * exp(-nap.slope * (λ - nap.λ_ref))
end

function scattering(nap::NonAlgalParticles{FT, BabinBricaud2003}, λ::Real) where {FT}
    return nap.b_ref * (nap.λ_ref / λ)^nap.spectral_exponent
end

# -----------------------------------------------------------------------------
# FellMineral dispatch: c = c_ref·(440/λ), split by ω_nap
# -----------------------------------------------------------------------------

function absorption(nap::NonAlgalParticles{FT, FellMineral}, λ::Real) where {FT}
    c = nap.c_ref * (nap.λ_ref / λ)^nap.spectral_exponent
    return (one(FT) - nap.ω_nap) * c
end

function scattering(nap::NonAlgalParticles{FT, FellMineral}, λ::Real) where {FT}
    c = nap.c_ref * (nap.λ_ref / λ)^nap.spectral_exponent
    return nap.ω_nap * c
end

# -----------------------------------------------------------------------------
# Shared vectorized forms + backscattering
# -----------------------------------------------------------------------------

absorption(nap::NonAlgalParticles, λ::AbstractVector) = map(Base.Fix1(absorption, nap), λ)
scattering(nap::NonAlgalParticles, λ::AbstractVector) = map(Base.Fix1(scattering, nap), λ)

backscattering(nap::NonAlgalParticles, λ::Real) =
    scattering(nap, λ) * backscatter_fraction(nap.phase)
backscattering(nap::NonAlgalParticles, λ::AbstractVector) =
    map(Base.Fix1(backscattering, nap), λ)

function Base.show(io::IO, nap::NonAlgalParticles{FT, M}) where {FT, M}
    print(io, "NonAlgalParticles{$FT}(", nameof(M), ", B=",
              backscatter_fraction(nap.phase), ")")
end
