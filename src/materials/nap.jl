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
# FellMineral       — crude beam-attenuation for mineral-dominated waters
#                     (Fell 1997 §5.1.4):  c(λ) = c_ref · (440/λ)
#                     Split into a = a_ref·(440/λ), b = (c_ref − a_ref)·(440/λ)
#                     using a fixed single-scattering albedo ϖ_NAP.
#
# Both parameterizations scatter into a Fournier-Forand phase function
# with a larger backscatter fraction than living phytoplankton (typical
# B ≈ 0.02 for detritus, 0.03 for minerals) — DESIGN §6.1 and Mobley
# (2002) justify FF over Petzold.
# =============================================================================

"""
$(TYPEDEF)

Babin et al. (2003) detritus + minerals absorption model with Morel &
Maritorena (2001)-style particulate scattering.
"""
struct BabinBricaud2003 <: AbstractNAPModel end

"""
$(TYPEDEF)

Fell (1997) §5.1.4 mineral-dominated beam-attenuation model `c = c_ref
· (440/λ)`, split into absorption and scattering by a fixed single-
scattering albedo.
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
    ω_nap  = 0.96,   # single-scattering albedo
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
    "Particulate-scattering phase function"
    phase::P
end

# BabinBricaud2003 constructor --------------------------------------------------

function NonAlgalParticles{FT, BabinBricaud2003}(;
        a_ref  = FT(0.05),
        b_ref  = FT(0.45),
        λ_ref  = FT(443),
        slope  = FT(0.011),
        B      = FT(0.02),
) where {FT<:AbstractFloat}
    a_ref ≥ zero(a_ref) ||
        throw(ArgumentError("NAP: a_ref must be non-negative, got $a_ref"))
    b_ref ≥ zero(b_ref) ||
        throw(ArgumentError("NAP: b_ref must be non-negative, got $b_ref"))
    (FT(0.005) ≤ slope ≤ FT(0.020)) ||
        @warn "NAP: BabinBricaud2003 slope $slope 1/nm outside typical 0.005-0.020 range"
    phase = FournierForandPhase{FT}(backscatter_fraction = FT(B))
    return NonAlgalParticles{FT, BabinBricaud2003, typeof(phase)}(
        BabinBricaud2003(), FT(a_ref), FT(b_ref), FT(λ_ref), FT(slope),
        zero(FT), zero(FT), phase)
end

# FellMineral constructor -------------------------------------------------------

function NonAlgalParticles{FT, FellMineral}(;
        c_ref  = FT(0.40),
        ω_nap  = FT(0.96),
        λ_ref  = FT(440),
        B      = FT(0.03),
) where {FT<:AbstractFloat}
    c_ref ≥ zero(c_ref) ||
        throw(ArgumentError("NAP: c_ref must be non-negative, got $c_ref"))
    (zero(ω_nap) ≤ ω_nap ≤ one(ω_nap)) ||
        throw(ArgumentError("NAP: single-scattering albedo ω_nap must be in [0, 1], got $ω_nap"))
    phase = FournierForandPhase{FT}(backscatter_fraction = FT(B))
    return NonAlgalParticles{FT, FellMineral, typeof(phase)}(
        FellMineral(), zero(FT), zero(FT), FT(λ_ref), zero(FT),
        FT(c_ref), FT(ω_nap), phase)
end

# Default model tag is BabinBricaud2003
NonAlgalParticles{FT}(; kwargs...) where {FT<:AbstractFloat} =
    NonAlgalParticles{FT, BabinBricaud2003}(; kwargs...)
NonAlgalParticles(; kwargs...) = NonAlgalParticles{Float64}(; kwargs...)

phase_function(nap::NonAlgalParticles) = nap.phase

# -----------------------------------------------------------------------------
# BabinBricaud2003 dispatch: exponential a, (440/λ) b
# -----------------------------------------------------------------------------

function absorption(nap::NonAlgalParticles{FT, BabinBricaud2003}, λ::Real) where {FT}
    return nap.a_ref * exp(-nap.slope * (λ - nap.λ_ref))
end

function scattering(nap::NonAlgalParticles{FT, BabinBricaud2003}, λ::Real) where {FT}
    return nap.b_ref * (nap.λ_ref / λ)
end

# -----------------------------------------------------------------------------
# FellMineral dispatch: c = c_ref·(440/λ), split by ω_nap
# -----------------------------------------------------------------------------

function absorption(nap::NonAlgalParticles{FT, FellMineral}, λ::Real) where {FT}
    c = nap.c_ref * (nap.λ_ref / λ)
    return (one(FT) - nap.ω_nap) * c
end

function scattering(nap::NonAlgalParticles{FT, FellMineral}, λ::Real) where {FT}
    c = nap.c_ref * (nap.λ_ref / λ)
    return nap.ω_nap * c
end

# -----------------------------------------------------------------------------
# Shared vectorized forms + backscattering
# -----------------------------------------------------------------------------

absorption(nap::NonAlgalParticles, λ::AbstractVector) = [absorption(nap, λi) for λi in λ]
scattering(nap::NonAlgalParticles, λ::AbstractVector) = [scattering(nap, λi) for λi in λ]

backscattering(nap::NonAlgalParticles, λ::Real) =
    scattering(nap, λ) * backscatter_fraction(nap.phase)
backscattering(nap::NonAlgalParticles, λ::AbstractVector) =
    [backscattering(nap, λi) for λi in λ]

function Base.show(io::IO, nap::NonAlgalParticles{FT, M}) where {FT, M}
    print(io, "NonAlgalParticles{$FT}(", nameof(M), ", B=",
              backscatter_fraction(nap.phase), ")")
end
