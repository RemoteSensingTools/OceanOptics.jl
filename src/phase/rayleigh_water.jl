# =============================================================================
# Depolarized Rayleigh phase function for pure water
# =============================================================================
#
# Pure water molecules scatter light nearly — but not exactly — as Rayleigh
# scatterers. The depolarization ratio ρ ≠ 0 accounts for anisotropic
# polarizability.
#
# The classical Morel (1974) value is ρ = 0.09, but more recent work
# (Farinato & Rowell 1976, and the Zhang et al. 2009 line) gives
# ρ ≈ 0.039 for pure seawater at 500 nm. We use ρ = 0.039 as the default;
# the user can override.
#
# The phase function in the intensity-only (scalar) form is:
#
#     p(μ) = 3 / [2 (2 + σ)] · [(1 + σ) + (1 - σ) μ²]
#
# where σ = 2ρ/(1 - ρ) is the Kahn & Arking "anisotropy factor". The
# normalization (1/2) ∫₋₁¹ p dμ = 1 is exact.
#
# Legendre moments are closed-form: β₀ = 1, β₂ = (1-σ)/(10·(...)),
# all other β_ℓ = 0 for ℓ ≠ 0, 2.
#
# For polarized RT, the full Greek-coefficient expansion (α, γ, δ, ε, ζ)
# has a closed form in terms of σ (Hansen & Travis 1974). A future
# `polarized_moments(pf) -> GreekCoefs` method can extend this.
# =============================================================================

"""
$(TYPEDEF)

Rayleigh-type phase function for pure water, with depolarization.

# Fields
$(TYPEDFIELDS)

Default `depolarization = 0.039` follows Farinato & Rowell (1976) and is
consistent with the Zhang, Hu & He (2009) pure-seawater model. The
Morel (1974) older value `ρ = 0.09` can be selected explicitly:

```julia
RayleighWaterPhase{Float64}(depolarization = 0.09)
```
"""
struct RayleighWaterPhase{FT} <: AbstractOceanPhaseFunction{FT}
    "Depolarization ratio ρ = I_⟂ / I_∥ at 90°"
    depolarization::FT

    function RayleighWaterPhase{FT}(depolarization::FT) where {FT<:AbstractFloat}
        (zero(FT) ≤ depolarization < FT(1)) ||
            throw(ArgumentError("depolarization must be in [0, 1), got $depolarization"))
        new{FT}(depolarization)
    end
end

RayleighWaterPhase{FT}(; depolarization = FT(0.039)) where {FT<:AbstractFloat} =
    RayleighWaterPhase{FT}(FT(depolarization))
RayleighWaterPhase(; kwargs...) = RayleighWaterPhase{Float64}(; kwargs...)

# "Anisotropy factor" σ = 2ρ/(1-ρ) — shows up everywhere
_rayleigh_σ(pf::RayleighWaterPhase{FT}) where {FT} = FT(2) * pf.depolarization / (one(FT) - pf.depolarization)

# -----------------------------------------------------------------------------

function phase_function_value(pf::RayleighWaterPhase{FT}, cosθ::Real) where {FT}
    μ = FT(cosθ)
    σ = _rayleigh_σ(pf)
    # Normalized so (1/2) ∫₋₁¹ p dμ = 1
    return FT(3) / (FT(2) * (FT(2) + σ)) * ((one(FT) + σ) + (one(FT) - σ) * μ^2)
end

function phase_function_moments(pf::RayleighWaterPhase{FT}, ℓ_max::Integer) where {FT}
    β = zeros(FT, ℓ_max + 1)
    σ = _rayleigh_σ(pf)
    β[1] = one(FT)                                    # β₀ = 1 (normalization)
    # The Rayleigh-with-depolarization phase function
    #   p(μ) = 3/(2(2+σ)) · [(1+σ) + (1-σ)μ²]
    # is a degree-2 polynomial in μ, so only β₀ and β₂ are non-zero in the
    # convention p(μ) = Σ_ℓ (2ℓ+1) β_ℓ P_ℓ(μ), β_ℓ = (1/2)∫₋₁¹ p·P_ℓ dμ.
    # Evaluating the integral directly yields β₂ = (1-σ)/(5(2+σ)); at σ=0
    # this recovers the pure-Rayleigh value β₂ = 1/10.
    if ℓ_max ≥ 2
        β[3] = (one(FT) - σ) / (FT(5) * (FT(2) + σ))
    end
    return β
end

function backscatter_fraction(pf::RayleighWaterPhase{FT}) where {FT}
    # Rayleigh is forward/back symmetric, so B = 0.5 exactly for any σ.
    return FT(0.5)
end

asymmetry_parameter(::RayleighWaterPhase{FT}) where {FT} = zero(FT)

is_polarizable(::RayleighWaterPhase) = true   # full Greek expansion is known; see notes above
