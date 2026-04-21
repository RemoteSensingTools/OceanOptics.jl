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
# For polarized RT, the full Greek-coefficient expansion (α, γ, δ, ϵ, ζ)
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


# -----------------------------------------------------------------------------

# Proutière depolarization factors used throughout the Greek expansion
# (Sanghavi 2014 Eq. 16; de Haan 1987; cf. vSmartMOM get_greek_rayleigh).
_rayleigh_dpl_p(ρ) = (1 - ρ) / (1 + ρ / 2)
_rayleigh_dpl_r(ρ) = (1 - 2ρ) / (1 - ρ)

function phase_function_value(pf::RayleighWaterPhase{FT}, cosθ::Real) where {FT}
    # De Haan convention: p(μ) = Σ β_ℓ P_ℓ(μ). For Rayleigh only β_0 = 1
    # and β_2 = ½·dpl_p are non-zero, giving
    #   p(μ) = 1 + β₂ · (3μ² − 1)/2
    μ      = FT(cosθ)
    β₂     = FT(0.5) * _rayleigh_dpl_p(pf.depolarization)
    return one(FT) + β₂ * (FT(3) * μ^2 - one(FT)) / FT(2)
end

function phase_function_moments(pf::RayleighWaterPhase{FT}, ℓ_max::Integer) where {FT}
    # Rayleigh is a degree-2 polynomial in μ, so only β₀ and β₂ are non-zero
    # in the de Haan convention (`p(μ) = Σ β_ℓ P_ℓ(μ)`). The closed-form
    # values come from Hansen-Travis (1974) with Proutière's depolarization
    # factor dpl_p = (1 − ρ)/(1 + ρ/2); at ρ = 0 this recovers pure
    # Rayleigh with β₂ = ½.
    β = zeros(FT, ℓ_max + 1)
    β[1] = one(FT)
    if ℓ_max ≥ 2
        β[3] = FT(0.5) * _rayleigh_dpl_p(pf.depolarization)
    end
    return β
end

"""
    phase_matrix_moments(pf::RayleighWaterPhase, ℓ_max) -> NamedTuple

Closed-form polarized Greek expansion of the depolarized Rayleigh phase
matrix (Hansen & Travis 1974; Sanghavi 2014 Eq. 16). Only ℓ ∈ {0, 1, 2}
carries non-zero moments:

```
β[0] = 1                     α[2] = 3·dpl_p
β[2] = ½·dpl_p               γ[2] = √(3/2)·dpl_p
δ[1] = (3/2)·dpl_p·dpl_r     ϵ[ℓ] = ζ[ℓ] = 0

dpl_p = (1 − ρ)/(1 + ρ/2)
dpl_r = (1 − 2ρ)/(1 − ρ)
```

with `ρ = pf.depolarization` the Cabannes depolarization ratio.
"""
function phase_matrix_moments(pf::RayleighWaterPhase{FT}, ℓ_max::Integer) where {FT}
    ℓ_max ≥ 0 || throw(ArgumentError("ℓ_max must be non-negative, got $ℓ_max"))
    N = ℓ_max + 1
    α = zeros(FT, N); β = zeros(FT, N); γ = zeros(FT, N)
    δ = zeros(FT, N); ϵ = zeros(FT, N); ζ = zeros(FT, N)

    ρ     = pf.depolarization
    dpl_p = _rayleigh_dpl_p(ρ)
    dpl_r = _rayleigh_dpl_r(ρ)

    β[1] = one(FT)                                 # ℓ = 0
    if ℓ_max ≥ 1
        δ[2] = FT(1.5) * dpl_p * dpl_r             # ℓ = 1
    end
    if ℓ_max ≥ 2
        α[3] = FT(3) * dpl_p                       # ℓ = 2
        β[3] = FT(0.5) * dpl_p
        γ[3] = dpl_p * sqrt(FT(1.5))
    end
    return (α = α, β = β, γ = γ, δ = δ, ϵ = ϵ, ζ = ζ)
end

function backscatter_fraction(pf::RayleighWaterPhase{FT}) where {FT}
    # Rayleigh is forward/back symmetric, so B = 0.5 exactly for any ρ.
    return FT(0.5)
end

asymmetry_parameter(::RayleighWaterPhase{FT}) where {FT} = zero(FT)

is_polarizable(::RayleighWaterPhase) = true
