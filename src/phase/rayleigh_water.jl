# =============================================================================
# Depolarized Rayleigh phase function for pure water
# =============================================================================
#
# Molecular scattering by water is Rayleigh scattering with a depolarization
# ratio ρ > 0 (anisotropic polarizability); ρ is the ratio for unpolarized
# incident light observed at 90°.
#
# Scalar phase function, normalized to (1/2) ∫₋₁¹ p dμ = 1:
#
#     p(μ) = 1 + β₂ P₂(μ),   β₂ = Δ/2,   Δ = (1 − ρ)/(1 + ρ/2),
#     P₂(μ) = (3μ² − 1)/2
#
# which is the Cabannes form (3/(4(1+2γ))) [(1+3γ) + (1−γ) μ²] with
# γ = ρ/(2−ρ). All other Legendre moments vanish.
#
# The phase matrix has the closed-form Greek expansion of Hansen & Travis
# (1974), available through `greek_expansion(pf, ℓ_max; polarized = true)`.
# =============================================================================

"Depolarization ratio of pure water, Farinato & Rowell (1976), as used by Zhang, Hu & He (2009)."
const _WATER_DEPOLARIZATION = 0.039

"""
$(TYPEDEF)

Rayleigh-type phase function for pure water, with depolarization.

# Fields
$(TYPEDFIELDS)

The phase function is `p(μ) = 1 + β₂ P₂(μ)` with `β₂ = (1 − ρ)/(2 + ρ)`, `ρ` the
depolarization ratio for unpolarized incident light at 90°. The default `ρ = 0.039`
follows Farinato & Rowell (1976) and is the value of the Zhang, Hu & He (2009)
pure-seawater model. The older value `ρ = 0.09` (Morel 1974) can be selected
explicitly:

```julia
RayleighWaterPhase{Float64}(depolarization = 0.09)
```
"""
struct RayleighWaterPhase{FT} <: AbstractOceanPhaseFunction{FT}
    "Depolarization ratio ρ = I_⟂ / I_∥ at 90° for unpolarized incident light"
    depolarization::FT

    function RayleighWaterPhase{FT}(depolarization::FT) where {FT<:Real}
        _check_parameter_type(FT)
        (_within(depolarization, 0, 1) && depolarization < FT(1)) ||
            throw(ArgumentError("depolarization must be in [0, 1), got $depolarization"))
        new{FT}(depolarization)
    end
end

RayleighWaterPhase{FT}(; depolarization = _WATER_DEPOLARIZATION) where {FT<:Real} =
    (_check_parameter_type(FT); RayleighWaterPhase{FT}(FT(depolarization)))
RayleighWaterPhase(; kwargs...) = RayleighWaterPhase{_parameter_type(values(kwargs))}(; kwargs...)


# -----------------------------------------------------------------------------

# Proutière depolarization factors used throughout the Greek expansion
# (Sanghavi 2014 Eq. 16; de Haan 1987; cf. vSmartMOM get_greek_rayleigh).
_rayleigh_dpl_p(ρ) = (1 - ρ) / (1 + ρ / 2)
_rayleigh_dpl_r(ρ) = (1 - 2ρ) / (1 - ρ)

function phase_function_value(pf::RayleighWaterPhase{FT}, cosθ::Real) where {FT}
    # De Haan convention: p(μ) = Σ β_ℓ P_ℓ(μ). For Rayleigh only β_0 = 1
    # and β_2 = ½·dpl_p are non-zero, giving
    #   p(μ) = 1 + β₂ · (3μ² − 1)/2
    β₂     = FT(0.5) * _rayleigh_dpl_p(pf.depolarization)
    return one(FT) + β₂ * (FT(3) * cosθ^2 - one(FT)) / FT(2)
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

# Closed-form polarized Greek expansion of the depolarized Rayleigh phase
# matrix (Hansen & Travis 1974; Sanghavi 2014 Eq. 16), identical to
# vSmartMOM's `get_greek_rayleigh`. Only ℓ ∈ {0, 1, 2} are non-zero:
#
#   β₀ = 1,  β₂ = ½·dpl_p,  α₂ = 3·dpl_p,  γ₂ = √(3/2)·dpl_p,
#   δ₁ = (3/2)·dpl_p·dpl_r,  ϵ = ζ = 0,
#   dpl_p = (1 − ρ)/(1 + ρ/2),  dpl_r = (1 − 2ρ)/(1 − ρ).
function _polarized_greek(pf::RayleighWaterPhase{FT}, ℓ_max::Integer) where {FT}
    ℓ_max ≥ 0 || throw(ArgumentError("ℓ_max must be non-negative, got $ℓ_max"))
    N = ℓ_max + 1
    α = zeros(FT, N); β = zeros(FT, N); γ = zeros(FT, N)
    δ = zeros(FT, N); ϵ = zeros(FT, N); ζ = zeros(FT, N)
    dpl_p = _rayleigh_dpl_p(pf.depolarization)
    dpl_r = _rayleigh_dpl_r(pf.depolarization)
    β[1] = one(FT)
    ℓ_max ≥ 1 && (δ[2] = FT(1.5) * dpl_p * dpl_r)
    if ℓ_max ≥ 2
        α[3] = FT(3) * dpl_p
        β[3] = FT(0.5) * dpl_p
        γ[3] = dpl_p * sqrt(FT(1.5))
    end
    return GreekCoefficients{FT}(α, β, γ, δ, ϵ, ζ)
end

function backscatter_fraction(pf::RayleighWaterPhase{FT}) where {FT}
    # Rayleigh is forward/back symmetric, so B = 0.5 exactly for any ρ.
    return FT(0.5)
end

asymmetry_parameter(::RayleighWaterPhase{FT}) where {FT} = zero(FT)

is_polarizable(::RayleighWaterPhase) = true
