# =============================================================================
# Henyey–Greenstein phase function
# =============================================================================
#
# The Henyey–Greenstein form (Henyey & Greenstein 1941, Astrophys. J. 93, 70)
# is an analytic one-parameter approximation with especially simple
# Legendre and backscatter properties:
#
#     p(μ) = (1 - g²) / (1 + g² - 2 g μ)^{3/2},    normalized so (1/2)∫p dμ = 1
#
# Properties
# ----------
#   * asymmetry parameter   g  = <μ>
#   * Legendre moments      β_ℓ = g^ℓ     (closed form, no quadrature)
#   * g = 0 → isotropic; g → +1 → forward peak; g → −1 → back peak
#
# HG is a poor fit for real ocean particles (a single Mie kernel is not
# well approximated by one parameter), but it remains useful:
#   * as a legacy / teaching phase function,
#   * as a smooth analytic test target for the Legendre machinery,
#   * as the kernel of the widely-used two-term TTHG form
#     `p_TT = α p_HG(g₁) + (1-α) p_HG(g₂)` for ocean particles
#     (Haltrin 2002, Appl. Opt. 41, 1022).
# =============================================================================

"""
$(TYPEDEF)

Single-term Henyey–Greenstein phase function with asymmetry parameter `g`.

# Fields
$(TYPEDFIELDS)
"""
struct HenyeyGreensteinPhase{FT} <: AbstractOceanPhaseFunction{FT}
    "Asymmetry parameter `g = ⟨μ⟩ = β₁`, must satisfy `|g| < 1`"
    g::FT

    function HenyeyGreensteinPhase{FT}(g::FT) where {FT<:AbstractFloat}
        abs(g) < one(FT) ||
            throw(ArgumentError("HenyeyGreenstein: |g| must be strictly less than 1, got g=$g"))
        new{FT}(g)
    end
end

HenyeyGreensteinPhase{FT}(; g = FT(0.85)) where {FT<:AbstractFloat} =
    HenyeyGreensteinPhase{FT}(FT(g))
HenyeyGreensteinPhase(; kwargs...) = HenyeyGreensteinPhase{Float64}(; kwargs...)

# -----------------------------------------------------------------------------
# Phase function value and moments
# -----------------------------------------------------------------------------

# Un-annotated in its second argument so `cosθ::Dual{...}` propagates.
# The denominator is bounded away from zero for `|g| < 1` and `μ ∈ [-1, 1]`
# (minimum at μ=1 gives (1-g)² > 0), so no epsilon floor is needed.
function phase_function_value(pf::HenyeyGreensteinPhase{FT}, cosθ::Real) where {FT}
    g = pf.g
    return (1 - g^2) / (1 + g^2 - 2 * g * cosθ)^FT(1.5)
end

# Closed-form moments: β_ℓ = g^ℓ. Overrides the numerical-quadrature
# fallback in phase/moments.jl, which would otherwise recompute them at
# quadrature cost O(ℓ_max²) per call.
function phase_function_moments(pf::HenyeyGreensteinPhase{FT}, ℓ_max::Integer) where {FT}
    ℓ_max ≥ 0 || throw(ArgumentError("ℓ_max must be non-negative, got $ℓ_max"))
    β = Vector{FT}(undef, ℓ_max + 1)
    g = pf.g
    @inbounds for ℓ in 0:ℓ_max
        β[ℓ + 1] = g^ℓ
    end
    return β
end

# -----------------------------------------------------------------------------
# Closed-form summary diagnostics
# -----------------------------------------------------------------------------

asymmetry_parameter(pf::HenyeyGreensteinPhase) = pf.g

"""
    backscatter_fraction(pf::HenyeyGreensteinPhase) -> FT

Closed form derived from direct integration of `(1/2) ∫₋₁⁰ p(μ) dμ` with
substitution `u = 1 + g² - 2gμ`:

```
B(g) = (1-g²)/(2g) · [ 1/√(1+g²) - 1/(1+g) ]
```

At `g = 0` the expression is 0/0; we short-circuit to `B = 1/2` (isotropic).
"""
function backscatter_fraction(pf::HenyeyGreensteinPhase{FT}) where {FT}
    g = pf.g
    # Isotropic-limit short-circuit. Comparing to a tight epsilon of the
    # float type keeps this branch inert for any `g` a user can realistically
    # construct (the |g| < 1 constructor check excludes the other boundary).
    abs(g) < 4 * eps(FT) && return FT(0.5)
    return (1 - g^2) / (2 * g) * (1 / sqrt(1 + g^2) - 1 / (1 + g))
end

is_polarizable(::HenyeyGreensteinPhase) = false
