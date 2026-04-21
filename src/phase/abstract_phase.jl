# =============================================================================
# Phase functions for ocean scattering
# =============================================================================
#
# Design: a phase function is a complete, self-contained object that knows
# how to evaluate itself at any scattering angle and how to produce its
# Legendre-moment expansion. The mathematics doesn't care whether the
# scatterer is water, phytoplankton, or a mineral particle; what matters is
# the shape of p(μ).
#
# Three standard forms are shipped:
#
#   RayleighWaterPhase           — molecular scattering with depolarization
#   FournierForandPhase          — the modern standard for ocean particles
#   PetzoldPhase                 — Petzold (1972) "San Diego Harbor" tabulated
#   HenyeyGreensteinPhase        — analytical, legacy; supports two-term form
#
# Anything else can be added by subtyping AbstractOceanPhaseFunction and
# implementing `phase_function_value` (single-angle) OR
# `phase_function_moments` (Legendre expansion) — the other is derived by
# numerical integration in `moments.jl`.
# =============================================================================

"""
$(TYPEDEF)

Abstract base type for scattering phase functions used in the ocean.

A phase function `p(μ)` describes the angular distribution of scattered
radiation, with `μ = cos(θ)` the cosine of the scattering angle.
Conventional normalization:

```
(1/2) ∫₋₁¹ p(μ) dμ = 1
```

The Legendre expansion convention in this package matches de Haan,
Bosma & Hovenier (1987) A&A 183, 371 and Sanghavi (2014) — the same
convention as HydroLight and vSmartMOM.jl:

```
p(μ) = Σ_ℓ β_ℓ P_ℓ(μ),   β_ℓ = (2ℓ+1)/2 · ∫₋₁¹ p(μ) P_ℓ(μ) dμ
```

This gives `β_0 = 1` for a normalized phase function and, for example,
`β_2 = 1/2` for pure Rayleigh (compare Hansen & Travis 1974).

# Required interface

A concrete subtype must implement *at least one* of:

- `phase_function_value(pf, μ) -> FT` — evaluate `p(μ)` at a single angle
- `phase_function_moments(pf, ℓ_max) -> Vector{FT}` — Legendre expansion
  coefficients `β_ℓ` of length `ℓ_max + 1` (index 1 = ℓ=0)

If only one is implemented, the other is obtained by numerical integration.
For expensive closed-form expansions, overloading both avoids redundant work.

# Optional interface

- `backscatter_fraction(pf) -> FT` — `B = (1/2)∫₋₁⁰ p(μ) dμ`. A default
  numerical-integration implementation is provided; closed-form phase
  functions (Fournier-Forand, Henyey-Greenstein) should override it.
- `asymmetry_parameter(pf) -> FT` — `g = (1/2)∫₋₁¹ μ p(μ) dμ = β₁ / 3`
- `is_polarizable(pf) -> Bool` — whether the full polarized Greek
  expansion `(α, β, γ, δ, ϵ, ζ)` is available via
  [`phase_matrix_moments`](@ref). Default `false`.
- `phase_matrix_moments(pf, ℓ_max) -> NamedTuple{(:α, :β, :γ, :δ, :ϵ, :ζ)}`
  — full polarized expansion. Default: β from `phase_function_moments`
  and zeros elsewhere (scalar fallback).
"""
abstract type AbstractOceanPhaseFunction{FT<:AbstractFloat} end

"""
    phase_function_value(pf, μ) -> FT

Evaluate the phase function at `μ = cos(θ)`.

The default implementation reconstructs the value from the Legendre
expansion — inefficient but correct for any subtype that implements
`phase_function_moments`.
"""
function phase_function_value(pf::AbstractOceanPhaseFunction{FT}, μ::Real) where {FT}
    # Default: reconstruct from Legendre moments, ℓ_max = 128 for accuracy.
    # Subtypes with closed forms should override.
    β = phase_function_moments(pf, 128)
    return _reconstruct_from_moments(β, FT(μ))
end

"""
    phase_function_moments(pf, ℓ_max) -> Vector{FT}

Return the Legendre expansion coefficients `β_ℓ` for ℓ = 0, 1, …, ℓ_max
satisfying

```
p(μ) = Σₗ βₗ Pₗ(μ),   with β₀ = 1 for normalized phase functions
```

This is the de Haan 1987 / HydroLight convention (see
[`AbstractOceanPhaseFunction`](@ref)). The default implementation
computes the moments by Gauss-Legendre quadrature from
`phase_function_value`. Subtypes with closed-form moments —
Henyey-Greenstein has `βₗ = (2ℓ+1) · gᴸ`, Rayleigh with
depolarization has a two-term expansion — should override.
"""
function phase_function_moments end

"""
    phase_matrix_moments(pf, ℓ_max) -> NamedTuple

Return the polarized scattering-matrix Greek coefficients
`(α, β, γ, δ, ϵ, ζ)` as a `NamedTuple` of length-`ℓ_max+1` vectors.
The mapping into the 4×4 phase matrix `B` is

```
B = [β  γ   0   0
     γ  α   0   0
     0  0   ζ   ϵ
     0  0  -ϵ   δ]
```

(Sanghavi 2014 Eq. 16; de Haan 1987 §2). The default implementation
promotes the scalar `phase_function_moments` into `β` and zeros the
other five — the right behavior for phase functions that are scalar-
only (Fournier-Forand, Henyey-Greenstein, Petzold). Subtypes with
closed-form polarized expansions (e.g. `RayleighWaterPhase`) should
override.
"""
function phase_matrix_moments(pf::AbstractOceanPhaseFunction{FT},
                              ℓ_max::Integer) where {FT}
    β = phase_function_moments(pf, ℓ_max)
    N = length(β)
    z = zeros(FT, N)
    return (α = copy(z), β = β, γ = copy(z),
            δ = copy(z), ϵ = copy(z), ζ = copy(z))
end

"""
    backscatter_fraction(pf) -> FT

Fraction of scattered light going into the back hemisphere (μ ∈ [-1, 0]).

Default: numerical integration. Closed-form phase functions should override.
"""
function backscatter_fraction(pf::AbstractOceanPhaseFunction{FT}) where {FT}
    # (1/2) ∫₋₁⁰ p(μ) dμ via Gauss-Legendre on [-1, 0]
    nq = 64
    μ_nodes, w_nodes = gausslegendre(nq)
    # Map [-1,1] → [-1,0]: μ' = (μ-1)/2, weight scaled by 1/2
    B = zero(FT)
    @inbounds for i in 1:nq
        μ′ = FT((μ_nodes[i] - 1) / 2)
        w′ = FT(w_nodes[i] / 2)
        B += w′ * phase_function_value(pf, μ′)
    end
    return B / 2   # the outer 1/2 from the definition
end

"""
    asymmetry_parameter(pf) -> FT

Mean cosine of the scattering angle, `g = <μ> = β₁ / 3`. In the de Haan
convention `β_1 = 3g`, so the division by 3 recovers `g`.
"""
function asymmetry_parameter(pf::AbstractOceanPhaseFunction{FT}) where {FT}
    β = phase_function_moments(pf, 1)
    return β[2] / FT(3)
end

"""
    is_polarizable(pf) -> Bool

Whether the phase function provides the full polarized scattering-matrix
expansion (α, β, γ, δ, ϵ, ζ). Default `false` — only β available.
"""
is_polarizable(::AbstractOceanPhaseFunction) = false

# -----------------------------------------------------------------------------
# Internal utilities
# -----------------------------------------------------------------------------

"""
    _reconstruct_from_moments(β, μ) -> FT

Evaluate `p(μ) = Σₗ βₗ Pₗ(μ)` (de Haan convention) using the three-term
Legendre recurrence. Internal; use `phase_function_value` instead.
"""
function _reconstruct_from_moments(β::AbstractVector{FT}, μ::FT) where {FT}
    n = length(β) - 1
    p = zero(FT)
    P_prev = zero(FT)
    P_curr = one(FT)                 # P₀(μ) = 1
    p += β[1] * P_curr
    if n >= 1
        P_next = μ                   # P₁(μ) = μ
        p += β[2] * P_next
        P_prev, P_curr = P_curr, P_next
        for ℓ in 2:n
            # (ℓ+1) P_{ℓ+1} = (2ℓ+1) μ P_ℓ − ℓ P_{ℓ-1}
            P_next = ((2ℓ - 1) * μ * P_curr - (ℓ - 1) * P_prev) / ℓ
            p += β[ℓ + 1] * P_next
            P_prev, P_curr = P_curr, P_next
        end
    end
    return p
end
