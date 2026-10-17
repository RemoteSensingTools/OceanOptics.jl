# =============================================================================
# Phase functions for ocean scattering — interface
# =============================================================================
#
# A phase function is a self-contained object that can evaluate p(μ) and
# produce its Greek expansion. Concrete forms:
#
#   RayleighWaterPhase       molecular scattering with depolarization (phase matrix)
#   FournierForandPhase      Fournier & Jonasz (1999) particle phase function
#   HenyeyGreensteinPhase    analytic one-parameter form
#   PetzoldPhase             Petzold (1972) San Diego Harbor, tabulated
#   TabulatedGreekPhase      given Greek coefficients (phase function or matrix)
#   TabulatedMuellerPhase    given angular scattering matrix (phase matrix)
#   WavelengthDependentPhase one of the above per tabulated wavelength
# =============================================================================

"""
$(TYPEDEF)

Abstract base type for scattering phase functions and phase matrices.

A phase function `p(μ)`, with `μ` the cosine of the scattering angle, is
normalized as `(1/2) ∫₋₁¹ p(μ) dμ = 1`. Expansions follow
`p(μ) = Σ_ℓ β_ℓ P_ℓ(μ)` with `β₀ = 1`; phase matrices use the layout and
signs of [`GreekCoefficients`](@ref).

# Interface

- `phase_function_value(pf, μ)` — `p(μ)`.
- [`phase_function_moments`](@ref)`(pf, ℓ_max)` — `β_ℓ`, `ℓ = 0:ℓ_max`.
  The default integrates `p` numerically in the scattering angle until
  converged; closed forms override it.
- [`greek_expansion`](@ref)`(pf, ℓ_max; polarized = false)` — the
  expansion as [`GreekCoefficients`](@ref). `polarized = true` requires
  [`is_polarizable`](@ref)`(pf)`: a scalar phase function does not
  determine a phase matrix.
- [`backscatter_fraction`](@ref)`(pf)`, [`asymmetry_parameter`](@ref)`(pf)`.
- [`is_wavelength_dependent`](@ref)`(pf)` — if true, every query takes a
  wavelength as an extra argument: `phase_function_moments(pf, ℓ_max, λ)`,
  `greek_expansion(pf, ℓ_max, λ; …)`, `backscatter_fraction(pf, λ)`.
  Wavelength-independent phase functions accept and ignore it.

New phase functions implement `phase_function_value` (or override
`phase_function_moments`); phase matrices also implement
`is_polarizable(pf) = true` and `_polarized_greek(pf, ℓ_max)`.
"""
abstract type AbstractOceanPhaseFunction{FT<:AbstractFloat} end

"""
    phase_function_value(pf, μ) -> FT

Phase function `p(μ)` at scattering cosine `μ`, normalized as
`(1/2) ∫₋₁¹ p dμ = 1`.
"""
function phase_function_value end

# Evaluation at scattering angle θ [rad]. Phase functions that are singular
# or peaked at θ → 0 override this to avoid forming cos θ near 1.
_phase_at_angle(pf::AbstractOceanPhaseFunction, θ::Real) = phase_function_value(pf, cos(θ))

# Kinks of tabulated phase functions (panel boundaries for the quadrature),
# whether p may be (integrably) singular at θ = 0, and the exponent a of its
# forward power law p ∝ θ^a (a > −2), used for the scattering below the
# quadrature floor; `nothing` if p is bounded there.
_angular_breaks(::AbstractOceanPhaseFunction) = ()
_forward_singular(::AbstractOceanPhaseFunction) = true
_forward_exponent(::AbstractOceanPhaseFunction) = nothing

# β_ℓ contribution of θ < θ_lo, where P_ℓ(cos θ) = 1: (2ℓ+1)/2 ∫₀^θ_lo p sin θ dθ
# with p = p(θ_lo) (θ/θ_lo)^a.
function _forward_tail(pf, θ_lo, ℓ_max)
    a = _forward_exponent(pf)
    (a === nothing || θ_lo == 0) && return zeros(ℓ_max + 1)
    m = _phase_at_angle(pf, θ_lo) * θ_lo^2 / (a + 2)
    return [(2ℓ + 1) / 2 * m for ℓ in 0:ℓ_max]
end

# β₀ tolerance: quadrature accuracy, or the rounding of a table stored in FT.
_normalization_tolerance(::Type{FT}) where {FT} = max(1e-8, 8 * eps(FT))

"""
    phase_function_moments(pf, ℓ_max[, λ]) -> Vector{FT}

Legendre expansion coefficients `β_ℓ`, `ℓ = 0:ℓ_max` (index 1 is `ℓ = 0`),
with `p(μ) = Σ_ℓ β_ℓ P_ℓ(μ)` and `β₀ = 1`.

The default evaluates `β_ℓ = (2ℓ+1)/2 ∫₀^π p(θ) P_ℓ(cos θ) sin θ dθ` by
composite Gauss–Legendre quadrature in `θ`, graded geometrically toward
the forward direction (with the analytic power-law tail below `1e-100`
rad where the phase function declares one), refining until `β_ℓ/(2ℓ+1)`
changes by less than `1e-10`. It throws if the integral does not converge
or if `β₀` differs from 1 by more than `max(1e-8, 8 eps(FT))` (the phase
function is not normalized); it never renormalizes silently.
"""
function phase_function_moments(pf::AbstractOceanPhaseFunction{FT}, ℓ_max::Integer) where {FT}
    ℓ_max ≥ 0 || throw(ArgumentError("ℓ_max must be non-negative, got $ℓ_max"))
    _requires_wavelength(pf)
    breaks   = _angular_breaks(pf)
    singular = _forward_singular(pf)
    β = _converged(ℓ_max; what = "Legendre expansion of $(nameof(typeof(pf)))") do level
        θ, w, θ_lo = _angular_rule(ℓ_max, level; breaks, forward_singular = singular)
        _legendre_moments([_phase_at_angle(pf, t) for t in θ], θ, w, ℓ_max) .+
            _forward_tail(pf, θ_lo, ℓ_max)
    end
    abs(β[1] - 1) ≤ _normalization_tolerance(FT) ||
        throw(ArgumentError("$(nameof(typeof(pf))) is not normalized: (1/2)∫p dμ = $(β[1])"))
    return FT.(β ./ β[1])
end

phase_function_moments(pf::AbstractOceanPhaseFunction, ℓ_max::Integer, λ::Real) =
    phase_function_moments(pf, ℓ_max)

"""
    greek_expansion(pf, ℓ_max[, λ]; polarized = false) -> GreekCoefficients

Expansion of `pf` up to `ℓ_max`. With `polarized = false` this is the
phase function only (`β`); with `polarized = true` it is the full phase
matrix, available only when [`is_polarizable`](@ref)`(pf)`. For
scalar-only phase functions the completely depolarizing approximation is
available explicitly through [`depolarizing_greek_expansion`](@ref).
"""
function greek_expansion(pf::AbstractOceanPhaseFunction, ℓ_max::Integer;
                            polarized::Bool = false)
    polarized || return GreekCoefficients(phase_function_moments(pf, ℓ_max))
    is_polarizable(pf) || throw(_scalar_only_error(pf))
    return _polarized_greek(pf, ℓ_max)
end

function greek_expansion(pf::AbstractOceanPhaseFunction, ℓ_max::Integer, λ::Real;
                            polarized::Bool = false)
    return greek_expansion(pf, ℓ_max; polarized)
end

"""
    depolarizing_greek_expansion(pf, ℓ_max[, λ]) -> GreekCoefficients

The phase function of `pf` with the *completely depolarizing* phase matrix
`diag(p(Θ), 0, 0, 0)` (see [`depolarized`](@ref)) — an explicit
approximation for scatterers whose phase matrix is unknown, not their
actual polarization behavior.
"""
depolarizing_greek_expansion(pf::AbstractOceanPhaseFunction, ℓ_max::Integer) =
    depolarized(greek_expansion(pf, ℓ_max))
depolarizing_greek_expansion(pf::AbstractOceanPhaseFunction, ℓ_max::Integer, λ::Real) =
    depolarized(greek_expansion(pf, ℓ_max, λ))

_polarized_greek(pf::AbstractOceanPhaseFunction, ℓ_max::Integer) =
    throw(ArgumentError("$(nameof(typeof(pf))) declares `is_polarizable` but does not " *
                        "implement `_polarized_greek`"))

_scalar_only_error(pf::AbstractOceanPhaseFunction) =
    ArgumentError("$(nameof(typeof(pf))) is scalar-only: its phase function does not " *
                  "determine the polarized Greek coefficients α, γ, δ, ϵ, ζ. Supply a " *
                  "phase matrix (`TabulatedMuellerPhase`, `TabulatedGreekPhase`) or " *
                  "select the completely depolarizing approximation explicitly.")

"""
    backscatter_fraction(pf[, λ]) -> FT

Fraction of scattering into the back hemisphere, `B = (1/2)∫₋₁⁰ p(μ) dμ`.
The default integrates `p` numerically; closed forms override it.
"""
function backscatter_fraction(pf::AbstractOceanPhaseFunction{FT}) where {FT}
    _requires_wavelength(pf)
    breaks = _angular_breaks(pf)
    B = _converged(16; what = "backscatter integral of $(nameof(typeof(pf)))") do level
        θ, w, _ = _angular_rule(16, level; breaks, forward_singular = false)
        [sum(w[i] * _phase_at_angle(pf, θ[i]) for i in eachindex(θ) if θ[i] > π / 2) / 2]
    end
    return FT(only(B))
end

backscatter_fraction(pf::AbstractOceanPhaseFunction, λ::Real) = backscatter_fraction(pf)

"""
    asymmetry_parameter(pf[, λ]) -> FT

Mean cosine of the scattering angle, `g = ⟨μ⟩ = β₁ / 3`.
"""
asymmetry_parameter(pf::AbstractOceanPhaseFunction{FT}) where {FT} =
    phase_function_moments(pf, 1)[2] / FT(3)

asymmetry_parameter(pf::AbstractOceanPhaseFunction, λ::Real) = asymmetry_parameter(pf)

"""
    is_polarizable(pf) -> Bool

Whether `pf` provides a full phase matrix, so that
`greek_expansion(pf, ℓ_max; polarized = true)` is defined.
"""
is_polarizable(::AbstractOceanPhaseFunction) = false

"""
    is_wavelength_dependent(pf) -> Bool

Whether `pf` differs between wavelengths; such phase functions must be
queried with an explicit wavelength. Default `false`.
"""
is_wavelength_dependent(::AbstractOceanPhaseFunction) = false

# Key under which the expansion of a wavelength-dependent phase function at λ
# can be cached: equal keys ⇒ identical expansions.
_phase_key(::AbstractOceanPhaseFunction, λ::Real) = λ

_requires_wavelength(pf) = is_wavelength_dependent(pf) &&
    throw(ArgumentError("$(nameof(typeof(pf))) depends on wavelength; pass the " *
                        "wavelength [nm] as an extra argument"))

"""
    _reconstruct_from_moments(β, μ) -> FT

`p(μ) = Σₗ βₗ Pₗ(μ)` by the three-term Legendre recurrence.
"""
function _reconstruct_from_moments(β::AbstractVector{FT}, μ::Real) where {FT}
    x = FT(μ)
    p = β[1]
    length(β) == 1 && return p
    P_prev, P = one(FT), x
    p += β[2] * P
    for ℓ in 2:length(β)-1
        P_prev, P = P, ((2ℓ - 1) * x * P - (ℓ - 1) * P_prev) / ℓ
        p += β[ℓ+1] * P
    end
    return p
end
