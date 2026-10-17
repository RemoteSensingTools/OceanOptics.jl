# =============================================================================
# Petzold (1972) San Diego Harbor phase function
# =============================================================================
#
# Tabulated volume scattering function at 55 angles, 0.1°–180°. Between
# nodes it is interpolated as a power law (linear in log p vs log θ), and
# below 0.1° the first segment's power law is continued into the forward
# peak (p ∝ θ^a with a ≈ −1.3, integrable). The table is normalized with
# the same converged angular quadrature used for its moments, so
# (1/2)∫p dμ = 1 holds for the interpolant as implemented; the forward
# extrapolation is a modelling choice that affects that normalization.
# =============================================================================

"""
$(TYPEDEF)

Petzold (1972) "San Diego Harbor" particulate phase function (Mobley 1994,
Table 3.10), a scalar phase function.

# Fields
$(TYPEDFIELDS)

Power-law interpolation between the 55 tabulated angles; power-law
extrapolation below 0.1° using the first two nodes; normalized as
`(1/2)∫p dμ = 1` for this interpolant.
"""
struct PetzoldPhase{FT} <: AbstractOceanPhaseFunction{FT}
    "Tabulated scattering angles [rad], ascending, from 0.1° to 180°"
    θ::Vector{FT}
    "Normalized phase function at those angles"
    p::Vector{FT}
    "Scattering coefficient implied by the tabulated VSF and its extrapolation, `4π·(1/2)∫β dμ` `[1/m]`"
    b::FT
end

const _PETZOLD_SDH_TABLE = load_spectral_csv(data_path("petzold_sdh.csv"))   # θ_deg, β

function PetzoldPhase{FT}() where {FT<:AbstractFloat}
    θ  = deg2rad.(_PETZOLD_SDH_TABLE.λ)
    β  = _PETZOLD_SDH_TABLE.y
    raw = PetzoldPhase{Float64}(θ, β, NaN)
    a  = _petzold_forward_exponent(raw)
    a > -2 || throw(ErrorException("Petzold forward extrapolation is not integrable (a = $a)"))
    I = only(_converged(16; what = "Petzold normalization") do level
        t, w, θ_lo = _angular_rule(16, level; breaks = θ)
        [sum(w[i] * _phase_at_angle(raw, t[i]) for i in eachindex(t)) / 2 +
         _forward_tail(raw, θ_lo, 0)[1]]
    end)
    return PetzoldPhase{FT}(FT.(θ), FT.(β ./ I), FT(4π * I))
end

PetzoldPhase() = PetzoldPhase{Float64}()

_petzold_forward_exponent(pf::PetzoldPhase) =
    log(pf.p[2] / pf.p[1]) / log(pf.θ[2] / pf.θ[1])

function _phase_at_angle(pf::PetzoldPhase, θ::Real)
    θs, ps = pf.θ, pf.p
    if θ ≤ first(θs)        # power-law forward peak; finite at θ = 0
        return first(ps) * (max(θ, 1e-300) / first(θs))^_petzold_forward_exponent(pf)
    end
    θ ≥ last(θs) && return last(ps)
    i = searchsortedlast(θs, θ)
    t = log(θ / θs[i]) / log(θs[i+1] / θs[i])
    return ps[i] * (ps[i+1] / ps[i])^t
end

phase_function_value(pf::PetzoldPhase{FT}, μ::Real) where {FT} =
    FT(_phase_at_angle(pf, acos(clamp(μ, -1, 1))))

_angular_breaks(pf::PetzoldPhase) = pf.θ
_forward_exponent(pf::PetzoldPhase) = _petzold_forward_exponent(pf)
