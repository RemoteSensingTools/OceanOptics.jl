# =============================================================================
# Tabulated angular scattering matrix
# =============================================================================
#
# Phase matrix given as F(Θ) on an angular grid covering 0°–180°, as in
# published benchmark supplements and Mie or T-matrix output, for
# macroscopically isotropic, mirror-symmetric media:
#
#   F = [F11 F12  0    0
#        F12 F22  0    0
#        0   0    F33  F34
#        0   0   -F34  F44]
#
# in the standard convention of Hansen & Travis (1974) and Mishchenko et al.
# (2002): F₁₂ < 0 for Rayleigh scattering at 90°. (Chowdhary et al. 2020
# tabulate F21 = F12 and F43 = −F34.)
#
# Between nodes log F11 and the ratios Fij/F11 are interpolated with
# shape-preserving monotone cubics (Fritsch & Carlson 1980, "PCHIP"): each
# interpolant stays within the range of its two neighbouring nodes, so
# F11 > 0 and |Fij| ≤ F11 hold everywhere. Pointwise errors on smooth
# tables are O(h²) near extrema and smaller elsewhere; the expansion
# coefficients converge faster (≈1e-8 for a 1° Rayleigh table). The expansion integrals are evaluated on
# Gauss-Legendre panels aligned with the table nodes (`angular.jl`) and
# mapped to `GreekCoefficients`.
# =============================================================================

"""
$(TYPEDEF)

Phase matrix given by its angular scattering matrix `F(Θ)`.

# Fields
$(TYPEDFIELDS)

# Construction

    TabulatedMuellerPhase(θ_deg; F11, F12, F33, F34, F22 = F11, F44 = F33,
                          normalization_tolerance = 1e-3, normalize = false,
                          source = "")

`θ_deg` is strictly increasing from 0 to 180. Elements are passed by name
so they cannot be swapped by position; `F22 = F11` and `F44 = F33` hold for
spheres. All elements use the standard convention above (`F12 < 0` for
Rayleigh scattering) and the normalization `(1/2)∫F11 dμ = 1`, i.e.
`(1/4π)∫F11 dΩ = 1`. Tables of `F/(4π)` normalized as `∫ F/(4π) dΩ = 1`
(e.g. Chowdhary et al. 2020, Table S1) must be multiplied by `4π` first.

The integral `(1/2)∫F11 dμ` of the interpolated table is stored as
`normalization`. The table is divided by it so that `β₀ = 1`; this is
allowed only when `|normalization − 1| ≤ normalization_tolerance`, which
covers rounding of published values, or with `normalize = true`.

Every element must satisfy `|Fij| ≤ F11` (a necessary condition for a
physical scattering matrix), within `1e-6` relative.

# Example

```julia
# Rayleigh scattering (no depolarization), standard convention
θ  = 0.0:1.0:180.0
μ  = cosd.(θ)
pf = TabulatedMuellerPhase(θ; F11 = 0.75 .* (1 .+ μ .^ 2), F12 = -0.75 .* (1 .- μ .^ 2),
                           F33 = 1.5 .* μ, F34 = zero(μ))
greek_expansion(pf, 4; polarized = true).γ[3]   # ≈ √(3/2), as get_greek_rayleigh(0)
```
"""
struct TabulatedMuellerPhase{FT} <: AbstractOceanPhaseFunction{FT}
    "Scattering angles [rad], strictly increasing from 0 to π"
    θ::Vector{FT}
    "F11 (phase function), normalized so that (1/2)∫F11 dμ = 1"
    F11::Vector{FT}
    "F12 = F21"
    F12::Vector{FT}
    "F22"
    F22::Vector{FT}
    "F33"
    F33::Vector{FT}
    "F34 = −F43"
    F34::Vector{FT}
    "F44"
    F44::Vector{FT}
    "(1/2)∫F11 dμ of the table as given, before normalization"
    normalization::FT
    "Monotone-cubic slopes of log F11 and of F12/F11, …, F44/F11 at the nodes"
    slopes::NTuple{6, Vector{FT}}
    "Provenance: citation, table, and any transformation applied"
    source::String
end

function TabulatedMuellerPhase(θ_deg::AbstractVector; F11, F12, F33, F34, F22 = F11, F44 = F33,
                               normalization_tolerance::Real = 1e-3, normalize::Bool = false,
                               source::AbstractString = "")
    who = "TabulatedMuellerPhase"
    FT  = float(promote_type(eltype(θ_deg), eltype(F11)))
    θ   = deg2rad.(collect(FT, θ_deg))
    n   = length(θ)
    n ≥ 2 || throw(ArgumentError("$who: need at least two angles"))
    _is_strictly_increasing(θ) || throw(ArgumentError("$who: angles must be strictly increasing"))
    abs(first(θ_deg)) ≤ 1e-9 && abs(last(θ_deg) - 180) ≤ 1e-9 ||
        throw(ArgumentError("$who: the table must cover 0° to 180°, got $(first(θ_deg))°–$(last(θ_deg))°"))
    els = (F11 = F11, F12 = F12, F22 = F22, F33 = F33, F34 = F34, F44 = F44)
    for (k, v) in pairs(els)
        length(v) == n || throw(DimensionMismatch("$who: $k has $(length(v)) values for $n angles"))
        all(isfinite, v) || throw(ArgumentError("$who: $k contains non-finite values"))
    end
    all(>(0), F11) || throw(ArgumentError("$who: F11 must be positive"))
    for k in (:F12, :F22, :F33, :F34, :F44), i in 1:n
        abs(els[k][i]) ≤ F11[i] * (1 + 1e-6) ||
            throw(ArgumentError("$who: |$k| > F11 at $(θ_deg[i])°; not a physical scattering " *
                                "matrix (wrong element, units or normalization?)"))
    end
    raw = TabulatedMuellerPhase{FT}(θ, (collect(FT, els[k]) for k in keys(els))..., one(FT),
                                    _mueller_slopes(θ, els), "")
    N = only(_converged(16; what = "$who normalization") do level
        t, w, _ = _angular_rule(16, level; breaks = θ, forward_singular = false)
        [sum(w[i] * _phase_at_angle(raw, t[i]) for i in eachindex(t)) / 2]
    end)
    normalize || abs(N - 1) ≤ normalization_tolerance ||
        throw(ArgumentError("$who: (1/2)∫F11 dμ = $N differs from 1 by more than " *
                            "$normalization_tolerance; check the normalization (tables of " *
                            "F/(4π) need a factor 4π) or pass normalize = true"))
    scale = v -> collect(FT, v) ./ N
    return TabulatedMuellerPhase{FT}(θ, (scale(els[k]) for k in keys(els))..., FT(N),
                                     raw.slopes, String(source))
end

# Interpolated quantities: log F11 and the ratios Fij/F11 (scale invariant).
function _mueller_slopes(θ, els)
    logF11 = log.(els.F11)
    ratios = (els[k] ./ els.F11 for k in (:F12, :F22, :F33, :F34, :F44))
    return (_pchip_slopes(θ, logF11), (_pchip_slopes(θ, r) for r in ratios)...)
end

"""
    _pchip_slopes(x, y) -> Vector

Node derivatives of the Fritsch–Carlson monotone piecewise-cubic Hermite
interpolant (as in MATLAB's `pchip`): zero at local extrema, a weighted
harmonic mean of the adjacent secants elsewhere, and shape-preserving
three-point formulas at the ends.
"""
function _pchip_slopes(x::AbstractVector, y::AbstractVector)
    n = length(x)
    h = diff(x)
    δ = diff(y) ./ h
    d = zeros(eltype(δ), n)
    n == 2 && (d .= δ[1]; return d)
    for k in 2:n-1
        if δ[k-1] * δ[k] > 0
            w1, w2 = 2h[k] + h[k-1], h[k] + 2h[k-1]
            d[k] = (w1 + w2) / (w1 / δ[k-1] + w2 / δ[k])
        end
    end
    d[1] = _pchip_end(h[1], h[2], δ[1], δ[2])
    d[n] = _pchip_end(h[n-1], h[n-2], δ[n-1], δ[n-2])
    return d
end

function _pchip_end(h1, h2, δ1, δ2)
    d = ((2h1 + h2) * δ1 - h1 * δ2) / (h1 + h2)
    sign(d) != sign(δ1) && return zero(d)
    sign(δ1) != sign(δ2) && abs(d) > abs(3δ1) && return 3δ1
    return d
end

@inline function _hermite(y0, y1, d0, d1, h, t)
    t2, t3 = t^2, t^3
    return (2t3 - 3t2 + 1) * y0 + (t3 - 2t2 + t) * h * d0 + (-2t3 + 3t2) * y1 + (t3 - t2) * h * d1
end

# Interpolation: log-linear F11, linear ratios Fij/F11.
function _mueller_interval(pf::TabulatedMuellerPhase, θ::Real)
    i = clamp(searchsortedlast(pf.θ, θ), 1, length(pf.θ) - 1)
    t = clamp((θ - pf.θ[i]) / (pf.θ[i+1] - pf.θ[i]), 0, 1)
    return i, t
end

function _phase_at_angle(pf::TabulatedMuellerPhase, θ::Real)
    i, t = _mueller_interval(pf, θ)
    h = pf.θ[i+1] - pf.θ[i]
    return exp(_hermite(log(pf.F11[i]), log(pf.F11[i+1]), pf.slopes[1][i], pf.slopes[1][i+1], h, t))
end

function _mueller_at(pf::TabulatedMuellerPhase, θs::AbstractVector)
    m = length(θs)
    F = (F11 = zeros(m), F12 = zeros(m), F22 = zeros(m), F33 = zeros(m), F34 = zeros(m), F44 = zeros(m))
    for (j, θ) in enumerate(θs)
        i, t = _mueller_interval(pf, θ)
        h   = pf.θ[i+1] - pf.θ[i]
        f11 = exp(_hermite(log(pf.F11[i]), log(pf.F11[i+1]), pf.slopes[1][i], pf.slopes[1][i+1], h, t))
        F.F11[j] = f11
        for (s, k) in enumerate((:F12, :F22, :F33, :F34, :F44))
            v = getfield(pf, k)
            r = _hermite(v[i] / pf.F11[i], v[i+1] / pf.F11[i+1],
                         pf.slopes[s+1][i], pf.slopes[s+1][i+1], h, t)
            getfield(F, k)[j] = r * f11
        end
    end
    return F
end

phase_function_value(pf::TabulatedMuellerPhase{FT}, μ::Real) where {FT} =
    FT(_phase_at_angle(pf, acos(clamp(μ, -1, 1))))

_angular_breaks(pf::TabulatedMuellerPhase) = pf.θ
_forward_singular(::TabulatedMuellerPhase) = false
is_polarizable(::TabulatedMuellerPhase) = true

function _polarized_greek(pf::TabulatedMuellerPhase{FT}, ℓ_max::Integer) where {FT}
    ℓ_max ≥ 0 || throw(ArgumentError("ℓ_max must be non-negative, got $ℓ_max"))
    m = _converged(ℓ_max; what = "Greek expansion of TabulatedMuellerPhase") do level
        θ, w, _ = _angular_rule(ℓ_max, level; breaks = pf.θ, forward_singular = false)
        _mueller_moments(_mueller_at(pf, θ), θ, w, ℓ_max)
    end
    β₀ = m[1][1]
    abs(β₀ - 1) ≤ _normalization_tolerance(FT) ||
        throw(ErrorException("TabulatedMuellerPhase: expansion gives β₀ = $β₀ for a table " *
                             "normalized to 1"))
    return _greek_from_mishchenko(FT, (x ./ β₀ for x in m)...)   # remove round-off in β₀
end

function Base.show(io::IO, pf::TabulatedMuellerPhase{FT}) where {FT}
    print(io, "TabulatedMuellerPhase{$FT}(", length(pf.θ), " angles, normalization ",
          round(pf.normalization; sigdigits = 6),
          isempty(pf.source) ? "" : ", \"$(pf.source)\"", ")")
end
