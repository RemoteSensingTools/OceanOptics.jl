# =============================================================================
# Fournier-Forand phase function
# =============================================================================
#
# Fournier & Jonasz (1999) closed-form phase function for an ensemble of
# particles with a Junge (hyperbolic) size distribution, each scattering in
# the anomalous-diffraction approximation. Two parameters:
#
#   n     — real refractive index of the particles relative to water
#   μ_j   — Junge slope of the size distribution, 3 < μ_j < 5
#
# The closed form is normalized exactly, ∫ p dΩ / (4π) = 1. Its forward
# peak is an integrable singularity, p ~ θ^(2|ν|−2) with ν = (3 − μ_j)/2, so
# it is evaluated through s = sin²(θ/2) rather than through cos θ, and its
# expansion is integrated in θ (see `angular.jl`).
#
# Given only a backscatter fraction B, (n, μ_j) follow the one-parameter
# family of Mobley, Sundman & Boss (2002), n = 1.01 + 0.1542 (μ_j − 3),
# which also underlies HydroLight's Fournier-Forand option; n = 1.10,
# μ_j = 3.5835 gives B = 0.0183 (Petzold average particle).
#
# References
# - Fournier & Forand, Proc. SPIE 2258, 194 (1994)
# - Fournier & Jonasz, Proc. SPIE 3761, 62 (1999)
# - Mobley, Sundman & Boss, Appl. Opt. 41, 1035 (2002)
# =============================================================================

"""
$(TYPEDEF)

Fournier-Forand phase function for ocean particles.

# Fields
$(TYPEDFIELDS)

# Construction

```julia
FournierForandPhase{Float64}(n = 1.10, μ_j = 3.5835)       # B ≈ 0.0183
FournierForandPhase{Float64}(backscatter_fraction = 0.01)  # Mobley et al. (2002) family
```

With `backscatter_fraction = B`, `(n, μ_j)` are found on the family
`n = 1.01 + 0.1542 (μ_j − 3)` of Mobley, Sundman & Boss (2002); `B` must
lie in the range that family spans for `3 < μ_j < 5`. The field is called
`μ_j` (Junge slope) because `μ` denotes `cos θ` throughout this package.
"""
struct FournierForandPhase{FT} <: AbstractOceanPhaseFunction{FT}
    "Real refractive index of the particles relative to water"
    n::FT
    "Junge slope of the particle size distribution, 3 < μ_j < 5"
    μ_j::FT

    function FournierForandPhase{FT}(n::FT, μ_j::FT) where {FT<:AbstractFloat}
        # n < 1.5 keeps clear of the poles of the closed form at δ₉₀ = 1
        # (n ≈ 1.82) and δ₁₈₀ = 1 (n ≈ 2.15); ocean particles have n ≲ 1.3.
        1 < n < 1.5 || throw(ArgumentError("FournierForand: n must lie in (1, 1.5), got $n"))
        FT(3) < μ_j < FT(5) ||
            throw(ArgumentError("FournierForand: Junge slope must lie in (3, 5), got $μ_j"))
        new{FT}(n, μ_j)
    end
end

function FournierForandPhase{FT}(; n = nothing, μ_j = nothing,
                                   backscatter_fraction = nothing) where {FT<:AbstractFloat}
    if backscatter_fraction !== nothing
        (n === nothing && μ_j === nothing) ||
            throw(ArgumentError("FournierForand: give either (n, μ_j) or backscatter_fraction"))
        nn, mj = _ff_mobley_parameters(Float64(backscatter_fraction))
        return FournierForandPhase{FT}(FT(nn), FT(mj))
    end
    (n === nothing || μ_j === nothing) &&
        throw(ArgumentError("FournierForand: give (n, μ_j) or backscatter_fraction"))
    return FournierForandPhase{FT}(FT(n), FT(μ_j))
end

FournierForandPhase(; kwargs...) = FournierForandPhase{Float64}(; kwargs...)

# p as a function of s = sin²(θ/2), normalized as (1/2)∫p dμ = 1.
function _ff_value(n, μ_j, s)
    s     = max(s, floatmin(Float64))                    # finite at θ = 0
    ν     = (3 - μ_j) / 2
    c     = 3 * (n - 1)^2 / 4                            # δ = s / c
    δ     = s / c
    δ_ν   = δ^ν
    x     = 1 - δ
    main  = if abs(x) < 1e-2
        _ff_main_series(ν, c, x) / (4π * δ_ν)            # removable 0/0 at δ = 1
    else
        (ν * x - (1 - δ_ν) + (δ * (1 - δ_ν) - ν * x) / s) / (4π * x^2 * δ_ν)
    end
    δ_180 = 1 / c
    back  = (1 - δ_180^ν) / (16π * (δ_180 - 1) * δ_180^ν) * (3 * (1 - 2s)^2 - 1)
    return 4π * (main + back)
end

# Numerator of the main term divided by (1 − δ)², as a power series in
# x = 1 − δ: with 1 − (1−x)^ν = Σ bₖ xᵏ,
#   N / x² = Σ_{k≥2} [−bₖ + (bₖ − ν)/c] x^(k−2),
# which avoids the cancellation of the closed form at δ → 1 (θ ≈ 10° for
# n = 1.1). Twelve terms are exact to round-off for |x| < 1e-2.
function _ff_main_series(ν, c, x)
    C   = one(x)                      # generalized binomial C(ν, k)
    sum = zero(x)
    xp  = one(x)
    for k in 1:12
        C *= (ν - k + 1) / k
        b  = -C * (-1)^k
        k ≥ 2 || continue
        sum += (-b + (b - ν) / c) * xp
        xp  *= x
    end
    return sum
end

"""
    phase_function_value(pf::FournierForandPhase, μ)

Fournier & Jonasz (1999) phase function, normalized as
`(1/2)∫₋₁¹ p dμ = 1`. Diverges (integrably) at `μ → 1`; values there are
large but finite.
"""
phase_function_value(pf::FournierForandPhase, μ::Real) = _ff_value(pf.n, pf.μ_j, (1 - μ) / 2)

_phase_at_angle(pf::FournierForandPhase, θ::Real) = _ff_value(pf.n, pf.μ_j, sin(θ / 2)^2)

# p ∝ θ^(μ_j − 5) as θ → 0 (forward tail below the quadrature floor).
_forward_exponent(pf::FournierForandPhase) = pf.μ_j - 5

"""
    backscatter_fraction(pf::FournierForandPhase) -> FT

Closed form (Fournier & Jonasz 1999; Mobley et al. 2002 Eq. 10):
`B = 1 − [1 − δ₉₀^(ν+1) − ½(1 − δ₉₀^ν)] / [(1 − δ₉₀) δ₉₀^ν]`.
"""
backscatter_fraction(pf::FournierForandPhase{FT}) where {FT} = FT(_ff_backscatter(pf.n, pf.μ_j))

function _ff_backscatter(n, μ_j)
    ν    = (3 - μ_j) / 2
    δ_90 = (4 / 3) * (1 / 2) / (n - 1)^2
    return 1 - (1 - δ_90^(ν + 1) - (1 - δ_90^ν) / 2) / ((1 - δ_90) * δ_90^ν)
end

# Mobley, Sundman & Boss (2002) parameter family.
_ff_mobley_n(μ_j) = 1.01 + 0.1542 * (μ_j - 3)

"""
    _ff_mobley_parameters(B) -> (n, μ_j)

Solve `B(n(μ_j), μ_j) = B` on the Mobley et al. (2002) family by bisection
in `μ_j ∈ (3, 5)`; `B` is monotone along it.
"""
function _ff_mobley_parameters(B::Float64)
    f(μ) = _ff_backscatter(_ff_mobley_n(μ), μ) - B
    lo, hi = 3.0 + 1e-9, 5.0 - 1e-9
    flo, fhi = f(lo), f(hi)
    flo < 0 < fhi ||
        throw(ArgumentError("FournierForand: backscatter fraction $B is outside the range " *
                            "$(B + flo)–$(B + fhi) spanned by the Mobley et al. (2002) " *
                            "family; give (n, μ_j) explicitly"))
    for _ in 1:200
        mid = (lo + hi) / 2
        f(mid) < 0 ? (lo = mid) : (hi = mid)
        hi - lo < 1e-14 && break
    end
    μ = (lo + hi) / 2
    return _ff_mobley_n(μ), μ
end
