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

# Open interval of the Junge slope μ_j over which the phase function is defined.
const _FF_JUNGE_SLOPE_RANGE = (3.0, 5.0)

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
lie in the range that family spans for `3 < μ_j < 5`. A dual-number `B` gives dual
`(n, μ_j)` with the derivative of that inversion, which the phase function, its
backscatter fraction and its Legendre expansion carry. The field is called
`μ_j` (Junge slope) because `μ` denotes `cos θ` throughout this package.
"""
struct FournierForandPhase{FT} <: AbstractOceanPhaseFunction{FT}
    "Real refractive index of the particles relative to water"
    n::FT
    "Junge slope of the particle size distribution, 3 < μ_j < 5"
    μ_j::FT

    function FournierForandPhase{FT}(n::FT, μ_j::FT) where {FT<:Real}
        _check_parameter_type(FT)
        # n < 1.5 keeps clear of the poles of the closed form at δ₉₀ = 1
        # (n ≈ 1.82) and δ₁₈₀ = 1 (n ≈ 2.15); ocean particles have n ≲ 1.3.
        1 < n < 1.5 || throw(ArgumentError("FournierForand: n must lie in (1, 1.5), got $n"))
        FT(first(_FF_JUNGE_SLOPE_RANGE)) < μ_j < FT(last(_FF_JUNGE_SLOPE_RANGE)) ||
            throw(ArgumentError("FournierForand: Junge slope must lie in $(_FF_JUNGE_SLOPE_RANGE), got $μ_j"))
        new{FT}(n, μ_j)
    end
end

function FournierForandPhase{FT}(; n = nothing, μ_j = nothing,
                                   backscatter_fraction = nothing) where {FT<:Real}
    _check_parameter_type(FT)
    if backscatter_fraction !== nothing
        (n === nothing && μ_j === nothing) ||
            throw(ArgumentError("FournierForand: give either (n, μ_j) or backscatter_fraction"))
        nn, mj = _ff_mobley_parameters(backscatter_fraction)
        return FournierForandPhase{FT}(FT(nn), FT(mj))
    end
    (n === nothing || μ_j === nothing) &&
        throw(ArgumentError("FournierForand: give (n, μ_j) or backscatter_fraction"))
    return FournierForandPhase{FT}(FT(n), FT(μ_j))
end

FournierForandPhase(; kwargs...) = FournierForandPhase{_parameter_type(values(kwargs))}(; kwargs...)

# p as a function of s = sin²(θ/2), normalized as (1/2)∫p dμ = 1.
function _ff_value(n, μ_j, s)
    s     = max(s, floatmin(Float64))                    # finite at θ = 0
    ν     = (3 - μ_j) / 2
    c     = 3 * (n - 1)^2 / 4                            # δ = s / c
    δ     = s / c
    δ_ν   = _ff_power(δ, s, c, ν)
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

# δ^ν for δ = s/c. Floating-point arguments use `^`. Any other number (a dual) takes
# exp(ν (log s − log c)): ForwardDiff's rule for x^y forms x^(y−1) and, in a second
# derivative, x^(y−2), which overflows at the forward floor s ≈ 1e-200 and gives
# ∞·0 = NaN, although every derivative of δ^ν is of the size of δ^ν log^k δ.
_ff_power(δ::AbstractFloat, s, c::AbstractFloat, ν::AbstractFloat) = δ^ν
_ff_power(δ, s, c, ν) = exp(ν * (log(s) - log(c)))

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

Return the phase function of Fournier & Jonasz (1999) at scattering cosine `μ`,
normalized as `(1/2)∫₋₁¹ p dμ = 1`.

It diverges (integrably) for `μ → 1`; values there are large but finite.
"""
phase_function_value(pf::FournierForandPhase{FT}, μ::Real) where {FT} =
    convert(_result_type(FT, μ), _ff_value(_wide(pf.n), _wide(pf.μ_j), (1 - _wide(μ)) / 2))

_phase_at_angle(pf::FournierForandPhase, θ::Real) = _ff_value(pf.n, pf.μ_j, sin(θ / 2)^2)

# p ∝ θ^(μ_j − 5) as θ → 0 (forward tail below the quadrature floor).
_forward_exponent(pf::FournierForandPhase) = pf.μ_j - 5

"""
    backscatter_fraction(pf::FournierForandPhase) -> FT

Return the backscatter fraction in closed form (Fournier & Jonasz 1999; Mobley et al.
2002 Eq. 10): `B = 1 − [1 − δ₉₀^(ν+1) − ½(1 − δ₉₀^ν)] / [(1 − δ₉₀) δ₉₀^ν]`.
"""
backscatter_fraction(pf::FournierForandPhase{FT}) where {FT} = FT(_ff_backscatter(pf.n, pf.μ_j))

function _ff_backscatter(n, μ_j)
    ν    = (3 - μ_j) / 2
    δ_90 = (4 / 3) * (1 / 2) / (n - 1)^2
    return 1 - (1 - δ_90^(ν + 1) - (1 - δ_90^ν) / 2) / ((1 - δ_90) * δ_90^ν)
end

# Mobley, Sundman & Boss (2002) parameter family n = n₀ + s (μ_j − μ₀), μ_j ∈ (μ₀, 5).
const _FF_MOBLEY_2002_N0     = 1.01
const _FF_MOBLEY_2002_SLOPE  = 0.1542
const _FF_MOBLEY_2002_μ0     = first(_FF_JUNGE_SLOPE_RANGE)   # the family starts at μ_j = 3
_ff_mobley_n(μ_j) = _FF_MOBLEY_2002_N0 + _FF_MOBLEY_2002_SLOPE * (μ_j - _FF_MOBLEY_2002_μ0)

"""
    _ff_mobley_parameters(B) -> (n, μ_j)

Solve `B(n(μ_j), μ_j) = B` on the Mobley et al. (2002) family; `B` is monotone along it.

The root `μ₀` is found by bisection in `μ_j ∈ (3, 5)`. For a `B` that is not a
floating-point number (a dual number) the root is refined in the number type of `B` by
the chord iteration `μ ← μ − (B(μ) − B) / B′(μ₀)`, started at `μ₀`, with `B′(μ₀)`
from the complex-step derivative. The value stays `μ₀` (to round-off); each step makes
one more order of derivative exact, so the `_FF_CHORD_STEPS` steps give exact
derivatives of `μ_j` with respect to `B` up to that order (nested duals, Hessians).
"""
function _ff_mobley_parameters(B::AbstractFloat)
    μ = _ff_mobley_root(Float64(B))
    return _ff_mobley_n(μ), μ
end

function _ff_mobley_parameters(B::Real)
    μ₀ = _ff_mobley_root(B)
    B₀ = _ff_backscatter(_ff_mobley_n(μ₀), μ₀)
    # B′(μ₀) along the family by the complex-step derivative, exact to round-off.
    z  = complex(μ₀, _COMPLEX_STEP)
    dB = imag(_ff_backscatter(_ff_mobley_n(z), z)) / _COMPLEX_STEP
    μ  = μ₀ + (B - B₀) / dB
    for _ in 2:_FF_CHORD_STEPS
        μ -= (_ff_backscatter(_ff_mobley_n(μ), μ) - B) / dB
    end
    return _ff_mobley_n(μ), μ
end

"""
Number of chord steps of `_ff_mobley_parameters(B::Real)`; derivatives of `μ_j(B)` up
to this order are exact, higher ones are not.
"""
const _FF_CHORD_STEPS = 4

"Imaginary increment of the complex-step derivative `f′(x) ≈ Im f(x + ih) / h`."
const _COMPLEX_STEP = 1e-20

# Root μ₀::Float64 of B(n(μ), μ) = B by bisection; only the value of `B` is compared.
function _ff_mobley_root(B::Real)
    f(μ) = _ff_backscatter(_ff_mobley_n(μ), μ) - B
    lo, hi = first(_FF_JUNGE_SLOPE_RANGE) + 1e-9, last(_FF_JUNGE_SLOPE_RANGE) - 1e-9
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
    return (lo + hi) / 2
end
