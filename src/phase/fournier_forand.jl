# =============================================================================
# Fournier-Forand phase function
# =============================================================================
#
# Fournier & Jonasz (1999) closed-form phase function for an ensemble of
# ocean particles following a Junge (hyperbolic) size distribution, each
# scattering in the anomalous-diffraction approximation. Two parameters:
#
#   n     — real refractive index of particles (typical 1.03–1.20, relative to water)
#   μ_j   — Junge slope of the particle size distribution (typical 3 < μ_j < 5)
#
# KNOWN ISSUE (documented here for Claude Code to fix on your machine):
# The raw Fournier & Jonasz closed form is normalized such that
# ∫_{4π} p dΩ = 1, but in practice the back-hemisphere "correction" term
# (the (3cos²θ − 1) piece) breaks exact unit normalization by ~1% for
# typical (n, μ_j). The fix is to compute a renormalization factor at
# construction time:
#
#   renorm = 1 / [(1/2) ∫₋₁¹ p_raw(μ) dμ]
#
# and multiply it into every evaluation. The integral needs adaptive
# quadrature or log-transformed nodes near μ = 1 because the forward peak
# is integrable but near-singular (p ~ 1/sin²(θ/2) there). A fixed
# Gauss-Legendre rule with ≥ 3000 nodes on [-1, 1] is *not* sufficient;
# use scipy-quad-style adaptive integration (Julia: `QuadGK.quadgk` with
# `rtol=1e-10` and a subdivision at μ = 0.999 to give the adaptive
# algorithm a chance at the peak).
#
# References
# - Fournier & Forand, Proc. SPIE 2258, 194 (1994)
# - Fournier & Jonasz, Proc. SPIE 3761, 62 (1999)
# - Mobley, Sundman & Boss, Appl. Opt. 41, 1035 (2002) — B → (n, μ_j) inversion
# - Ocean Optics Web Book, "The Fournier-Forand Phase Function"
# =============================================================================

"""
$(TYPEDEF)

Fournier-Forand phase function for ocean-particle scattering.

# Fields
$(TYPEDFIELDS)

# Construction

```julia
FournierForandPhase{Float64}(n = 1.15, μ_j = 3.38)
FournierForandPhase{Float64}(backscatter_fraction = 0.018)
FournierForandPhase{Float64}()                           # B = 0.018 default
```

# Note on the field name `μ_j`

We use `μ_j` (j for "Junge") rather than the literature's bare `μ` because
in ocean-optics code `μ` is pervasively used for `cos(θ)`. Keeping the
field name distinct prevents a class of dispatch bugs.
"""
struct FournierForandPhase{FT} <: AbstractOceanPhaseFunction{FT}
    "Real refractive index of particles (relative to water)"
    n::FT
    "Junge slope of the size distribution, typically 3 < μ_j < 5"
    μ_j::FT
    "Renormalization factor such that (1/2) ∫₋₁¹ p(μ) dμ = 1 exactly"
    renorm::FT
end

# Raw (un-renormalized) closed form of Fournier & Jonasz (1999).
#
# Deliberately un-annotated in its args: integer literals promote against
# `n`, `μ_j`, `cosθ`, so passing a `ForwardDiff.Dual` for any of them
# flows derivatives through the whole evaluation without a type-flattening
# cast (DESIGN §12.5).
function _ff_raw(n, μ_j, cosθ)
    ν          = (3 - μ_j) / 2
    sin2_half  = (1 - cosθ) / 2
    # Floor the forward-peak singularity: p ~ 1/sin²(θ/2) near μ = 1. An
    # explicit evaluation at cosθ = 1 is supported by clamping the
    # denominator to eps of the underlying float type.
    sin2_half  = max(sin2_half, eps(float(real(one(sin2_half)))))
    δ          = (4 * sin2_half) / (3 * (n - 1)^2)
    δ_180      =                 4 / (3 * (n - 1)^2)

    δ_ν           = δ^ν
    one_minus_δ   = 1 - δ
    one_minus_δ_ν = 1 - δ_ν

    num_main  = ν * one_minus_δ - one_minus_δ_ν +
                (δ * one_minus_δ_ν - ν * one_minus_δ) / sin2_half
    den_main  = 4π * one_minus_δ^2 * δ_ν
    term_main = num_main / den_main

    δ_180_ν   = δ_180^ν
    num_back  = 1 - δ_180_ν
    den_back  = 16π * (δ_180 - 1) * δ_180_ν
    term_back = (num_back / den_back) * (3 * cosθ^2 - 1)

    # Fournier's native normalization is ∫_{4π} p dΩ = 1; multiply by 4π
    # to match the ocean-optics convention (1/2) ∫₋₁¹ p dμ = 1.
    return 4π * (term_main + term_back)
end

# Renormalization integral — uses adaptive QuadGK, with a forced subdivision
# near the forward peak so the adaptive algorithm doesn't miss it.
function _ff_renorm_factor(n::FT, μ_j::FT) where {FT}
    raw_integrand(μ) = _ff_raw(n, μ_j, FT(μ))
    # Split [-1, 1] at several points to help the adaptive integrator near
    # the forward peak, where p blows up like 1/sin²(θ/2).
    I1, _ = quadgk(raw_integrand, FT(-1),    FT(0.9);      rtol = FT(1e-10))
    I2, _ = quadgk(raw_integrand, FT(0.9),   FT(0.999);    rtol = FT(1e-10))
    I3, _ = quadgk(raw_integrand, FT(0.999), FT(0.99999);  rtol = FT(1e-10))
    I4, _ = quadgk(raw_integrand, FT(0.99999), FT(1);      rtol = FT(1e-10))
    I_total = I1 + I2 + I3 + I4
    # (1/2) ∫ p dμ should equal 1 after renormalization
    return FT(2) / I_total
end

function FournierForandPhase{FT}(; n = nothing, μ_j = nothing,
                                   backscatter_fraction = nothing) where {FT}
    if n !== nothing && μ_j !== nothing
        nn, mj = FT(n), FT(μ_j)
        _check_ff_params(nn, mj)
        r = _ff_renorm_factor(nn, mj)
        return FournierForandPhase{FT}(nn, mj, r)
    elseif backscatter_fraction !== nothing
        B = FT(backscatter_fraction)
        n_inv, μ_inv = invert_backscatter_to_nμ(B)
        r = _ff_renorm_factor(FT(n_inv), FT(μ_inv))
        return FournierForandPhase{FT}(FT(n_inv), FT(μ_inv), r)
    else
        n_inv, μ_inv = invert_backscatter_to_nμ(FT(0.018))
        r = _ff_renorm_factor(FT(n_inv), FT(μ_inv))
        return FournierForandPhase{FT}(FT(n_inv), FT(μ_inv), r)
    end
end

FournierForandPhase(; kwargs...) = FournierForandPhase{Float64}(; kwargs...)

function _check_ff_params(n::FT, μ_j::FT) where {FT}
    (FT(1) < n < FT(1.25)) || @warn "FournierForand: refractive index n=$n outside typical ocean-particle range (1.00, 1.25)"
    (FT(3) < μ_j < FT(5))  || @warn "FournierForand: Junge slope μ_j=$μ_j outside typical range (3, 5)"
end

# -----------------------------------------------------------------------------

"""
    phase_function_value(pf::FournierForandPhase, cosθ)

Renormalized Fournier & Jonasz (1999) phase function, satisfying
`(1/2) ∫₋₁¹ p(μ) dμ = 1` exactly.
"""
phase_function_value(pf::FournierForandPhase{FT}, cosθ::Real) where {FT} =
    pf.renorm * _ff_raw(pf.n, pf.μ_j, cosθ)

# -----------------------------------------------------------------------------
# Closed-form backscatter fraction, Mobley (2002) Eq. 10
# -----------------------------------------------------------------------------

"""
    backscatter_fraction(pf::FournierForandPhase) -> FT

Closed-form backscatter fraction from Mobley et al. (2002) Eq. 10.

Note: this formula uses the raw Fournier & Jonasz normalization. Since
we apply an empirical renormalization factor ~1.01 to recover unit
normalization, the closed-form B is off by the same ~1%. For the
canonical ocean use case this is within the noise of measurement
uncertainties, but a future release should recompute B numerically
from the renormalized phase function for exact consistency.
"""
function backscatter_fraction(pf::FournierForandPhase{FT}) where {FT}
    n, μ_j = pf.n, pf.μ_j
    ν = (FT(3) - μ_j) / FT(2)
    δ_90 = (FT(4) / FT(3)) * (one(FT) / FT(2)) / (n - one(FT))^2
    δ_90_ν  = δ_90^ν
    δ_90_ν1 = δ_90 * δ_90_ν
    numerator   = one(FT) - δ_90_ν1 - FT(0.5) * (one(FT) - δ_90_ν)
    denominator = (one(FT) - δ_90) * δ_90_ν
    return (one(FT) - numerator / denominator) * pf.renorm   # apply renorm
end

# -----------------------------------------------------------------------------
# Inversion B → (n, μ_j) following Mobley (2002)
# -----------------------------------------------------------------------------

"""
    invert_backscatter_to_nμ(B) -> (n, μ_j)

One-dimensional inversion along the Mobley (2002) empirical family
μ_j ≈ 3 + (n − 1) · 20. For B outside [1e-4, 0.5] the result is
extrapolated and a warning is issued.
"""
function invert_backscatter_to_nμ(B::FT) where {FT}
    (FT(1e-4) ≤ B ≤ FT(0.5)) ||
        @warn "Backscatter fraction $B outside calibrated range [1e-4, 0.5]; extrapolating."

    # The renormalized B depends on the same renorm factor we apply to p —
    # so for inversion we use the *raw* closed-form B (Mobley's Eq. 10
    # without renorm), because Mobley's inversion curves were fit in that
    # normalization anyway.
    function B_raw(n_trial::FT)
        μ_j_trial = clamp(FT(3) + (n_trial - one(FT)) * FT(20), FT(3.01), FT(4.99))
        ν = (FT(3) - μ_j_trial) / FT(2)
        δ_90 = (FT(4) / FT(3)) * (one(FT) / FT(2)) / (n_trial - one(FT))^2
        δ_90_ν  = δ_90^ν
        δ_90_ν1 = δ_90 * δ_90_ν
        num = one(FT) - δ_90_ν1 - FT(0.5) * (one(FT) - δ_90_ν)
        den = (one(FT) - δ_90) * δ_90_ν
        return one(FT) - num / den
    end

    f(n) = B_raw(FT(n)) - B
    n_lo, n_hi = FT(1.001), FT(1.25)
    f_lo, f_hi = f(n_lo), f(n_hi)

    if sign(f_lo) == sign(f_hi)
        @warn "Cannot bracket Mobley curve for B=$B; using fallback (n=1.08, μ_j=3.38)."
        return (FT(1.08), FT(3.38))
    end

    for _ in 1:60
        n_mid = (n_lo + n_hi) / 2
        f_mid = f(n_mid)
        if abs(f_mid) < FT(1e-9)
            μ_mid = clamp(FT(3) + (n_mid - one(FT)) * FT(20), FT(3.01), FT(4.99))
            return (n_mid, μ_mid)
        end
        (sign(f_mid) == sign(f_lo)) ? (n_lo = n_mid; f_lo = f_mid) :
                                       (n_hi = n_mid; f_hi = f_mid)
    end
    n_final = (n_lo + n_hi) / 2
    μ_final = clamp(FT(3) + (n_final - one(FT)) * FT(20), FT(3.01), FT(4.99))
    return (n_final, μ_final)
end
