# =============================================================================
# Angular quadrature and generalized-spherical-function expansions
# =============================================================================
#
# Expansion integrals are evaluated in the scattering angle θ, not in
# μ = cos θ. Ocean particle phase functions are strongly forward peaked —
# Fournier-Forand diverges like θ^(2|ν|-2) — and an appreciable part of
# their scattering lies at angles far below what μ can resolve near 1 in
# floating point (for B = 0.001, 0.25 % of it lies within 1e-8 rad).
#
# The rule is composite Gauss-Legendre in θ: uniform panels wide enough to
# resolve P_ℓ for ℓ ≤ ℓ_max, the first panel replaced by geometrically
# graded panels down to 1e-100 rad (hp-quadrature for the integrable
# forward singularity), and optional extra break points (table nodes,
# π/2). Results are refined until successive levels agree.
#
# Expansions follow Mishchenko, Travis & Lacis (2002) §4.11 with Wigner
# d-functions d^s_{mn}(θ); the coefficients are then mapped to the
# vSmartMOM / Hovenier convention of `GreekCoefficients`:
#
#   β = α₁,  α = α₂,  ζ = α₃,  δ = α₄,  γ = −β₁,  ϵ = −β₂
#
# — the same sign flips vSmartMOM applies to VLIDORT-format coefficients.
# Angular scattering matrices are in the standard convention of Hansen &
# Travis (1974) and Mishchenko et al. (2002): F₁₂ < 0 for Rayleigh
# scattering at 90°.
# =============================================================================

const _GEOMETRIC_RATIO = 0.15
const _θ_FLOOR         = 1e-100

"""
    _angular_rule(ℓ_max, level; breaks = (), forward_singular = true) -> (θ, w, θ_lo)

Nodes `θ ∈ (θ_lo, π)` and weights such that `∫_{θ_lo}^π f(θ) sin θ dθ ≈
Σ wᵢ f(θᵢ)`, and the lower bound `θ_lo` (below `1e-100` with
`forward_singular`, else 0). `level` doubles the panel count and adds nodes
per panel. `breaks` are additional panel boundaries in radians (kinks of
tabulated functions).
"""
function _angular_rule(ℓ_max::Integer, level::Integer; breaks = (),
                       forward_singular::Bool = true)
    n_panel = 16 + 8level
    n_geo   = 10 + 4level
    width   = π / (max(16, cld(ℓ_max, 4)) * 2^level)

    edges = sort!(unique!(Float64[0.0; π / 2; π; [b for b in breaks if 0 < b < π]...]))
    bounds = Float64[]
    for i in 1:length(edges)-1
        a, b = edges[i], edges[i+1]
        m = max(1, ceil(Int, (b - a) / width))
        append!(bounds, range(a, b; length = m + 1)[1:end-1])
    end
    push!(bounds, π)

    x, wx = gausslegendre(n_panel)
    θ, w  = Float64[], Float64[]
    first_upper = bounds[2]
    θ_lo = 0.0
    if forward_singular
        # Geometric panels [h σ^(k+1), h σ^k] replace the first panel.
        xg, wg = gausslegendre(n_geo)
        hi = first_upper
        while hi > _θ_FLOOR
            lo = hi * _GEOMETRIC_RATIO
            _append_panel!(θ, w, xg, wg, lo, hi)
            hi = lo
        end
        θ_lo = hi
    else
        _append_panel!(θ, w, x, wx, 0.0, first_upper)
    end
    for i in 2:length(bounds)-1
        _append_panel!(θ, w, x, wx, bounds[i], bounds[i+1])
    end
    return θ, w, θ_lo
end

function _append_panel!(θ, w, x, wx, a, b)
    half, mid = (b - a) / 2, (b + a) / 2
    for (xi, wi) in zip(x, wx)
        t = mid + half * xi
        push!(θ, t)
        push!(w, half * wi * sin(t))
    end
    return nothing
end

"""
    _converged(f, ℓ_max; tol = 1e-10, levels = 0:4)

Evaluate `f(level)` (returning a vector of normalized moments) at
increasing refinement until two successive levels agree to `tol`; throws
if they never do.
"""
function _converged(f, ℓ_max::Integer; tol = 1e-10, levels = 0:4, what = "expansion")
    prev = f(first(levels))
    for level in levels[2:end]
        cur  = f(level)
        diff = maximum(abs, _scaled(cur, ℓ_max) .- _scaled(prev, ℓ_max))
        diff ≤ tol && return cur
        prev = cur
    end
    throw(ErrorException("$what did not converge to $tol within $(length(levels)) " *
                         "refinement levels (ℓ_max = $ℓ_max)"))
end

# Moments divided by (2ℓ+1), so all lie in [-1, 1] for physical functions.
_scaled(v::AbstractVector, ℓ_max) = v ./ (2 .* (0:length(v)-1) .+ 1)
_scaled(v::Tuple, ℓ_max) = reduce(vcat, (_scaled(x, ℓ_max) for x in v))

# -----------------------------------------------------------------------------
# Legendre moments of a scalar function of θ
# -----------------------------------------------------------------------------

"""
    _legendre_moments(fvals, θ, w, ℓ_max) -> Vector{Float64}

`β_ℓ = (2ℓ+1)/2 ∫₀^π f(θ) P_ℓ(cos θ) sin θ dθ` for `ℓ = 0:ℓ_max`.
"""
function _legendre_moments(fvals::AbstractVector, θ::AbstractVector, w::AbstractVector,
                           ℓ_max::Integer)
    β = zeros(ℓ_max + 1)
    @inbounds for i in eachindex(θ)
        x  = cos(θ[i])
        fw = fvals[i] * w[i]
        P_prev, P = 0.0, 1.0
        β[1] += fw
        for ℓ in 1:ℓ_max
            P_prev, P = P, ((2ℓ - 1) * x * P - (ℓ - 1) * P_prev) / ℓ
            β[ℓ+1] += fw * P
        end
    end
    for ℓ in 0:ℓ_max
        β[ℓ+1] *= (2ℓ + 1) / 2
    end
    return β
end

# -----------------------------------------------------------------------------
# Wigner d-functions d^s_{mn}(θ) for (m, n) ∈ {(0,2), (2,2), (2,-2)}
# -----------------------------------------------------------------------------
#
# Recurrence (Mishchenko et al. 2002, Eq. B.22):
#   d^{s+1}_{mn} = [(2s+1)(s(s+1)x − mn) d^s_{mn}
#                   − (s+1)√(s²−m²)√(s²−n²) d^{s−1}_{mn}]
#                  / [s √((s+1)²−m²) √((s+1)²−n²)]
# from d^{1} = 0 and d^2_{02} = (√6/4) sin²θ, d^2_{22} = cos⁴(θ/2),
# d^2_{2,−2} = sin⁴(θ/2), written with half angles for accuracy near 0, π.

@inline function _wigner_start(θ)
    s2, c2 = sin(θ / 2)^2, cos(θ / 2)^2
    return sqrt(6.0) * s2 * c2, c2^2, s2^2          # d02, d22, d2-2 at s = 2
end

@inline function _wigner_step(s, x, mn, m2, n2, d, d_prev)
    a = (2s + 1) * (s * (s + 1) * x - mn) * d
    b = s ≥ 2 ? (s + 1) * sqrt((s^2 - m2) * (s^2 - n2)) * d_prev : 0.0
    return (a - b) / (s * sqrt(((s + 1)^2 - m2) * ((s + 1)^2 - n2)))
end

"""
    _mueller_moments(F, θ, w, ℓ_max) -> NTuple{6, Vector{Float64}}

Mishchenko-convention coefficients `(α₁, α₂, α₃, α₄, β₁, β₂)` of the
scattering matrix with elements `F = (F11, F12, F22, F33, F34, F44)`
sampled at the nodes `θ`.
"""
function _mueller_moments(F, θ::AbstractVector, w::AbstractVector, ℓ_max::Integer)
    α₁ = zeros(ℓ_max + 1); α₄ = zeros(ℓ_max + 1)
    β₁ = zeros(ℓ_max + 1); β₂ = zeros(ℓ_max + 1)
    sp = zeros(ℓ_max + 1); sm = zeros(ℓ_max + 1)       # α₂ ± α₃
    @inbounds for i in eachindex(θ)
        x  = cos(θ[i]); wi = w[i]
        f11, f12, f22, f33, f34, f44 = (F.F11[i], F.F12[i], F.F22[i],
                                        F.F33[i], F.F34[i], F.F44[i])
        P_prev, P = 0.0, 1.0
        α₁[1] += wi * f11; α₄[1] += wi * f44
        if ℓ_max ≥ 1
            P_prev, P = P, x
            α₁[2] += wi * f11 * P; α₄[2] += wi * f44 * P
        end
        d02, d22, d2m2 = _wigner_start(θ[i])
        d02p = d22p = d2m2p = 0.0
        for s in 2:ℓ_max
            P_prev, P = P, ((2s - 1) * x * P - (s - 1) * P_prev) / s
            α₁[s+1] += wi * f11 * P;          α₄[s+1] += wi * f44 * P
            β₁[s+1] += wi * f12 * d02;        β₂[s+1] += wi * f34 * d02
            sp[s+1] += wi * (f22 + f33) * d22; sm[s+1] += wi * (f22 - f33) * d2m2
            if s < ℓ_max
                d02,  d02p  = _wigner_step(s, x,  0, 0, 4, d02,  d02p),  d02
                d22,  d22p  = _wigner_step(s, x,  4, 4, 4, d22,  d22p),  d22
                d2m2, d2m2p = _wigner_step(s, x, -4, 4, 4, d2m2, d2m2p), d2m2
            end
        end
    end
    c  = [(2s + 1) / 2 for s in 0:ℓ_max]
    α₂ = c .* (sp .+ sm) ./ 2
    α₃ = c .* (sp .- sm) ./ 2
    return c .* α₁, α₂, α₃, c .* α₄, c .* β₁, c .* β₂
end

# Mishchenko (α₁, α₂, α₃, α₄, β₁, β₂) → GreekCoefficients (vSmartMOM convention).
function _greek_from_mishchenko(::Type{FT}, α₁, α₂, α₃, α₄, β₁, β₂) where {FT}
    return GreekCoefficients{FT}(FT.(α₂), FT.(α₁), FT.(.-β₁), FT.(α₄), FT.(.-β₂), FT.(α₃))
end

"""
    evaluate_scattering_matrix(g::GreekCoefficients, μ) -> NamedTuple

Angular scattering matrix reconstructed from the expansion at scattering
cosines `μ`, in the standard convention (Hansen & Travis 1974;
Mishchenko et al. 2002: `F₁₂ < 0` for Rayleigh scattering at 90°):

```
F = [F11 F12  0    0
     F12 F22  0    0
     0   0    F33  F34
     0   0   -F34  F44]
```

normalized so that `(1/2)∫ F11 dμ = 1`. Returns `(; F11)` for a scalar
expansion and `(; F11, F12, F22, F33, F34, F44)` otherwise.
"""
function evaluate_scattering_matrix(g::GreekCoefficients{FT}, μ::AbstractVector{<:Real}) where {FT}
    L   = n_moments(g) - 1
    pol = is_polarized(g)
    out = [zeros(FT, length(μ)) for _ in 1:6]
    for (i, μi) in enumerate(μ)
        x = clamp(μi, -1, 1)
        θ = acos(x)
        P_prev, P = zero(FT), one(FT)
        f11 = g.β[1]; f44 = pol ? g.δ[1] : zero(FT)
        if L ≥ 1
            P_prev, P = P, FT(x)
            f11 += g.β[2] * P; pol && (f44 += g.δ[2] * P)
        end
        f12 = f34 = sp = sm = zero(FT)
        d02, d22, d2m2 = _wigner_start(θ)
        d02p = d22p = d2m2p = 0.0
        for s in 2:L
            P_prev, P = P, ((2s - 1) * x * P - (s - 1) * P_prev) / s
            f11 += g.β[s+1] * P
            if pol
                f44 += g.δ[s+1] * P
                f12 -= g.γ[s+1] * d02;  f34 -= g.ϵ[s+1] * d02
                sp  += (g.α[s+1] + g.ζ[s+1]) * d22
                sm  += (g.α[s+1] - g.ζ[s+1]) * d2m2
            end
            if s < L
                d02,  d02p  = _wigner_step(s, x,  0, 0, 4, d02,  d02p),  d02
                d22,  d22p  = _wigner_step(s, x,  4, 4, 4, d22,  d22p),  d22
                d2m2, d2m2p = _wigner_step(s, x, -4, 4, 4, d2m2, d2m2p), d2m2
            end
        end
        out[1][i] = f11; out[2][i] = f12; out[3][i] = (sp + sm) / 2
        out[4][i] = (sp - sm) / 2; out[5][i] = f34; out[6][i] = f44
    end
    pol || return (; F11 = out[1])
    return (; F11 = out[1], F12 = out[2], F22 = out[3], F33 = out[4], F34 = out[5], F44 = out[6])
end
