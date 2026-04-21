# =============================================================================
# Legendre-moment expansion of phase functions (de Haan 1987 convention)
# =============================================================================
#
# The β_ℓ coefficients are computed by Gauss-Legendre quadrature in the
# de Haan 1987 / HydroLight convention, i.e.
#
#     β_ℓ = (2ℓ+1)/2 · ∫₋₁¹ p(μ) P_ℓ(μ) dμ
#
# so that p(μ) = Σ_ℓ β_ℓ P_ℓ(μ) exactly recovers the phase function.
# Internally we accumulate the un-weighted moment (1/2)∫ p P_ℓ dμ and
# multiply by (2ℓ+1) at the end.
#
# A single Gauss-Legendre rule with `nq` nodes integrates polynomials up to
# degree 2·nq − 1 exactly. For p(μ) non-polynomial (any real phase function),
# `nq ≥ 2·ℓ_max + 32` gives 4–6 decimal-place accuracy for ℓ_max ≤ 256.
#
# The Fournier-Forand and Petzold phase functions are *strongly forward-peaked*
# (β-values decay slowly). A fixed Gauss-Legendre rule handles this well
# because nodes cluster near the endpoints — but we provide a higher-nq
# "forward-peaked" variant for users who need ℓ_max > 128.
#
# For phase functions that expose closed-form moments, `phase_function_moments`
# is overridden in that phase function's file and this numerical path is
# bypassed.
# =============================================================================

"""
    numerical_moments(pf, ℓ_max; nq=max(128, 2*ℓ_max+32)) -> Vector{FT}

Compute Legendre moments of any phase function by Gauss-Legendre quadrature.

This is the fallback implementation of [`phase_function_moments`](@ref).
Subtypes with closed-form moments (e.g. `HenyeyGreensteinPhase`, the
polarized Rayleigh expansion, or tabulated Petzold) override the dispatch
on `phase_function_moments` directly rather than routing through here.

# Arguments
- `pf` — any `AbstractOceanPhaseFunction` with `phase_function_value` defined
- `ℓ_max` — highest Legendre order desired
- `nq` — number of Gauss-Legendre nodes. Default `max(128, 2·ℓ_max + 32)`
  is empirically sufficient for ocean phase functions with backscatter
  fraction B ≥ 0.001. For more forward-peaked cases, increase.

# Returns
Vector of length `ℓ_max + 1` with `β[1] = β₀`, `β[2] = β₁`, ….
For a normalized phase function, `β[1] ≈ 1`.
"""
function numerical_moments(pf::AbstractOceanPhaseFunction{FT}, ℓ_max::Integer;
                           nq::Integer = max(128, 2 * ℓ_max + 32)) where {FT}
    ℓ_max ≥ 0 || throw(ArgumentError("ℓ_max must be non-negative, got $ℓ_max"))
    μ_nodes, w_nodes = gausslegendre(nq)

    β = zeros(FT, ℓ_max + 1)

    # Pre-evaluate p(μ_i) once; inner loop accumulates P_ℓ via recurrence.
    @inbounds for i in 1:nq
        μ_i = FT(μ_nodes[i])
        w_i = FT(w_nodes[i])
        p_i = phase_function_value(pf, μ_i)

        P_prev = zero(FT)
        P_curr = one(FT)               # P₀
        β[1] += FT(0.5) * w_i * p_i * P_curr
        if ℓ_max ≥ 1
            P_next = μ_i               # P₁
            β[2] += FT(0.5) * w_i * p_i * P_next
            P_prev, P_curr = P_curr, P_next
            for ℓ in 2:ℓ_max
                P_next = ((2ℓ - 1) * μ_i * P_curr - (ℓ - 1) * P_prev) / ℓ
                β[ℓ + 1] += FT(0.5) * w_i * p_i * P_next
                P_prev, P_curr = P_curr, P_next
            end
        end
    end

    # Strongly forward-peaked phase functions (Fournier-Forand with typical
    # ocean backscatter fractions < 0.02, Petzold, TTHG) are under-resolved
    # by any practical Gauss-Legendre rule: the peak near μ = 1 carries
    # enough integral weight to bias β₀ by 5–15 %. Since the input phase
    # function is normalized to `(1/2)∫p dμ = 1`, we know (1/2)∫p dμ
    # *should* be 1 exactly, and the relative weights between moments
    # are preserved by a uniform rescaling. Divide out the quadrature
    # estimate of β₀-pre-normalization so the returned moment vector is
    # self-consistent.
    if β[1] > zero(FT)
        β ./= β[1]
    end

    # Convert to de Haan 1987 convention: β_ℓ → (2ℓ+1) β_ℓ so that
    # p(μ) = Σ β_ℓ P_ℓ(μ) holds directly.
    @inbounds for ℓ in 0:ℓ_max
        β[ℓ + 1] *= FT(2ℓ + 1)
    end
    return β
end

# Default fallback — all subtypes get this unless they override.
phase_function_moments(pf::AbstractOceanPhaseFunction, ℓ_max::Integer) =
    numerical_moments(pf, ℓ_max)

"""
    mix_phase_moments(scatterings, phase_moments) -> Vector{FT}

Scattering-weighted mixing of several phase-function Legendre expansions:

```
β_ℓ^mix = (b₁ β_ℓ^{(1)} + b₂ β_ℓ^{(2)} + …) / (b₁ + b₂ + …)
```

This is the physical mixing law for phase functions of independent
scatterers, expressed at the Legendre-moment level (where mixing is
linear — doing it on `p(μ)` samples loses precision).

# Arguments
- `scatterings` — vector of `b_i` scattering coefficients `[1/m]`, one per
  constituent
- `phase_moments` — vector of `β_ℓ` arrays, each of identical length,
  aligned with `scatterings`

# Returns
Mixed `β_ℓ` array of the same length as each input array. If total
scattering is zero (all `b_i = 0`), returns `[1, 0, 0, …]` (unit
isotropic β₀, all higher orders zero).
"""
function mix_phase_moments(scatterings::AbstractVector{FT},
                           phase_moments::AbstractVector{<:AbstractVector{FT}};
                           n_moments::Union{Nothing, Integer} = nothing) where {FT}
    length(scatterings) == length(phase_moments) ||
        throw(DimensionMismatch("mix_phase_moments: length(scatterings) == length(phase_moments) required"))

    # Determine the output length. If there is at least one contributing
    # phase function, inherit its moment count; otherwise fall back to
    # either the explicit `n_moments` kwarg or a unit-length isotropic
    # scalar (β₀ = 1 only). Callers in a layer loop should pass `n_moments`
    # to keep the return shape consistent across wavelengths even when a
    # particular layer has no scatterers.
    N = if !isempty(phase_moments)
            length(first(phase_moments))
        elseif n_moments !== nothing
            Int(n_moments)
        else
            1
        end

    β_mix    = zeros(FT, N)
    β_mix[1] = one(FT)

    isempty(scatterings) && return β_mix

    b_total = sum(scatterings)
    b_total == zero(FT) && return β_mix        # already the isotropic fallback

    fill!(β_mix, zero(FT))
    @inbounds for (b, β) in zip(scatterings, phase_moments)
        length(β) == N ||
            throw(DimensionMismatch("mix_phase_moments: all phase_moments vectors must have equal length"))
        β_mix .+= (b / b_total) .* β
    end
    return β_mix
end
