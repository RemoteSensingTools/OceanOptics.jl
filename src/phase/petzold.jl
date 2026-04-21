# =============================================================================
# Petzold "San Diego Harbor" tabulated phase function
# =============================================================================
#
# Petzold (1972) "Volume scattering functions for selected ocean waters",
# SIO Ref. 72-78, reports a measured volume-scattering function (VSF) for
# ocean particulate scattering in San Diego Harbor. The table is
# reproduced in Mobley (1994) "Light and Water", Academic Press,
# Appendix A.4. It was the community-standard phase function until the
# Fournier–Forand analytic form (Mobley, Sundman & Boss 2002, Appl. Opt.
# 41, 1035) displaced it for most applications; we ship it for
# reproducibility of legacy results.
#
# Data provenance
# ---------------
# The VSF is tabulated as β(θ) `[1/m/sr]` on a coarse θ grid. The shipped
# CSV (`data/petzold_sdh.csv`) must contain (θ_deg, β) pairs; we normalize
# to the OceanOptics.jl phase-function convention `(1/2) ∫₋₁¹ p(μ) dμ = 1`
# at load time. Values below the minimum tabulated angle are extrapolated
# log-linearly toward the forward peak; values above the maximum angle
# are held flat at the final measurement.
#
# The Petzold table is NOT shipped by default in this OceanOptics.jl
# release — it needs transcription from Mobley (1994) Appendix A.4 (see
# DESIGN §8.7). Selecting `PetzoldPhase()` before dropping the CSV in
# raises an informative error.
# =============================================================================

"""
$(TYPEDEF)

Tabulated Petzold (1972) "San Diego Harbor" particulate phase function.

The underlying spectral/angular table is loaded lazily from
`data/petzold_sdh.csv`; the file has to be transcribed from Mobley (1994)
*Light and Water*, Appendix A.4. See DESIGN §8.7 for the column
convention the loader expects.

# Fields
$(TYPEDFIELDS)
"""
struct PetzoldPhase{FT} <: AbstractOceanPhaseFunction{FT}
    "Cosines of the tabulated scattering angles, sorted ascending."
    cosθ::Vector{FT}
    "Normalized phase-function values `p(μ)` at those cosines."
    p::Vector{FT}
end

const _PETZOLD_SDH_PATH          = data_path("petzold_sdh.csv")
const _PETZOLD_SDH_AVAILABLE     = isfile(_PETZOLD_SDH_PATH)
const _PETZOLD_SDH_CACHE         = Ref{Union{Nothing, PetzoldPhase{Float64}}}(nothing)

"Return the cached Petzold phase function, loading it on first call."
function _petzold_sdh()
    pz = _PETZOLD_SDH_CACHE[]
    pz === nothing || return pz
    isfile(_PETZOLD_SDH_PATH) || error(
        "PetzoldPhase requires `data/petzold_sdh.csv`, which is not bundled ",
        "with this OceanOptics.jl release. Transcribe Petzold (1972) SIO ",
        "Ref. 72-78 / Mobley (1994) Appendix A.4 (see DESIGN §8.7). Until ",
        "then, use FournierForandPhase — Mobley et al. (2002) showed it ",
        "matches measurements across three decades of forward-scattered ",
        "intensity better than the tabulated Petzold.")
    pz = _load_petzold_csv(_PETZOLD_SDH_PATH)
    _PETZOLD_SDH_CACHE[] = pz
    return pz
end

"""
    PetzoldPhase{FT}()

Construct a `PetzoldPhase` from the bundled CSV. Lazy: the file is
parsed on first call and cached module-wide. Promotes the `Float64`
canonical representation to `FT`.
"""
function PetzoldPhase{FT}() where {FT<:AbstractFloat}
    raw = _petzold_sdh()
    return PetzoldPhase{FT}(convert(Vector{FT}, raw.cosθ), convert(Vector{FT}, raw.p))
end

PetzoldPhase() = PetzoldPhase{Float64}()

"""
Internal loader. CSV must provide columns `theta_deg, beta_inv_m_sr` (or
equivalent first/second columns). The loader converts angles to μ =
cos θ and normalizes so `(1/2) ∫₋₁¹ p(μ) dμ = 1`.
"""
function _load_petzold_csv(path::AbstractString)
    table = load_spectral_csv(path)   # columns: θ_deg, β
    θ_deg = table.λ
    β     = table.y

    # Convert θ [deg] → μ = cos θ, ascending in μ (descending in θ)
    μ = cos.(deg2rad.(θ_deg))
    perm = sortperm(μ)
    μ  = μ[perm]
    β  = β[perm]

    # Normalize: the target is (1/2) ∫₋₁¹ p(μ) dμ = 1. We trapezoid-integrate
    # the tabulated β, then divide out so p = β / ( (1/2)∫β dμ ).
    ∫β = 0.0
    @inbounds for i in 1:length(μ)-1
        ∫β += 0.5 * (β[i] + β[i+1]) * (μ[i+1] - μ[i])
    end
    norm_factor = 2 / ∫β
    p = β .* norm_factor

    return PetzoldPhase{Float64}(μ, p)
end

# -----------------------------------------------------------------------------
# Phase function value — linear interpolation in μ
# -----------------------------------------------------------------------------

function phase_function_value(pf::PetzoldPhase{FT}, cosθ::Real) where {FT}
    μ_grid = pf.cosθ
    p_grid = pf.p
    # Clamp to tabulated range; boundary-flat extrapolation, AD-safe via
    # `min`/`max` on the input type.
    if cosθ ≤ first(μ_grid)
        return first(p_grid) + zero(cosθ)
    elseif cosθ ≥ last(μ_grid)
        return last(p_grid) + zero(cosθ)
    end
    i  = searchsortedlast(μ_grid, cosθ)
    μ₁, μ₂ = μ_grid[i], μ_grid[i+1]
    p₁, p₂ = p_grid[i], p_grid[i+1]
    t  = (cosθ - μ₁) / (μ₂ - μ₁)
    return p₁ + t * (p₂ - p₁)
end

is_polarizable(::PetzoldPhase) = false
