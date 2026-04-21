# =============================================================================
# Phytoplankton absorption and scattering
# =============================================================================
#
# The phytoplankton absorption spectrum a_φ(λ) is the single biggest
# moving piece of ocean optics in the visible. The canonical empirical
# parameterization is Bricaud's power-law:
#
#     a_φ*(λ) = A(λ) · Chl^(-E(λ))    [m²/mg]
#     a_φ(λ)  = a_φ*(λ) · Chl         [1/m]
#
# where `A(λ)` and `E(λ)` are spectral coefficients fit to globally
# distributed HPLC + spectrophotometer samples. Two generations ship:
# the original Bricaud et al. 1995 JGR fit (visible only; NOT bundled),
# and the 1998 update (Bricaud, Morel, Babin, Allali, Claustre JGR 103,
# 31033) extended to UV by Morrison & Nelson (2004) and Vasilkov et al.
# (2005). The 1998 update is the community-standard replacement and is
# the version bundled here as `data/bricaud_1998_{A,E}.csv`.
#
# Particulate scattering uses the Morel & Maritorena (2001, JGR 106,
# 7163) Case-1 form:
#
#     b_p(λ) = 0.30 · Chl^0.62 · (550/λ)          [1/m]
#
# The angular distribution is Fournier-Forand with a typical open-ocean
# backscatter fraction B ≈ 0.01; the caller can override at construction.
# =============================================================================

"""
$(TYPEDEF)

Bricaud et al. (1998, JGR 103, 31033) power-law phytoplankton absorption,
extended to UV by Morrison & Nelson (2004) and Vasilkov et al. (2005).
Bundled as `data/bricaud_1998_A.csv` and `data/bricaud_1998_E.csv`,
covering 300–1000 nm at 2 nm resolution.
"""
struct Bricaud1998 <: AbstractPhytoplanktonModel end

"""
$(TYPEDEF)

Gordon (1992) broadband-averaged phytoplankton absorption (simpler,
less accurate, but historical baseline). Reserved as a model tag;
implementation stubbed for future release.
"""
struct Gordon1992 <: AbstractPhytoplanktonModel end

"""
$(TYPEDEF)

Phytoplankton (living algal cells) as an absorbing scatterer.

`Chl` is the near-surface chlorophyll-a concentration `[mg/m³]`. The
model tag `M` selects the absorption parameterization; scattering
currently uses Morel & Maritorena (2001) irrespective of `M`, and the
phase function defaults to `FournierForandPhase` with backscatter
fraction 0.01 (typical open-ocean living phytoplankton).

# Fields
$(TYPEDFIELDS)

# Examples

```julia
# Open-ocean Case-1 water at Chl = 0.5 mg/m³
phyto = Phytoplankton{Float64, Bricaud1998}(Chl = 0.5)

# Custom backscatter fraction (coastal, larger cells)
phyto = Phytoplankton{Float64, Bricaud1998}(Chl = 2.0, B = 0.015)
```
"""
struct Phytoplankton{FT, M<:AbstractPhytoplanktonModel,
                        P<:AbstractOceanPhaseFunction{FT}} <: AbstractAbsorbingScatterer{FT}
    "Parameterization-model tag (e.g. `Bricaud1998()`)"
    model::M
    "Chlorophyll-a concentration `[mg/m³]`"
    Chl::FT
    "Particulate-scattering phase function"
    phase::P
end

"""
    Phytoplankton{FT, M}(; Chl, B=0.01)

Construct a `Phytoplankton` instance with Fournier-Forand phase function
parameterized by backscatter fraction `B`.
"""
function Phytoplankton{FT, M}(; Chl, B = FT(0.01)) where {FT<:AbstractFloat, M<:AbstractPhytoplanktonModel}
    Chl ≥ zero(Chl) ||
        throw(ArgumentError("Phytoplankton: Chl must be non-negative, got $Chl"))
    phase = FournierForandPhase{FT}(backscatter_fraction = FT(B))
    return Phytoplankton{FT, M, typeof(phase)}(M(), FT(Chl), phase)
end

Phytoplankton{FT}(; kwargs...) where {FT<:AbstractFloat} =
    Phytoplankton{FT, Bricaud1998}(; kwargs...)
Phytoplankton(; kwargs...) = Phytoplankton{Float64}(; kwargs...)

phase_function(p::Phytoplankton) = p.phase

# -----------------------------------------------------------------------------
# Absorption — Bricaud (1998) + UV extensions
# -----------------------------------------------------------------------------

"A(λ) spectral pre-factor, Bricaud 1998 + UV extensions, 300–1000 nm."
const _BRICAUD_1998_A_TABLE = load_spectral_csv(data_path("bricaud_1998_A.csv"))

"E(λ) spectral exponent, Bricaud 1998 + UV extensions, 300–1000 nm."
const _BRICAUD_1998_E_TABLE = load_spectral_csv(data_path("bricaud_1998_E.csv"))

function absorption(p::Phytoplankton{FT, Bricaud1998}, λ::Real) where {FT}
    # a_φ*(λ) = A(λ) · Chl^{-E(λ)}, a_φ(λ) = a_φ*(λ) · Chl = A · Chl^{1-E}
    A_λ = linterp(_BRICAUD_1998_A_TABLE, λ)
    E_λ = linterp(_BRICAUD_1998_E_TABLE, λ)
    return A_λ * p.Chl^(one(FT) - E_λ)
end

function absorption(p::Phytoplankton{FT, Gordon1992}, λ::Real) where {FT}
    error("Phytoplankton{Gordon1992} absorption is not yet implemented. ",
          "Use Phytoplankton{$(nameof(FT)), Bricaud1998} with the bundled ",
          "Bricaud 1998 + UV-extended tables, or wait for a future release.")
end

absorption(p::Phytoplankton, λ::AbstractVector) = [absorption(p, λi) for λi in λ]

# -----------------------------------------------------------------------------
# Scattering — Morel & Maritorena (2001) Case-1 particulate form
# -----------------------------------------------------------------------------
#
#   b_p(λ) = 0.30 · Chl^0.62 · (550/λ)        [1/m]
#
# Integrated over living-cell size distribution; intended for open-ocean
# Case-1 waters. In coastal / turbid waters this can underestimate the
# particulate scattering by a factor of 2 or more and should be paired
# with an explicit `NonAlgalParticles` constituent.

const _MORELMARITORENA_COEFF    = 0.30
const _MORELMARITORENA_EXPONENT = 0.62
const _MORELMARITORENA_λ0       = 550.0       # [nm]

function scattering(p::Phytoplankton{FT}, λ::Real) where {FT}
    # All-literal arithmetic promotes correctly for `ForwardDiff.Dual` inputs.
    return _MORELMARITORENA_COEFF * p.Chl^_MORELMARITORENA_EXPONENT *
           (_MORELMARITORENA_λ0 / λ)
end

scattering(p::Phytoplankton, λ::AbstractVector) = [scattering(p, λi) for λi in λ]

# Backscatter — phase-function-determined. Delegates to the configured
# Fournier-Forand's closed-form `backscatter_fraction`.
function backscattering(p::Phytoplankton{FT}, λ::Real) where {FT}
    return scattering(p, λ) * backscatter_fraction(p.phase)
end

backscattering(p::Phytoplankton, λ::AbstractVector) = [backscattering(p, λi) for λi in λ]

function Base.show(io::IO, p::Phytoplankton{FT, M}) where {FT, M}
    print(io, "Phytoplankton{$FT}(", nameof(M), ", Chl=$(p.Chl) mg/m³, ",
              "B=", backscatter_fraction(p.phase), ")")
end
