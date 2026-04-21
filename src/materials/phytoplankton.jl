# =============================================================================
# Phytoplankton absorption and scattering
# =============================================================================
#
# The phytoplankton absorption spectrum a_φ(λ) is the single biggest
# moving piece of ocean optics in the visible. The canonical empirical
# parameterization is Bricaud et al. (1995, 1998):
#
#     a_φ*(λ) = A(λ) · Chl^(-E(λ))    [m²/mg]
#     a_φ(λ)  = a_φ*(λ) · Chl         [1/m]
#
# where `A(λ)` and `E(λ)` are 1 nm-resolution spectra obtained by fitting
# 815 globally distributed HPLC+spectrophotometer samples (Bricaud,
# Babin, Morel, Claustre 1995, JGR 100, 13321; extension in Bricaud,
# Morel, Babin, Allali, Claustre 1998, JGR 103, 31033).
#
# Particulate scattering uses the Morel & Maritorena (2001, JGR 106,
# 7163) Case-1 form:
#
#     b_p(λ) = 0.30 · Chl^0.62 · (550/λ)          [1/m]
#
# The angular distribution is Fournier-Forand with a typical open-ocean
# backscatter fraction B ≈ 0.01; the caller can override at construction.
#
# Data provenance
# ---------------
# The Bricaud (1995) A(λ) / E(λ) tables are NOT shipped with this release
# (they require transcription from JGR 100, 13321 Table 2; see DESIGN §8.3).
# Drop a two-column CSV `data/bricaud_1995.csv` into place to activate.
# Columns: `lambda_nm, A_m2_per_mg` and a second file `data/bricaud_1995_E.csv`
# with `lambda_nm, E`.
# =============================================================================

"""
$(TYPEDEF)

Bricaud et al. (1995, 1998) power-law phytoplankton absorption. Requires
the `A(λ)` and `E(λ)` spectral tables from JGR 100, 13321 Table 2.
"""
struct Bricaud1995 <: AbstractPhytoplanktonModel end

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
phyto = Phytoplankton{Float64, Bricaud1995}(Chl = 0.5)

# Custom backscatter fraction (coastal, larger cells)
phyto = Phytoplankton{Float64, Bricaud1995}(Chl = 2.0, B = 0.015)
```
"""
struct Phytoplankton{FT, M<:AbstractPhytoplanktonModel,
                        P<:AbstractOceanPhaseFunction{FT}} <: AbstractAbsorbingScatterer{FT}
    "Parameterization-model tag (e.g. `Bricaud1995()`)"
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
    Phytoplankton{FT, Bricaud1995}(; kwargs...)
Phytoplankton(; kwargs...) = Phytoplankton{Float64}(; kwargs...)

phase_function(p::Phytoplankton) = p.phase

# -----------------------------------------------------------------------------
# Absorption — Bricaud (1995)
# -----------------------------------------------------------------------------

const _BRICAUD_1995_A_PATH      = data_path("bricaud_1995_A.csv")
const _BRICAUD_1995_E_PATH      = data_path("bricaud_1995_E.csv")
const _BRICAUD_1995_AVAILABLE   = isfile(_BRICAUD_1995_A_PATH) && isfile(_BRICAUD_1995_E_PATH)
const _BRICAUD_1995_A_CACHE     = Ref{Union{Nothing, SpectralTable{Float64}}}(nothing)
const _BRICAUD_1995_E_CACHE     = Ref{Union{Nothing, SpectralTable{Float64}}}(nothing)

function _bricaud_1995_tables()
    A = _BRICAUD_1995_A_CACHE[]
    E = _BRICAUD_1995_E_CACHE[]
    (A !== nothing && E !== nothing) && return (A = A, E = E)
    isfile(_BRICAUD_1995_A_PATH) && isfile(_BRICAUD_1995_E_PATH) || error(
        "Phytoplankton{Bricaud1995} absorption requires both ",
        "`data/bricaud_1995_A.csv` and `data/bricaud_1995_E.csv`, which ",
        "are not bundled with this OceanOptics.jl release. Transcribe ",
        "Table 2 of Bricaud et al. (1995) JGR 100, 13321 into those ",
        "paths (see DESIGN §8.3); each file is (lambda_nm, value).")
    A = load_spectral_csv(_BRICAUD_1995_A_PATH)
    E = load_spectral_csv(_BRICAUD_1995_E_PATH)
    _BRICAUD_1995_A_CACHE[] = A
    _BRICAUD_1995_E_CACHE[] = E
    return (A = A, E = E)
end

function absorption(p::Phytoplankton{FT, Bricaud1995}, λ::Real) where {FT}
    tbls = _bricaud_1995_tables()
    # a_φ*(λ) = A(λ) · Chl^{-E(λ)}, a_φ(λ) = a_φ*(λ) · Chl
    A_λ  = linterp(tbls.A, λ)
    E_λ  = linterp(tbls.E, λ)
    return A_λ * p.Chl^(one(FT) - E_λ)    # = A·Chl^(-E)·Chl = A·Chl^(1-E)
end

function absorption(p::Phytoplankton{FT, Gordon1992}, λ::Real) where {FT}
    error("Phytoplankton{Gordon1992} absorption is not yet implemented. ",
          "Use Phytoplankton{$(nameof(FT)), Bricaud1995} with the Bricaud ",
          "tables, or wait for a future release.")
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
