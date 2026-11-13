# =============================================================================
# Phytoplankton-bearing Case-1 water: particulate absorption and scattering
# =============================================================================
#
# Absorption is the chlorophyll power law of Bricaud et al. (1998, JGR 103,
# 31033) for the particulate absorption of Case-1 water, phytoplankton plus
# co-varying detritus:
#
#     a_p(λ) = A(λ) · Chl^E(λ)        [1/m],  Chl in mg/m³
#
# with A(440) = 0.0520 and E(440) = 0.6350. `data/bricaud_1998_{A,E}.csv` hold
# A(λ) and E(λ) as tabulated for HydroLight's Case-1 model (`AEmidUVabs.txt`),
# which extends the 1998 fit to the UV after Morrison & Nelson (2004) and
# Vasilkov et al. (2005).
#
# Scattering is the total particulate scattering of Case-1 water of Gordon &
# Morel (1983), as cited by Fell (1997), Eq. 5.5:
#
#     b_p(λ) = 0.30 · Chl^0.62 · (550/λ)          [1/m]
#
# The angular distribution is Fournier-Forand with a typical open-ocean
# backscatter fraction B ≈ 0.01; the caller can override at construction.
# =============================================================================

"""
$(TYPEDEF)

Particulate absorption of Case-1 water after Bricaud et al. (1998, JGR 103, 31033).

The absorption is ``a_p(λ) = A(λ) Chl^{E(λ)}`` in `1/m` with `Chl` in `mg/m³`
(`A(440) = 0.0520`, `E(440) = 0.6350`). It is the absorption of all particles,
phytoplankton plus the detritus that co-varies with it, not of phytoplankton alone:
adding a separate [`NonAlgalParticles`](@ref) constituent counts detrital absorption
twice, unless a detrital or mineral load beyond the Case-1 background is intended.

`A(λ)` and `E(λ)` are the tables of HydroLight's Case-1 model (`AEmidUVabs.txt`),
bundled as `data/bricaud_1998_A.csv` and `data/bricaud_1998_E.csv`: 300–1000 nm at
2 nm, with the UV part after Morrison & Nelson (2004) and Vasilkov et al. (2005).
`Chl = 0` gives zero absorption.
"""
struct Bricaud1998 <: AbstractPhytoplanktonModel end

"""
$(TYPEDEF)

Scattering law of the particles of a [`Phytoplankton`](@ref) constituent.
"""
abstract type AbstractParticleScatteringModel end

"""
$(TYPEDEF)

Total particulate scattering of Case-1 water after Gordon & Morel (1983), as cited by
Fell (1997), Eq. 5.5:

``b_p(λ) = b_\\mathrm{ref}\\, Chl^{m}\\, (λ_\\mathrm{ref}/λ)^{γ}`` `[1/m]`, `Chl` in `mg/m³`,

with `b_ref = 0.30` m⁻¹ at `λ_ref = 550` nm, `m = 0.62` and `γ = 1`. Fell gives
`b_ref = 0.30 ± 0.15`. It is the scattering of all particles that co-vary with
chlorophyll, not of living cells alone; in coastal or turbid water the particulate
scattering can be larger by a factor of 2 or more.

# Fields
$(TYPEDFIELDS)
"""
struct GordonMorel1983{FT<:Real} <: AbstractParticleScatteringModel
    "Scattering coefficient at `λ_ref` for `Chl = 1 mg/m³` `[1/m]`"
    b_ref::FT
    "Exponent `m` of the chlorophyll concentration"
    chl_exponent::FT
    "Reference wavelength `[nm]`"
    λ_ref::FT
    "Exponent `γ` of the spectral law `(λ_ref/λ)^γ`"
    spectral_exponent::FT

    function GordonMorel1983{FT}(b_ref, chl_exponent, λ_ref, spectral_exponent) where {FT<:Real}
        _check_parameter_type(FT)
        b_ref ≥ 0 || throw(ArgumentError("GordonMorel1983: b_ref must be non-negative, got $b_ref"))
        λ_ref > 0 || throw(ArgumentError("GordonMorel1983: λ_ref must be positive, got $λ_ref"))
        new{FT}(b_ref, chl_exponent, λ_ref, spectral_exponent)
    end
end

GordonMorel1983{FT}(; b_ref = 0.30, chl_exponent = 0.62, λ_ref = 550,
                    spectral_exponent = 1) where {FT<:Real} =
    GordonMorel1983{FT}(b_ref, chl_exponent, λ_ref, spectral_exponent)
GordonMorel1983(; kwargs...) = GordonMorel1983{_parameter_type(values(kwargs))}(; kwargs...)

# A Phytoplankton constructed without an explicit element type promotes with the element
# type of its scattering model, as with its real parameters (a dual model gives a dual
# constituent).
_promote_reals(::Type{T}, m::GordonMorel1983{S}, rest...) where {T, S} =
    _promote_reals(promote_type(T, S), rest...)

# A floating-point scattering model is evaluated in double precision and converted to
# the constituent's type; a model of any other element type (a dual) must be reachable
# from the constituent's type by promotion, or its derivative could not be stored.
_check_model_precision(::Type, ::AbstractParticleScatteringModel) = nothing
function _check_model_precision(::Type{FT}, m::GordonMorel1983{S}) where {FT, S}
    S <: AbstractFloat || promote_type(FT, S) === FT ||
        throw(ArgumentError("Phytoplankton{$FT}: the scattering model has element type $S, " *
                            "which $FT cannot hold; construct Phytoplankton{$(promote_type(FT, S))} " *
                            "or omit the element type so that it follows the model"))
    return nothing
end

"""
$(TYPEDEF)

Particles of Case-1 water (phytoplankton and co-varying detritus) as an absorbing
scatterer, parameterized by the chlorophyll-a concentration.

`Chl` is the chlorophyll-a concentration `[mg/m³]`. The model tag `M` selects the
absorption parameterization; with [`Bricaud1998`](@ref) it is the particulate
absorption ``a_p = A(λ) Chl^{E(λ)}``, which includes detritus (see there before
combining with [`NonAlgalParticles`](@ref)). Scattering is the total particulate
scattering of Case-1 water given by `scattering_model`, by default
[`GordonMorel1983`](@ref), ``b_p = 0.30 Chl^{0.62} (550/λ)``, and the phase function
defaults to `FournierForandPhase` with backscatter fraction `B = 0.01`.

An optional pigment partition ([`ProspectPigments`](@ref),
[`Bidigare1990Pigments`](@ref)) attributes the absorption to pigment
groups: [`absorption_parts`](@ref) then names the groups,
[`layer_optics`](@ref) reports one absorption component per group, and
[`chlorophyll_fluorescence`](@ref) can weight each group by its energy
transfer to chlorophyll. The total absorption is unchanged.

# Fields
$(TYPEDFIELDS)

# Examples

```julia
# Open-ocean Case-1 water at Chl = 0.5 mg/m³
phyto = Phytoplankton{Float64, Bricaud1998}(Chl = 0.5)

# Custom backscatter fraction (coastal, larger cells)
phyto = Phytoplankton{Float64, Bricaud1998}(Chl = 2.0, B = 0.015)

# Absorption attributed to chlorophyll and carotenoids
phyto = Phytoplankton(Chl = 0.5, pigments = ProspectPigments(carotenoid_ratio = 0.3))
absorption_parts(phyto)           # (:chlorophyll, :carotenoids, :unassigned)

phyto = Phytoplankton(Chl = 0.5, scattering_model = GordonMorel1983(b_ref = 0.45))
```
"""
struct Phytoplankton{FT, M<:AbstractPhytoplanktonModel, P<:AbstractOceanPhaseFunction{FT},
                     G<:Union{Nothing, AbstractPigmentPartition},
                     S<:AbstractParticleScatteringModel} <: AbstractAbsorbingScatterer{FT}
    "Parameterization-model tag (e.g. `Bricaud1998()`)"
    model::M
    "Chlorophyll-a concentration `[mg/m³]`"
    Chl::FT
    "Particulate-scattering phase function"
    phase::P
    "Pigment partition of the absorption, or `nothing` (a single part, `:total`)"
    pigments::G
    "Scattering law and its parameters (default [`GordonMorel1983`](@ref))"
    scattering_model::S
end

"""
    Phytoplankton{FT, M}(; Chl, B = 0.01, pigments = nothing,
                           scattering_model = GordonMorel1983())

Construct a [`Phytoplankton`](@ref) with chlorophyll-a concentration `Chl` `[mg/m³]`, a
Fournier-Forand phase function of backscatter fraction `B`, optionally a pigment
partition of its absorption, and the scattering law `scattering_model`.

The parameters of `scattering_model` are kept in the precision they are given in
(double by default); the scattering coefficient is evaluated in double precision and
returned in the precision of the constituent.
"""
function Phytoplankton{FT, M}(; Chl, B = 0.01,
                              pigments::Union{Nothing, AbstractPigmentPartition} = nothing,
                              scattering_model::AbstractParticleScatteringModel = GordonMorel1983()
                              ) where {FT<:Real, M<:AbstractPhytoplanktonModel}
    _check_parameter_type(FT)
    _check_model_precision(FT, scattering_model)
    Chl ≥ zero(Chl) ||
        throw(ArgumentError("Phytoplankton: Chl must be non-negative, got $Chl"))
    phase = FournierForandPhase{FT}(backscatter_fraction = FT(B))
    pg    = _partition_precision(FT, pigments)
    return Phytoplankton{FT, M, typeof(phase), typeof(pg), typeof(scattering_model)}(
        M(), FT(Chl), phase, pg, scattering_model)
end

# Floating-point partitions take the phytoplankton's precision, so that the
# absorption parts have the type of the total; dual numbers (AD) pass through.
_partition_precision(::Type{FT}, p) where {FT} = p
_partition_precision(::Type{FT}, p::ProspectPigments{<:AbstractFloat, X}) where {FT, X} =
    ProspectPigments{FT, X}(p.carotenoid_ratio, p.outside)
_partition_precision(::Type{FT}, p::Bidigare1990Pigments{<:AbstractFloat, X}) where {FT, X} =
    Bidigare1990Pigments{FT, X}(p.outside, _concentrations(p)...)

Phytoplankton{FT}(; kwargs...) where {FT<:Real} =
    Phytoplankton{FT, Bricaud1998}(; kwargs...)
Phytoplankton(; kwargs...) = Phytoplankton{_parameter_type(values(kwargs))}(; kwargs...)

scattering_phase(p::Phytoplankton) = p.phase

# -----------------------------------------------------------------------------
# Absorption — Bricaud et al. (1998) particulate absorption, HydroLight tables
# -----------------------------------------------------------------------------

"Coefficient `A(λ)` `[1/m]` of `a_p = A·Chl^E`, 300–1000 nm (HydroLight `AEmidUVabs.txt`)."
const _BRICAUD_1998_A_TABLE = load_spectral_csv(data_path("bricaud_1998_A.csv"))

"Exponent `E(λ)` of `a_p = A·Chl^E`, 300–1000 nm (HydroLight `AEmidUVabs.txt`)."
const _BRICAUD_1998_E_TABLE = load_spectral_csv(data_path("bricaud_1998_E.csv"))

function absorption(p::Phytoplankton{FT, Bricaud1998}, λ::Real) where {FT}
    # a_p(λ) = A(λ) · Chl^E(λ); E > 0 over the whole table, so Chl = 0 gives 0.
    x   = _wide(λ)
    A_λ = linterp(_BRICAUD_1998_A_TABLE, x)
    E_λ = linterp(_BRICAUD_1998_E_TABLE, x)
    return convert(_result_type(FT, λ), A_λ * _wide(p.Chl)^E_λ)
end

absorption(p::Phytoplankton, λ::AbstractVector) = map(Base.Fix1(absorption, p), λ)

absorption_parts(p::Phytoplankton) = p.pigments === nothing ? (:total,) : pigment_parts(p.pigments)

function absorption_by_part(p::Phytoplankton, λ::Real)
    a = absorption(p, λ)
    pig = p.pigments
    return pig === nothing ? (a,) : map(f -> f * a, pigment_fractions(pig, λ))
end

# -----------------------------------------------------------------------------
# Scattering — the law of the constituent's `scattering_model`
# -----------------------------------------------------------------------------

scattering(p::Phytoplankton{FT}, λ::Real) where {FT} =
    convert(_result_type(FT, λ), _particle_scattering(p.scattering_model, _wide(p.Chl), _wide(λ)))

# b_p(λ) = b_ref · Chl^m · (λ_ref/λ)^γ
_particle_scattering(m::GordonMorel1983, Chl, λ) =
    m.b_ref * Chl^m.chl_exponent * (m.λ_ref / λ)^m.spectral_exponent

scattering(p::Phytoplankton, λ::AbstractVector) = map(Base.Fix1(scattering, p), λ)

# Backscatter — phase-function-determined. Delegates to the configured
# Fournier-Forand's closed-form `backscatter_fraction`.
function backscattering(p::Phytoplankton{FT}, λ::Real) where {FT}
    return scattering(p, λ) * backscatter_fraction(p.phase)
end

backscattering(p::Phytoplankton, λ::AbstractVector) = map(Base.Fix1(backscattering, p), λ)

function Base.show(io::IO, p::Phytoplankton{FT, M}) where {FT, M}
    print(io, "Phytoplankton{$FT}(", nameof(M), ", Chl=$(p.Chl) mg/m³, ",
              "B=", backscatter_fraction(p.phase),
              p.pigments === nothing ? "" : ", $(p.pigments)", ")")
end
