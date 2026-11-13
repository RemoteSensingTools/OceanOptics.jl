# =============================================================================
# Pigment partitions of phytoplankton absorption
# =============================================================================
#
# Chlorophyll fluorescence is excited by the light that chlorophyll absorbs
# and by the share of light absorbed by accessory pigments that they pass on
# to chlorophyll (antenna); the rest is dissipated (photoprotection,
# quenching). To account for this, a pigment partition attributes a
# phytoplankton absorption coefficient a_ph(λ) to pigment groups,
#
#     a_k(λ) = f_k(λ) · a_ph(λ),      f_k(λ) = c_k a*_k(λ) / Σ_j c_j a*_j(λ),
#
# from in vivo specific absorption spectra a*_k of the groups and their
# relative contents c_k. The total a_ph — packaging included — still comes
# from the phytoplankton model; only its attribution is added. With
# `Bricaud1998` that total is the particulate absorption a_p, detritus
# included, and the shares apply to all of it. The shares
# come from unpackaged spectra, so packaging is assumed to flatten all
# groups alike at each wavelength.
#
# Outside the tabulated range of the pigment spectra the shares are those of
# the nearest tabulated wavelength (`HoldEdgeShares`, the default), so that no
# absorption is dropped there; `UnassignedOutsideTable` instead puts all of
# it in the part `:unassigned`. Where every group spectrum is zero the
# absorption is `:unassigned` in both cases. Fluorescence counts
# `:unassigned` with η = 0 unless told otherwise.
#
# Two partitions are bundled:
#   ProspectPigments    — land analogy: chlorophyll a+b and carotenoids with
#                         the PROSPECT-D in vivo spectra used by Fluspect.
#   Bidigare1990Pigments — ocean pigment groups (chlorophylls a, b, c,
#                         photosynthetic and photoprotective carotenoids,
#                         phycoerythrins) after Bidigare et al. (1990).
# A partition only needs `pigment_parts`, `pigment_class` and
# `pigment_fractions`.
# =============================================================================

"""
$(TYPEDEF)

Attribution of phytoplankton absorption to pigment groups. Subtypes
implement [`pigment_parts`](@ref), [`pigment_class`](@ref) and
[`pigment_fractions`](@ref).
"""
abstract type AbstractPigmentPartition end

"""
    pigment_parts(p::AbstractPigmentPartition) -> Tuple{Vararg{Symbol}}

Return the names of the pigment groups of `p`; the last is `:unassigned`.
"""
function pigment_parts end

"""
    pigment_class(p::AbstractPigmentPartition, part::Symbol) -> Symbol

Return the role of `part`: `:chlorophyll`, `:carotenoid` (photosynthetic or lumped
carotenoids), `:photoprotective_carotenoid`, `:phycobilin` or `:unassigned`.

The role sets the default energy-transfer efficiency of the part in
[`chlorophyll_fluorescence`](@ref).
"""
function pigment_class end

"""
    pigment_fractions(p::AbstractPigmentPartition, λ) -> Tuple

Return the share of phytoplankton absorption attributed to each part of
[`pigment_parts`](@ref)`(p)` at `λ` `[nm]`; the shares sum to one.

Outside the tabulated pigment spectra (both bundled partitions start at 400 nm) the
result depends on the partition's `outside` rule: with [`HoldEdgeShares`](@ref) (the
default) the shares of the nearest tabulated wavelength apply, so absorption in the UV
is attributed as at 400 nm; with [`UnassignedOutsideTable`](@ref) it is all
`:unassigned`, and [`chlorophyll_fluorescence`](@ref) then has no excitation there
(`:unassigned` has η = 0). Where all group spectra vanish the absorption is
`:unassigned` under both rules.
"""
function pigment_fractions end

# Shares c_k a*_k / Σ c_j a*_j, or everything unassigned where the sum is zero.
function _pigment_shares(s::Tuple{T, Vararg{T}}) where {T}
    S = sum(s)
    S > 0 || return (map(_ -> zero(T), s)..., one(T))
    return (map(x -> x / S, s)..., zero(T))
end

_unassigned_shares(::Val{N}, ::Type{T}) where {N, T} = (ntuple(_ -> zero(T), Val(N))..., one(T))

"""
    AbstractOutOfRangeShares

Rule for the pigment shares at wavelengths outside the tabulated pigment spectra:
[`HoldEdgeShares`](@ref) or [`UnassignedOutsideTable`](@ref).
"""
abstract type AbstractOutOfRangeShares end

"""
    HoldEdgeShares()

Outside the tabulated range the pigment shares are those of the nearest tabulated
wavelength (the default of the bundled partitions).

The partition then attributes all absorption at every wavelength; the shares below
400 nm are an extrapolation, not data.
"""
struct HoldEdgeShares <: AbstractOutOfRangeShares end

"""
    UnassignedOutsideTable()

Outside the tabulated range all absorption goes to the part `:unassigned`, which
[`chlorophyll_fluorescence`](@ref) counts with zero transfer efficiency by default.
"""
struct UnassignedOutsideTable <: AbstractOutOfRangeShares end

# Wavelength at which the share tables are read, or `nothing` for "all unassigned".
_share_wavelength(::HoldEdgeShares, λ, lo, hi) = clamp(λ, lo, hi)
_share_wavelength(::UnassignedOutsideTable, λ, lo, hi) = _within(λ, lo, hi) ? λ : nothing

# -----------------------------------------------------------------------------
# Land analogy: PROSPECT-D chlorophyll and carotenoid spectra
# -----------------------------------------------------------------------------

const _PROSPECT_D_TABLE = _read_csv_columns(data_path("prospect_d_pigment_absorption.csv"))
const _PROSPECT_KAB     = SpectralTable(_PROSPECT_D_TABLE["lambda_nm"], _PROSPECT_D_TABLE["Kab"])
const _PROSPECT_KCA     = SpectralTable(_PROSPECT_D_TABLE["lambda_nm"], _PROSPECT_D_TABLE["Kca"])

"""
$(TYPEDEF)

Land-analogous partition into chlorophyll (a + b) and carotenoids, with
the in vivo specific absorption spectra `Kab` and `Kca` of the leaf model
PROSPECT-D (Féret et al. 2017) that Fluspect uses
(`data/prospect_d_pigment_absorption.csv`, 400–800 nm). Parts:
`(:chlorophyll, :carotenoids, :unassigned)`. Outside 400–800 nm the shares at
400 nm or 800 nm apply (`outside = HoldEdgeShares()`, the default), or all
absorption is `:unassigned` (`outside = UnassignedOutsideTable()`).

# Fields
$(TYPEDFIELDS)

# Example
```julia
phyto = Phytoplankton(Chl = 0.5, pigments = ProspectPigments(carotenoid_ratio = 0.25))
absorption_by_part(phyto, 450.0)    # (a_chl, a_car, 0.0)
```
SCOPE's default leaf has 40 µg cm⁻² chlorophyll and 10 µg cm⁻²
carotenoids, a ratio of 0.25; phytoplankton ratios differ between groups
and with light history.
"""
struct ProspectPigments{FT<:Real, X<:AbstractOutOfRangeShares} <: AbstractPigmentPartition
    "Carotenoid-to-chlorophyll mass ratio"
    carotenoid_ratio::FT
    "Rule for the shares outside the tabulated 400–800 nm"
    outside::X

    function ProspectPigments{FT, X}(carotenoid_ratio, outside::X) where {FT<:Real, X}
        isfinite(carotenoid_ratio) && carotenoid_ratio ≥ 0 ||
            throw(ArgumentError("ProspectPigments: carotenoid_ratio must be finite and " *
                                "non-negative, got $carotenoid_ratio"))
        new{FT, X}(carotenoid_ratio, outside)
    end
end

ProspectPigments(; carotenoid_ratio::Real,
                 outside::AbstractOutOfRangeShares = HoldEdgeShares()) =
    ProspectPigments{float(typeof(carotenoid_ratio)), typeof(outside)}(carotenoid_ratio, outside)

pigment_parts(::ProspectPigments) = (:chlorophyll, :carotenoids, :unassigned)

pigment_class(p::ProspectPigments, part::Symbol) =
    part === :chlorophyll ? :chlorophyll :
    part === :carotenoids ? :carotenoid  :
    part === :unassigned  ? :unassigned  : _unknown_pigment(p, part)

function pigment_fractions(p::ProspectPigments, λ::Real)
    T = promote_type(typeof(p.carotenoid_ratio), typeof(float(λ)))
    x = _share_wavelength(p.outside, λ, first(_PROSPECT_KAB.λ), last(_PROSPECT_KAB.λ))
    x === nothing && return _unassigned_shares(Val(2), T)
    return _pigment_shares((T(linterp(_PROSPECT_KAB, x)),
                            T(p.carotenoid_ratio * linterp(_PROSPECT_KCA, x))))
end

function Base.show(io::IO, p::ProspectPigments)
    print(io, "ProspectPigments(carotenoid_ratio = ", p.carotenoid_ratio, ")")
end

# -----------------------------------------------------------------------------
# Ocean pigment groups: Bidigare et al. (1990)
# -----------------------------------------------------------------------------

const _BIDIGARE_TABLE   = _read_csv_columns(data_path("bidigare_1990_pigments.csv"))
const _BIDIGARE_COLUMNS = ("chl_a", "chl_b", "chl_c", "psc", "ppc", "pe_peb", "pe_pub")
const _BIDIGARE_SPECTRA = Tuple(SpectralTable(_BIDIGARE_TABLE["lambda_nm"], _BIDIGARE_TABLE[c])
                                for c in _BIDIGARE_COLUMNS)

"""
$(TYPEDEF)

Partition into the algal pigment groups of Bidigare et al. (1990), with
their in vivo weight-specific absorption spectra
(`data/bidigare_1990_pigments.csv`, 400–750 nm, from the Ocean Optics Web
Book). Fields are pigment concentrations `[mg m⁻³]`, e.g. from HPLC; only
their ratios matter, since the total comes from the phytoplankton model.

Parts: `(:chlorophyll_a, :chlorophyll_b, :chlorophyll_c,
:photosynthetic_carotenoids, :photoprotective_carotenoids,
:phycoerythrin_peb, :phycoerythrin_pub, :unassigned)`. Below 400 nm the
shares at 400 nm apply (`outside = HoldEdgeShares()`, the default) or all
absorption is `:unassigned` (`outside = UnassignedOutsideTable()`). Where all
group spectra vanish (above the chlorophyll a red band, ≈ 732 nm, and so
also above 750 nm) all absorption is `:unassigned` under both rules.

# Fields
$(TYPEDFIELDS)

# Example
```julia
pig   = Bidigare1990Pigments(chl_a = 0.5, chl_c = 0.1, psc = 0.3, ppc = 0.15)
phyto = Phytoplankton(Chl = 0.5, pigments = pig)
absorption_part(phyto, :photoprotective_carotenoids, 470.0)
```
"""
struct Bidigare1990Pigments{FT<:Real, X<:AbstractOutOfRangeShares} <: AbstractPigmentPartition
    "Chlorophyll a `[mg m⁻³]`"
    chl_a::FT
    "Chlorophyll b `[mg m⁻³]`"
    chl_b::FT
    "Chlorophyll c `[mg m⁻³]`"
    chl_c::FT
    "Photosynthetic carotenoids (e.g. fucoxanthin, peridinin) `[mg m⁻³]`"
    psc::FT
    "Photoprotective carotenoids (e.g. zeaxanthin, diadinoxanthin, β-carotene) `[mg m⁻³]`"
    ppc::FT
    "Phycoerythrin rich in phycoerythrobilin (Synechococcus WH7803 type) `[mg m⁻³]`"
    pe_peb::FT
    "Phycoerythrin rich in phycourobilin (Synechococcus WH8103 type) `[mg m⁻³]`"
    pe_pub::FT
    "Rule for the shares outside the tabulated 400–750 nm"
    outside::X

    function Bidigare1990Pigments{FT, X}(outside::X, c...) where {FT<:Real, X}
        all(x -> isfinite(x) && x ≥ 0, c) ||
            throw(ArgumentError("Bidigare1990Pigments: concentrations must be finite and " *
                                "non-negative, got $c"))
        any(>(0), c) ||
            throw(ArgumentError("Bidigare1990Pigments: at least one concentration must be positive"))
        new{FT, X}(c..., outside)
    end
end

function Bidigare1990Pigments(; chl_a::Real, chl_b::Real = 0, chl_c::Real = 0, psc::Real = 0,
                              ppc::Real = 0, pe_peb::Real = 0, pe_pub::Real = 0,
                              outside::AbstractOutOfRangeShares = HoldEdgeShares())
    c  = (chl_a, chl_b, chl_c, psc, ppc, pe_peb, pe_pub)
    FT = float(promote_type(map(typeof, c)...))
    return Bidigare1990Pigments{FT, typeof(outside)}(outside, map(FT, c)...)
end

pigment_parts(::Bidigare1990Pigments) =
    (:chlorophyll_a, :chlorophyll_b, :chlorophyll_c, :photosynthetic_carotenoids,
     :photoprotective_carotenoids, :phycoerythrin_peb, :phycoerythrin_pub, :unassigned)

function pigment_class(p::Bidigare1990Pigments, part::Symbol)
    part in (:chlorophyll_a, :chlorophyll_b, :chlorophyll_c) && return :chlorophyll
    part === :photosynthetic_carotenoids && return :carotenoid
    part === :photoprotective_carotenoids && return :photoprotective_carotenoid
    part in (:phycoerythrin_peb, :phycoerythrin_pub) && return :phycobilin
    part === :unassigned && return :unassigned
    return _unknown_pigment(p, part)
end

_concentrations(p::Bidigare1990Pigments) = (p.chl_a, p.chl_b, p.chl_c, p.psc, p.ppc, p.pe_peb, p.pe_pub)

function pigment_fractions(p::Bidigare1990Pigments{FT}, λ::Real) where {FT}
    T = promote_type(FT, typeof(float(λ)))
    tbl = first(_BIDIGARE_SPECTRA)
    x = _share_wavelength(p.outside, λ, first(tbl.λ), last(tbl.λ))
    x === nothing && return _unassigned_shares(Val(7), T)
    c = _concentrations(p)
    return _pigment_shares(ntuple(i -> T(c[i] * linterp(_BIDIGARE_SPECTRA[i], x)), Val(7)))
end

function Base.show(io::IO, p::Bidigare1990Pigments)
    names = ("chl_a", "chl_b", "chl_c", "psc", "ppc", "pe_peb", "pe_pub")
    set   = [string(n, " = ", c) for (n, c) in zip(names, _concentrations(p)) if c > 0]
    print(io, "Bidigare1990Pigments(", join(set, ", "), ")")
end

_unknown_pigment(p, part) =
    throw(ArgumentError("$(nameof(typeof(p))) has no pigment part :$part; its parts are " *
                        "$(pigment_parts(p))"))
