# =============================================================================
# Fluorescence processes
# =============================================================================
#
# A fluorescence process turns absorbed photons into emitted photons. With
# quantum yield φ and excitation pathways k,
#
#     a^I(λ') = φ · Σ_k η_k · a_k(λ')        [1/m]
#     f(λ', λ) = f(λ)                         [1/nm], ∫ f dλ = 1
#
# where a_k is the absorption of one absorber, or one pigment group of it
# (see `absorption_parts`), and η_k ∈ [0, 1] is the fraction of the photons
# it absorbs that reach the fluorescing chlorophyll. Chlorophyll has η = 1;
# accessory pigments acting as antenna pass on a fraction η ≤ 1; light
# absorbed for photoprotection, or otherwise quenched, has η = 0.
#
#   IsotropicFluorescence — generic isotropic fluorescence from any set of
#                           excitation pathways and emission spectrum.
#   chlorophyll_fluorescence — chlorophyll fluorescence of `Phytoplankton`,
#                           with the chloroplast-level emission of the land
#                           leaf model Fluspect by default and pathways
#                           from the phytoplankton's pigment partition.
#   CDOMFluorescence      — CDOM fluorescence, a separate type so solvers
#                           can dispatch on it.
#
# Quantum yields have no defensible universal value and are always given
# explicitly.
# =============================================================================

"""
$(TYPEDEF)

One excitation pathway of a fluorescence process: the photons absorbed by
part `part` of `absorber` (see [`absorption_parts`](@ref)), of which the
fraction `efficiency` is passed on to the fluorescing pigment.

# Fields
$(TYPEDFIELDS)

# Construction

    ExcitationPathway(absorber, part = :total, efficiency = 1)

`efficiency` (η) lies in `[0, 1]`: 1 for chlorophyll itself, the energy
transfer efficiency for an antenna pigment, 0 for absorption that is
dissipated (photoprotective pigments, quenching).
"""
struct ExcitationPathway{FT<:AbstractFloat, A<:AbstractAbsorber{FT}}
    "Absorbing constituent; must also be a constituent of every layer the process is in"
    absorber::A
    "Absorption part of `absorber`"
    part::Symbol
    "Fraction of absorbed photons passed to the fluorescing pigment (η)"
    efficiency::FT

    function ExcitationPathway{FT, A}(absorber::A, part::Symbol,
                                      efficiency) where {FT<:AbstractFloat, A<:AbstractAbsorber{FT}}
        part in absorption_parts(absorber) ||
            throw(ArgumentError("ExcitationPathway: $(nameof(A)) has no absorption part :$part; " *
                                "its parts are $(absorption_parts(absorber))"))
        isfinite(efficiency) && 0 ≤ efficiency ≤ 1 ||
            throw(ArgumentError("ExcitationPathway: efficiency must lie in [0, 1], got $efficiency"))
        new{FT, A}(absorber, part, efficiency)
    end
end

ExcitationPathway(absorber::AbstractAbsorber{FT}, part::Symbol = :total,
                  efficiency::Real = 1) where {FT} =
    ExcitationPathway{FT, typeof(absorber)}(absorber, part, FT(efficiency))

function Base.show(io::IO, x::ExcitationPathway)
    print(io, "ExcitationPathway(", nameof(typeof(x.absorber)), ", :", x.part, ", η = ",
          x.efficiency, ")")
end

"""
    excitation_pathways(p::AbstractOceanInelasticProcess) -> Tuple

Return the [`ExcitationPathway`](@ref)s whose absorbed photons drive `p`; the tuple is
empty for processes not driven by absorption (water Raman).
"""
excitation_pathways(::AbstractOceanInelasticProcess) = ()

_pathways(x::ExcitationPathway)     = (x,)
_pathways(x::AbstractAbsorber)      = (ExcitationPathway(x),)
_pathways(x::Union{Tuple, AbstractVector}) = Tuple(Iterators.flatten(map(_pathways, x)))

"""
$(TYPEDEF)

Isotropic fluorescence: photons absorbed through the excitation pathways
are re-emitted with quantum yield `quantum_yield` and spectral
distribution `emission`.

# Fields
$(TYPEDFIELDS)

# Construction

    IsotropicFluorescence(excitation; emission, quantum_yield, excitation_range)

`excitation` is an absorber (all of its absorption, η = 1), an
[`ExcitationPathway`](@ref), or a collection of them. The excitation
coefficient is

`excitation_absorption(p, λ′) = φ · Σ_k η_k · absorption_part(absorber_k, part_k, λ′)`.

# Example
```julia
phyto = Phytoplankton(Chl = 0.5, pigments = ProspectPigments(carotenoid_ratio = 0.3))
f = IsotropicFluorescence([ExcitationPathway(phyto, :chlorophyll, 1.0),
                           ExcitationPathway(phyto, :carotenoids, 0.5)];
                          emission = chloroplast_emission(), quantum_yield = 0.02,
                          excitation_range = (400.0, 700.0))
```
See [`chlorophyll_fluorescence`](@ref) for the same with defaults per
pigment group.
"""
struct IsotropicFluorescence{FT, X<:Tuple{Vararg{ExcitationPathway{FT}}},
                             E<:AbstractSpectralRedistribution{FT}} <: AbstractOceanInelasticProcess{FT}
    "Excitation pathways"
    excitation::X
    "Spectral redistribution `f(λ', λ)` `[1/nm]` of the emitted photons"
    emission::E
    "Fluorescence quantum yield φ: photons emitted per photon delivered to the fluorescing pigment"
    quantum_yield::FT
    "Excitation-wavelength integration window `[nm]`"
    excitation_range::Tuple{FT, FT}

    function IsotropicFluorescence{FT, X, E}(excitation::X, emission::E, quantum_yield,
                                             excitation_range) where {FT, X, E}
        isempty(excitation) &&
            throw(ArgumentError("IsotropicFluorescence: no excitation pathway"))
        _check_yield_and_range(quantum_yield, excitation_range, "IsotropicFluorescence")
        new{FT, X, E}(excitation, emission, quantum_yield, excitation_range)
    end
end

function IsotropicFluorescence(excitation; emission::AbstractSpectralRedistribution{FT},
                               quantum_yield::Real,
                               excitation_range::Tuple{<:Real, <:Real}) where {FT}
    x = _pathways(excitation)
    all(p -> p isa ExcitationPathway{FT}, x) ||
        throw(ArgumentError("IsotropicFluorescence: absorbers and emission must share the " *
                            "floating-point type $FT"))
    return IsotropicFluorescence{FT, typeof(x), typeof(emission)}(
        x, emission, FT(quantum_yield), (FT(excitation_range[1]), FT(excitation_range[2])))
end

function _check_yield_and_range(φ, range, who)
    isfinite(φ) && 0 ≤ φ ≤ 1 ||
        throw(ArgumentError("$who: quantum_yield must be in [0, 1], got $φ"))
    0 < range[1] < range[2] ||
        throw(ArgumentError("$who: excitation_range must be (lo, hi) with 0 < lo < hi, got $range"))
    return nothing
end

# a^I(λ') = φ · Σ_k η_k a_k(λ')
function excitation_absorption(p::IsotropicFluorescence, λ_prime::Real)
    return p.quantum_yield * sum(x.efficiency * absorption_part(x.absorber, x.part, λ_prime)
                                 for x in p.excitation)
end

emission(p::IsotropicFluorescence, λ_prime::Real, λ::Real) =
    redistribution(p.emission, λ_prime, λ)

excitation_range(p::IsotropicFluorescence)    = p.excitation_range
excitation_pathways(p::IsotropicFluorescence) = p.excitation
is_isotropic(::IsotropicFluorescence)         = true
conversion_units(::IsotropicFluorescence)     = PhotonConversion()
"""
    quantum_yield(p)

Return the fluorescence quantum yield of process `p`: photons emitted per photon
delivered to the fluorescing pigment.
"""
quantum_yield(p::IsotropicFluorescence)       = p.quantum_yield

function Base.show(io::IO, p::IsotropicFluorescence{FT}) where {FT}
    print(io, "IsotropicFluorescence{$FT}(φ = ", p.quantum_yield, ", ",
          join(("$(nameof(typeof(x.absorber))):$(x.part) η=$(x.efficiency)" for x in p.excitation), ", "),
          ", ", p.emission, ")")
end

# -----------------------------------------------------------------------------
# Chloroplast-level emission (land analogy, Fluspect)
# -----------------------------------------------------------------------------

const _FLUSPECT_EMISSION = _read_csv_columns(data_path("fluspect_chloroplast_emission.csv"))
const _FLUSPECT_SOURCE   = "Fluspect Optipar2017_ProspectD (SCOPE; Vilfan et al. 2016, 2018)"

"""
    chloroplast_emission([FT = Float64]; psi_share = nothing) -> TabulatedEmission

Return the chlorophyll fluorescence emission spectrum at the photosystem (chloroplast)
level, before any reabsorption.

The spectra are those used for leaves by Fluspect
(`data/fluspect_chloroplast_emission.csv`, 640–851 nm):

- `psi_share = nothing`: Fluspect's combined spectrum `phi` (the default
  of SCOPE 2), with maxima at 685 and 739 nm;
- `psi_share = w ∈ [0, 1]`: the photosystem mixture
  `w · phiI + (1 − w) · phiII` of the separate PSI (maximum at 728 nm) and
  PSII (682 nm) spectra, with `w` the share of emitted photons from PSI.

The spectra are photon densities per nm. Reabsorption by water and by the
pigments themselves is left to the radiative transfer: the emission
shape contains none.
"""
function chloroplast_emission(::Type{FT} = Float64;
                              psi_share::Union{Nothing, Real} = nothing) where {FT<:AbstractFloat}
    λ = _FLUSPECT_EMISSION["lambda_nm"]
    if psi_share === nothing
        # phi carries one numerical artefact, -4.5e-6 at 641 nm (see the data file)
        shape  = max.(_FLUSPECT_EMISSION["phi"], 0)
        source = _FLUSPECT_SOURCE * ", combined phi"
    else
        0 ≤ psi_share ≤ 1 ||
            throw(ArgumentError("chloroplast_emission: psi_share must lie in [0, 1], got $psi_share"))
        shape  = psi_share .* _FLUSPECT_EMISSION["phiI"] .+
                 (1 - psi_share) .* _FLUSPECT_EMISSION["phiII"]
        source = _FLUSPECT_SOURCE * ", PSI share $psi_share"
    end
    return TabulatedEmission(FT.(λ), FT.(shape); source)
end

"Limits of photosynthetically available radiation `[nm]`."
const _PAR_RANGE_NM = (400.0, 700.0)

"""
    chlorophyll_fluorescence(phyto::Phytoplankton; quantum_yield,
                             carotenoid_transfer = nothing, phycobilin_transfer = 1,
                             photoprotective_transfer = 0, transfer = (;), emission = chloroplast_emission(FT),
                             excitation_range = (400, 700)) -> IsotropicFluorescence

Build the chlorophyll-a fluorescence process of `phyto`, with one excitation pathway per
absorption part of `phyto`.

| part | default η |
|:-----|:----------|
| `:total` (no pigment partition) | 1: all phytoplankton absorption |
| chlorophylls | 1 |
| carotenoids (photosynthetic, or Prospect's lumped carotenoids) | `carotenoid_transfer`: required when the partition has them |
| photoprotective carotenoids | `photoprotective_transfer` (0) |
| phycobilins | `phycobilin_transfer` |
| `:unassigned` | 0 |

Phycobilisomes are a photosystem II antenna; the default `phycobilin_transfer = 1` is a
modelling choice. Photoprotective carotenoids dissipate the energy they absorb instead of
passing it on to chlorophyll, hence `photoprotective_transfer = 0`.

`transfer` overrides η for named parts, e.g.
`(photosynthetic_carotenoids = 0.8, photoprotective_carotenoids = 0.0)`
for [`Bidigare1990Pigments`](@ref). `carotenoid_transfer = 0` is
chlorophyll-only excitation (as in Fluspect, whose excitation is the
chlorophyll share of leaf absorption); `1` makes carotenoids a lossless
antenna; values between give the share that is passed on rather than
quenched.

The pigment spectra of the partitions start at 400 nm. With the default
[`HoldEdgeShares`](@ref) rule the shares at 400 nm apply below, so excitation in the UV
(`excitation_range` starting below 400 nm) is attributed to the pigments as at 400 nm
and nothing is dropped; a partition built with `outside = UnassignedOutsideTable()`
has no excitation outside its table. Without a partition, all absorption within
`excitation_range` excites. With
[`Bricaud1998`](@ref) the absorption is the particulate absorption, detritus included,
and all of it enters the excitation.

`quantum_yield` (φ) has no default: in situ values range over roughly
0.01–0.05 (Loisel et al. 2023, who use HydroLight's 0.02). The default emission is
Fluspect's chloroplast-level spectrum ([`chloroplast_emission`](@ref)), the default
excitation window PAR.
`phyto` must also be a constituent of every layer the process is in; see
[`layer_optics`](@ref).

# Example
```julia
phyto = Phytoplankton(Chl = 0.5, pigments = ProspectPigments(carotenoid_ratio = 0.3))
chl_only = chlorophyll_fluorescence(phyto; quantum_yield = 0.02, carotenoid_transfer = 0.0)
antenna  = chlorophyll_fluorescence(phyto; quantum_yield = 0.02, carotenoid_transfer = 0.6)
```
"""
function chlorophyll_fluorescence(phyto::Phytoplankton{FT};
        quantum_yield::Real,
        carotenoid_transfer::Union{Nothing, Real} = nothing,
        phycobilin_transfer::Real = 1,
        photoprotective_transfer::Real = 0,
        transfer::NamedTuple = NamedTuple(),
        emission::AbstractSpectralRedistribution{FT} = chloroplast_emission(FT),
        excitation_range::Tuple{<:Real, <:Real} = _PAR_RANGE_NM,
) where {FT<:AbstractFloat}
    parts = absorption_parts(phyto)
    for k in keys(transfer)
        k in parts || throw(ArgumentError("chlorophyll_fluorescence: `transfer` names :$k, " *
                                          "which is not an absorption part of the phytoplankton $parts"))
    end
    if phyto.pigments === nothing && carotenoid_transfer !== nothing
        throw(ArgumentError("chlorophyll_fluorescence: carotenoid_transfer needs a pigment " *
                            "partition; construct the Phytoplankton with `pigments = " *
                            "ProspectPigments(...)` or `Bidigare1990Pigments(...)`"))
    end
    function efficiency(part)
        haskey(transfer, part) && return transfer[part]
        part === :total && return 1
        class = pigment_class(phyto.pigments, part)
        class === :chlorophyll && return 1
        class === :phycobilin  && return phycobilin_transfer
        class === :photoprotective_carotenoid && return photoprotective_transfer
        class === :unassigned  && return 0
        carotenoid_transfer === nothing &&
            throw(ArgumentError("chlorophyll_fluorescence: choose the carotenoid energy " *
                                "transfer: `carotenoid_transfer = 0` excites by chlorophyll " *
                                "absorption only, `1` uses carotenoids as a lossless antenna, " *
                                "values between pass on that share (or set :$part in `transfer`)"))
        return carotenoid_transfer
    end
    pathways = Tuple(ExcitationPathway(phyto, part, efficiency(part)) for part in parts)
    return IsotropicFluorescence(pathways; emission, quantum_yield, excitation_range)
end

# =============================================================================
# CDOM fluorescence
# =============================================================================

"""
$(TYPEDEF)

CDOM fluorescence as an isotropic inelastic process. Structurally equivalent
to [`IsotropicFluorescence`](@ref) with a single [`CDOM`](@ref) pathway,
but carried as a distinct type so solvers can dispatch separately and
users can read model intent directly from the type.

Typical excitation/emission footprint per Coble (1996, 2007) and
Bricaud-Morel-Prieur (1981): broad CDOM emission around 420–450 nm for
UV excitation (300–380 nm), with quantum yield 1–2 % when integrated
over the band. Fell (1997) uses a simple Gaussian model.

# Fields
$(TYPEDFIELDS)

# Construction

    CDOMFluorescence(cdom; emission, quantum_yield, excitation_range = (300, 400))
"""
struct CDOMFluorescence{FT, C<:CDOM{FT}, E<:AbstractSpectralRedistribution{FT}} <: AbstractOceanInelasticProcess{FT}
    "CDOM absorber driving the fluorescence"
    cdom::C
    "Spectral redistribution `f(λ', λ)` `[1/nm]`"
    emission::E
    "Quantum yield `φ_g` (dimensionless, band-integrated)"
    quantum_yield::FT
    "Excitation-wavelength integration window `[nm]`"
    excitation_range::Tuple{FT, FT}

    function CDOMFluorescence{FT, C, E}(cdom::C, emission::E, quantum_yield,
                                        excitation_range) where {FT, C, E}
        _check_yield_and_range(quantum_yield, excitation_range, "CDOMFluorescence")
        new{FT, C, E}(cdom, emission, quantum_yield, excitation_range)
    end
end

function CDOMFluorescence(cdom::CDOM{FT}; emission::AbstractSpectralRedistribution{FT},
                          quantum_yield::Real,
                          excitation_range::Tuple{<:Real, <:Real} = (FT(300), FT(400))) where {FT}
    return CDOMFluorescence{FT, typeof(cdom), typeof(emission)}(
        cdom, emission, FT(quantum_yield), (FT(excitation_range[1]), FT(excitation_range[2])))
end

# Trait implementations
excitation_absorption(p::CDOMFluorescence, λ_prime::Real) =
    p.quantum_yield * absorption(p.cdom, λ_prime)

emission(p::CDOMFluorescence, λ_prime::Real, λ::Real) =
    redistribution(p.emission, λ_prime, λ)

excitation_range(p::CDOMFluorescence)    = p.excitation_range
excitation_pathways(p::CDOMFluorescence) = (ExcitationPathway(p.cdom),)
is_isotropic(::CDOMFluorescence)         = true
conversion_units(::CDOMFluorescence)     = PhotonConversion()
quantum_yield(p::CDOMFluorescence)       = p.quantum_yield

function Base.show(io::IO, p::CDOMFluorescence{FT}) where {FT}
    print(io, "CDOMFluorescence{$FT}(φ=", p.quantum_yield, ", ", p.emission, ")")
end
