# =============================================================================
# Fluorescence processes
# =============================================================================
#
# Ocean fluorescence sources covered here:
#
#   IsotropicFluorescence  — generic isotropic fluorescence parameterized by
#                            an `AbstractAbsorber` (phytoplankton, NAP, …),
#                            a spectral-redistribution profile, and a
#                            quantum yield. Chlorophyll-a fluorescence is
#                            the primary use case.
#   CDOMFluorescence       — convenience wrapper specialized on a `CDOM`
#                            absorber; same physics, narrower type.
#
# The chlorophyll-fluorescence parameterization matches Fell (1997) §2.2.9
# and Gordon (1979): Gaussian emission centered at 685 nm with σ = 10.6 nm
# (FWHM ≈ 25 nm), quantum yield φ_C ≈ 0.003 over the 400–700 nm
# photosynthetically-active excitation range.
# =============================================================================

"""
$(TYPEDEF)

Isotropic fluorescence from an `AbstractAbsorber` with a configurable
spectral redistribution and quantum yield.

# Fields
$(TYPEDFIELDS)

# Examples

```julia
phyto = Phytoplankton{Float64, Bricaud1995}(Chl = 0.5)
sif   = chlorophyll_fluorescence(phyto)                   # default config
sif   = IsotropicFluorescence{Float64}(phyto;
            emission        = GaussianEmission{Float64}(center_nm = 685.0,
                                                        σ_nm      = 10.6),
            quantum_yield   = 0.003,
            excitation_range = (400.0, 700.0))
```
"""
struct IsotropicFluorescence{FT, A<:AbstractAbsorber{FT},
                                 E<:AbstractSpectralRedistribution{FT}} <: AbstractOceanInelasticProcess{FT}
    "Absorbing constituent whose absorbed photons drive emission"
    absorber::A
    "Spectral redistribution `f(λ', λ)` `[1/nm]`"
    emission::E
    "Quantum yield `φ` (dimensionless, typically 0.001–0.01 for Chl-a)"
    quantum_yield::FT
    "Excitation-wavelength integration window `[nm]`"
    excitation_range::Tuple{FT, FT}

    function IsotropicFluorescence{FT, A, E}(absorber::A, emission::E,
                                              quantum_yield::FT,
                                              excitation_range::Tuple{FT, FT}) where {
            FT<:AbstractFloat, A<:AbstractAbsorber{FT},
            E<:AbstractSpectralRedistribution{FT}}
        zero(FT) ≤ quantum_yield ≤ one(FT) ||
            throw(ArgumentError("IsotropicFluorescence: quantum_yield must be in [0, 1], got $quantum_yield"))
        excitation_range[1] < excitation_range[2] ||
            throw(ArgumentError("IsotropicFluorescence: excitation_range must be (lo < hi), got $excitation_range"))
        new{FT, A, E}(absorber, emission, quantum_yield, excitation_range)
    end
end

"""
    IsotropicFluorescence{FT}(absorber; emission, quantum_yield, excitation_range)

Construct an `IsotropicFluorescence` process from an `AbstractAbsorber`.
`absorber` must be typed `<:AbstractAbsorber{FT}`, which — after the
hierarchy refactor in §4.1 — includes all `AbstractAbsorbingScatterer`
subtypes (`PureWater`, `Phytoplankton`, `NonAlgalParticles`).
"""
function IsotropicFluorescence{FT}(absorber::AbstractAbsorber{FT};
        emission::AbstractSpectralRedistribution{FT},
        quantum_yield                         = FT(0.003),
        excitation_range::Tuple{<:Real, <:Real} = (FT(400), FT(700)),
) where {FT<:AbstractFloat}
    return IsotropicFluorescence{FT, typeof(absorber), typeof(emission)}(
        absorber, emission, FT(quantum_yield),
        (FT(excitation_range[1]), FT(excitation_range[2])))
end

# -----------------------------------------------------------------------------
# Trait implementations
# -----------------------------------------------------------------------------

# a^I(λ') = φ · a_absorber(λ'): quantum-yield-weighted excitation absorption
excitation_absorption(p::IsotropicFluorescence, λ_prime::Real) =
    p.quantum_yield * absorption(p.absorber, λ_prime)

emission(p::IsotropicFluorescence, λ_prime::Real, λ::Real) =
    redistribution(p.emission, λ_prime, λ)

excitation_range(p::IsotropicFluorescence) = p.excitation_range
is_isotropic(::IsotropicFluorescence)      = true
quantum_yield(p::IsotropicFluorescence)    = p.quantum_yield

function Base.show(io::IO, p::IsotropicFluorescence{FT}) where {FT}
    print(io, "IsotropicFluorescence{$FT}(", typeof(p.absorber).name.name,
              ", φ=", p.quantum_yield, ", ", p.emission, ")")
end

# -----------------------------------------------------------------------------
# Convenience constructor: chlorophyll fluorescence per Fell (1997)
# -----------------------------------------------------------------------------

"""
    chlorophyll_fluorescence(phyto::Phytoplankton; peak_nm=685.0, σ_nm=10.6,
                             quantum_yield=0.003,
                             excitation_range=(400.0, 700.0))

Construct the standard oceanic chlorophyll-a fluorescence process: a
Gaussian emission band at 685 nm with σ ≈ 10.6 nm (FWHM ≈ 25 nm) and
quantum yield φ_C = 0.003 — the parameterization used in Fell (1997)
§2.2.9 (North Sea average; cites Fischer & Kronfeld 1990, Deep-Sea Res.
37, 865) and Gordon (1979, Appl. Opt. 18, 1161).
"""
function chlorophyll_fluorescence(phyto::Phytoplankton{FT};
        peak_nm                                 = FT(685.0),
        σ_nm                                    = FT(10.6),
        quantum_yield                           = FT(0.003),
        excitation_range::Tuple{<:Real, <:Real} = (FT(400), FT(700)),
) where {FT<:AbstractFloat}
    em = GaussianEmission{FT}(center_nm = FT(peak_nm), σ_nm = FT(σ_nm))
    return IsotropicFluorescence{FT}(phyto;
        emission         = em,
        quantum_yield    = quantum_yield,
        excitation_range = excitation_range)
end

# =============================================================================
# CDOM fluorescence
# =============================================================================

"""
$(TYPEDEF)

CDOM fluorescence as an isotropic inelastic process. Structurally equivalent
to [`IsotropicFluorescence`](@ref) specialized on a [`CDOM`](@ref) absorber,
but carried as a distinct type so solvers can dispatch separately and
users can read model intent directly from the type.

Typical excitation/emission footprint per Coble (1996, 2007) and
Bricaud-Morel-Prieur (1981): broad CDOM emission around 420–450 nm for
UV excitation (300–380 nm), with quantum yield 1–2 % when integrated
over the band. Fell (1997) uses a simple Gaussian model.

# Fields
$(TYPEDFIELDS)
"""
struct CDOMFluorescence{FT, E<:AbstractSpectralRedistribution{FT}} <: AbstractOceanInelasticProcess{FT}
    "CDOM absorber driving the fluorescence"
    cdom::CDOM{FT}
    "Spectral redistribution `f(λ', λ)` `[1/nm]`"
    emission::E
    "Quantum yield `φ_g` (dimensionless, band-integrated)"
    quantum_yield::FT
    "Excitation-wavelength integration window `[nm]`"
    excitation_range::Tuple{FT, FT}
end

"""
    CDOMFluorescence{FT}(cdom; emission, quantum_yield=0.01, excitation_range=(300,400))
"""
function CDOMFluorescence{FT}(cdom::CDOM{FT};
        emission::AbstractSpectralRedistribution{FT},
        quantum_yield                           = FT(0.01),
        excitation_range::Tuple{<:Real, <:Real} = (FT(300), FT(400)),
) where {FT<:AbstractFloat}
    zero(FT) ≤ FT(quantum_yield) ≤ one(FT) ||
        throw(ArgumentError("CDOMFluorescence: quantum_yield must be in [0, 1], got $quantum_yield"))
    excitation_range[1] < excitation_range[2] ||
        throw(ArgumentError("CDOMFluorescence: excitation_range must be (lo < hi)"))
    return CDOMFluorescence{FT, typeof(emission)}(
        cdom, emission, FT(quantum_yield),
        (FT(excitation_range[1]), FT(excitation_range[2])))
end

# Trait implementations
excitation_absorption(p::CDOMFluorescence, λ_prime::Real) =
    p.quantum_yield * absorption(p.cdom, λ_prime)

emission(p::CDOMFluorescence, λ_prime::Real, λ::Real) =
    redistribution(p.emission, λ_prime, λ)

excitation_range(p::CDOMFluorescence) = p.excitation_range
is_isotropic(::CDOMFluorescence)      = true
quantum_yield(p::CDOMFluorescence)    = p.quantum_yield

function Base.show(io::IO, p::CDOMFluorescence{FT}) where {FT}
    print(io, "CDOMFluorescence{$FT}(φ=", p.quantum_yield, ", ", p.emission, ")")
end
