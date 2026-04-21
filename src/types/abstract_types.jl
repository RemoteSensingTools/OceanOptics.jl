# =============================================================================
# Abstract type hierarchy
# =============================================================================
#
# Design philosophy
# -----------------
# Three orthogonal axes describe what a water constituent is:
#
#     OPTICAL ROLE      :  absorbs? scatters? emits (fluorescence)?
#     BIOGEOCHEMICAL ID :  is it water? phytoplankton? CDOM? a mineral?
#     MODEL CHOICE      :  which parameterization (Morel, Bricaud, Mobley-Case1, ...)?
#
# The plan in `docs/dev_notes/ocean.md` conflates these by subtyping a single
# `AbstractOceanConstituent`. That forces every new parameterization (say, a
# second phytoplankton absorption model) to become a sibling of
# `PureWater`, `CDOM`, `Phytoplankton`, which is not the physical structure.
#
# Instead we split:
#
#   AbstractOceanConstituent
#       ├── AbstractAbsorber        # defines absorption(·, λ)
#       ├── AbstractScatterer       # defines scattering(·, λ) + phase function
#       └── AbstractAbsorbingScatterer <: {AbstractAbsorber, AbstractScatterer}
#
#   AbstractPureWaterModel          # parameterization namespace
#       ├── SmithBaker1981
#       ├── PopeFry1997
#       └── MasonConeFry2016         # UV extension
#
#   AbstractPhytoplanktonModel
#       ├── Gordon1992
#       ├── Bricaud1998              # aph* power law
#       └── MobleyNewCase1            # HydroLight 5 default
#
# A concrete `Phytoplankton` struct carries a *model* as a type parameter, so
# absorption/scattering dispatch hits the right parameterization without
# branches or if-chains:
#
#     struct Phytoplankton{FT, M<:AbstractPhytoplanktonModel} <: AbstractAbsorbingScatterer
#         model::M
#         Chl::FT
#     end
#     absorption(p::Phytoplankton{FT, Bricaud1998}, λ)  # dispatches on M
#     absorption(p::Phytoplankton{FT, Gordon1992},   λ)  # different method
#
# This means swapping a parameterization is a single-symbol change at the
# call site, not a rewrite of the struct.
# =============================================================================

"""
$(TYPEDEF)

Root of the ocean-constituent type hierarchy. Any material present in the
water column that contributes to absorption, scattering, or emission subtypes
this.
"""
abstract type AbstractOceanConstituent{FT<:AbstractFloat} end

"""
$(TYPEDEF)

Trait for constituents that absorb light. Concrete types must implement
[`absorption`](@ref) `(c, λ) -> FT` returning an absorption coefficient
in `[1/m]`.
"""
abstract type AbstractAbsorber{FT} <: AbstractOceanConstituent{FT} end

"""
$(TYPEDEF)

Trait for constituents that scatter light. Concrete types must implement
[`scattering`](@ref) `(c, λ) -> FT` returning a scattering coefficient
in `[1/m]`, plus a phase function accessible via [`phase_function`](@ref).
"""
abstract type AbstractScatterer{FT} <: AbstractOceanConstituent{FT} end

"""
$(TYPEDEF)

Constituents that both absorb and scatter (phytoplankton, non-algal
particles, pure water). Must implement both [`absorption`](@ref) and
[`scattering`](@ref), and provide a [`phase_function`](@ref).

Subtypes [`AbstractAbsorber`](@ref): any `AbstractAbsorbingScatterer` is
dispatchable wherever an `AbstractAbsorber` is required. DESIGN §4.1
called for joint inheritance from both `AbstractAbsorber` and
`AbstractScatterer`, which Julia's single-inheritance model cannot
express; we pick `AbstractAbsorber` as the parent because the bulk of
downstream code (inelastic source operators in §7, constituent-mixing
helpers) dispatches through absorption-centric slots. Scatterer-centric
behavior is still reachable via specific methods on each concrete
absorbing-scatterer type.
"""
abstract type AbstractAbsorbingScatterer{FT} <: AbstractAbsorber{FT} end

# =============================================================================
# Parameterization namespaces — model tags
# =============================================================================
# These are empty singleton types used as type parameters on the concrete
# constituent structs. They exist solely to drive multiple dispatch so that
# `absorption(p::Phytoplankton{FT, Bricaud1998}, λ)` and
# `absorption(p::Phytoplankton{FT, Gordon1992},  λ)` resolve to different
# methods with zero runtime branching.

# Pure water carries two distinct model axes (absorption + scattering) and
# therefore lives under its own per-axis abstract types
# (`AbstractPureWaterAbsorptionModel`, `AbstractPureWaterScatteringModel`)
# declared in `src/materials/pure_water.jl`.

"Parameterization-model tag for phytoplankton IOPs."
abstract type AbstractPhytoplanktonModel end

"Parameterization-model tag for CDOM absorption."
abstract type AbstractCDOMModel end

"Parameterization-model tag for non-algal particle (NAP) IOPs."
abstract type AbstractNAPModel end

# =============================================================================
# Generic contract (traits): docstrings for the functions we expect subtypes
# to implement. Each concrete file supplies its method.
# =============================================================================

"""
    absorption(c::AbstractAbsorber, λ) -> FT

Spectral absorption coefficient `a(λ)` `[1/m]` contributed by constituent `c`.

`λ` may be a plain number (interpreted as nm), a `Unitful.Quantity` with
spectral units (nm, µm, cm⁻¹, …), or an `AbstractVector` of either.
"""
function absorption end

"""
    scattering(c::AbstractScatterer, λ) -> FT

Spectral scattering coefficient `b(λ)` `[1/m]`.
"""
function scattering end

"""
    backscattering(c::AbstractScatterer, λ) -> FT

Backscattering coefficient `b_b(λ)` `[1/m]` — the integral of the volume
scattering function over the back hemisphere. Equals `b(λ) * B` where `B`
is the backscatter fraction determined by the phase function.
"""
function backscattering end

"""
    attenuation(c::AbstractOceanConstituent, λ) -> FT

Beam attenuation `c(λ) = a(λ) + b(λ)` `[1/m]`. Default implementation sums
`absorption` and `scattering`; specific constituents (pure absorbers,
pure scatterers) should override.
"""
attenuation(c::AbstractAbsorbingScatterer, λ) = absorption(c, λ) + scattering(c, λ)
attenuation(c::AbstractAbsorber,           λ) = absorption(c, λ)
attenuation(c::AbstractScatterer,          λ) = scattering(c, λ)

# The default scattering/absorption for a pure absorber / pure scatterer is 0.
# Keeps the `iop()` summation branch-free.
scattering(::AbstractAbsorber{FT},  λ) where {FT} = zero(FT)
absorption(::AbstractScatterer{FT}, λ) where {FT} = zero(FT)
backscattering(::AbstractAbsorber{FT}, λ) where {FT} = zero(FT)

"""
    phase_function(c::AbstractScatterer) -> AbstractOceanPhaseFunction

The phase function associated with this scatterer. Returns an
[`AbstractOceanPhaseFunction`](@ref) — the phase function itself owns the
knowledge of what `p(μ)` and its Legendre moments are; the scatterer just
says which one to use.
"""
function phase_function end

# =============================================================================
# Inelastic processes (fluorescence, Raman) — abstract root
# =============================================================================
#
# Defined in the types spine so that `OceanLayer` can carry a
# `fluorophores::Vector{<:AbstractOceanInelasticProcess{FT}}` field at its
# own include time. Concrete subtypes and their kernel functions live in
# `src/inelastic/`.

"""
$(TYPEDEF)

Root of the inelastic-source hierarchy (DESIGN §4.4 / §7.1). Concrete
subtypes — `IsotropicFluorescence`, `WaterRaman`, `CDOMFluorescence` —
describe one excitation-to-emission channel apiece and implement the
Fell (1997) factorization primitives `excitation_absorption(p, λ')`,
`emission(p, λ', λ)`, `excitation_range(p)`, and `is_isotropic(p)`.
"""
abstract type AbstractOceanInelasticProcess{FT<:AbstractFloat} end
