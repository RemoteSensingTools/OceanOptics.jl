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
#   AbstractPureWaterAbsorptionModel   # parameterization namespace
#       ├── IOCCG2018
#       ├── PopeFry1997
#       ├── Kou1993
#       └── SmithBaker1981
#
#   AbstractPhytoplanktonModel
#       ├── Bricaud1998              # a_p = A·Chl^E
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
#     absorption(p::Phytoplankton{FT, OtherModel},   λ)  # different method
#
# This means swapping a parameterization is a single-symbol change at the
# call site, not a rewrite of the struct.
# =============================================================================

"""
$(TYPEDEF)

Root of the ocean-constituent type hierarchy. Any material present in the
water column that contributes to absorption, scattering, or emission subtypes
this.

`FT` is the element type of the constituent's parameters and results: a floating-point
type, or any other `Real` that promotes with `Float64` such as a `ForwardDiff.Dual`,
which carries derivatives with respect to the parameters through the scalar functions
and [`layer_optics`](@ref). The constituents of one layer share one `FT`.
"""
abstract type AbstractOceanConstituent{FT<:Real} end

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
in `[1/m]`, plus a phase function accessible via [`scattering_phase`](@ref).
"""
abstract type AbstractScatterer{FT} <: AbstractOceanConstituent{FT} end

"""
$(TYPEDEF)

Constituents that both absorb and scatter (phytoplankton, non-algal
particles, pure water). Must implement both [`absorption`](@ref) and
[`scattering`](@ref), and provide a [`scattering_phase`](@ref).

The type subtypes [`AbstractAbsorber`](@ref), so an absorbing scatterer is accepted
wherever an absorber is required. It is not an [`AbstractScatterer`](@ref): the
scattering side is reached through the methods of each concrete type.
"""
abstract type AbstractAbsorbingScatterer{FT} <: AbstractAbsorber{FT} end

# =============================================================================
# Parameterization namespaces — model tags
# =============================================================================
# These are empty singleton types used as type parameters on the concrete
# constituent structs. They exist solely to drive multiple dispatch so that
# `absorption(p::Phytoplankton{FT, Bricaud1998}, λ)` and
# `absorption(p::Phytoplankton{FT, OtherModel},  λ)` resolve to different
# methods with zero runtime branching.

# Pure water carries two distinct model axes (absorption + scattering) and
# therefore lives under its own per-axis abstract types
# (`AbstractPureWaterAbsorptionModel`, `AbstractPureWaterScatteringModel`)
# declared in `src/materials/pure_water.jl`.

"""
$(TYPEDEF)

Parameterization-model tag for phytoplankton IOPs.
"""
abstract type AbstractPhytoplanktonModel end

"""
$(TYPEDEF)

Parameterization-model tag for CDOM absorption.
"""
abstract type AbstractCDOMModel end

"""
$(TYPEDEF)

Parameterization-model tag for non-algal particle (NAP) IOPs.
"""
abstract type AbstractNAPModel end

# =============================================================================
# Generic contract (traits): docstrings for the functions we expect subtypes
# to implement. Each concrete file supplies its method.
# =============================================================================

"""
    absorption(c::AbstractAbsorber, λ) -> FT

Return the absorption coefficient `a(λ)` `[1/m]` of constituent `c`.

`λ` is a wavelength in nm, or a vector of wavelengths.

# Examples
```jldoctest
julia> absorption(CDOM(), 440.0)
0.03
```
"""
function absorption end

"""
    scattering(c::AbstractScatterer, λ) -> FT

Return the scattering coefficient `b(λ)` `[1/m]` of constituent `c`.

`λ` is a wavelength in nm, or a vector of wavelengths.
"""
function scattering end

"""
    backscattering(c::AbstractScatterer, λ) -> FT

Return the backscattering coefficient `b_b(λ)` `[1/m]` of constituent `c`.

It is the integral of the volume scattering function over the back hemisphere,
`b_b = b(λ) · B` with `B` the backscatter fraction of the phase function.
"""
function backscattering end

"""
    attenuation(c::AbstractOceanConstituent, λ) -> FT

Return the beam attenuation coefficient `c(λ) = a(λ) + b(λ)` `[1/m]`.

The method for absorbing scatterers sums [`absorption`](@ref) and [`scattering`](@ref).
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
    scattering_phase(c::AbstractScatterer) -> AbstractOceanPhaseFunction

Return the phase function of scatterer `c`.

The phase function owns `p(μ)` and its Legendre moments; the scatterer names which one
applies.
"""
function scattering_phase end

# -----------------------------------------------------------------------------
# Absorption parts
# -----------------------------------------------------------------------------
#
# An absorber may attribute its absorption to parts — phytoplankton with a
# pigment partition splits a_ph into chlorophyll, carotenoid, … parts — so
# that a solver can report the light absorbed by each and a fluorescence
# process can name the parts that excite it. `layer_optics` makes one
# absorption component per part.

"""
    absorption_parts(c::AbstractAbsorber) -> Tuple{Vararg{Symbol}}

Return the names of the parts into which `c` attributes its absorption.

The order is that of [`absorption_by_part`](@ref). The result is `(:total,)` unless `c`
splits its absorption, as [`Phytoplankton`](@ref) does with a pigment partition.
"""
absorption_parts(::AbstractAbsorber) = (:total,)

"""
    absorption_by_part(c::AbstractAbsorber, λ) -> Tuple

Return the absorption coefficient `[1/m]` of each part of `c` at `λ` `[nm]`.

The parts are those of [`absorption_parts`](@ref)`(c)`; they sum to
[`absorption`](@ref)`(c, λ)`.
"""
absorption_by_part(c::AbstractAbsorber, λ::Real) = (absorption(c, λ),)

"""
    absorption_part(c::AbstractAbsorber, part::Symbol, λ) -> FT

Return the absorption coefficient `[1/m]` of one part of `c` at `λ` `[nm]`.
"""
function absorption_part(c::AbstractAbsorber, part::Symbol, λ::Real)
    i = findfirst(==(part), absorption_parts(c))
    i === nothing &&
        throw(ArgumentError("$(nameof(typeof(c))) has no absorption part :$part; its parts " *
                            "are $(absorption_parts(c))"))
    return absorption_by_part(c, λ)[i]
end

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

Root of the inelastic-source hierarchy.

A concrete subtype ([`IsotropicFluorescence`](@ref), [`WaterRaman`](@ref),
[`CDOMFluorescence`](@ref)) describes one excitation-to-emission channel and implements
the primitives of the factorization of Fell (1997): [`excitation_absorption`](@ref)
`(p, λ')`, [`emission`](@ref)`(p, λ', λ)`, [`excitation_range`](@ref)`(p)`,
[`is_isotropic`](@ref)`(p)` and [`conversion_units`](@ref)`(p)`, which states whether the
excitation coefficient counts photons or energy.
"""
abstract type AbstractOceanInelasticProcess{FT<:Real} end

# =============================================================================
# Precision of scalar results
# =============================================================================
#
# A model holds its parameters in its element type `FT`; data tables and
# fitted constants are `Float64`. A scalar function whose formula contains
# such constants evaluates in at least double precision and rounds its
# result once to the type that the model and the query argument imply, so
# that a `Float32` model queried at a `Float32` wavelength returns `Float32`
# with the rounding error of one conversion.

"""
    _result_type(FT, x...)

Return the type of a scalar result of a model with element type `FT` queried at `x...`.

It is `FT` promoted with the floating-point type of each argument.
Dual numbers and floats wider than `FT` therefore pass through.
"""
_result_type(::Type{FT}, x) where {FT} = promote_type(FT, typeof(float(x)))
_result_type(::Type{FT}, x, y) where {FT} = promote_type(FT, typeof(float(x)), typeof(float(y)))

"""
    _wide(x)

Return `x` in at least double precision.

`Float16` and `Float32` are widened to `Float64` (exactly); any other number, a dual
number included, is returned as it is and promotes with the `Float64` constants of the
formula.
"""
_wide(x::Union{Float16, Float32}) = Float64(x)
_wide(x::Real) = x

"""
    _within(x, lo, hi) -> Bool

Return `true` when the value of `x` lies in the closed interval `[lo, hi]`.

For a floating-point `x` this is `lo ≤ x ≤ hi`. Any other real `x` (an integer, a dual
number) is compared with the floating-point neighbours `prevfloat(lo) < x < nextfloat(hi)`:
ForwardDiff 1.x orders dual numbers lexicographically by value and then partials, so
`Dual(35, 1) ≤ 35` is `false` and `≤` would reject a dual at an upper bound; a strict
comparison with the neighbour is decided by the value alone. For a dual whose value
has more precision than `lo` and `hi`, the interval is widened by one unit in the last
place of the bounds.
"""
_within(x::AbstractFloat, lo, hi) = lo ≤ x ≤ hi
_within(x::Real, lo, hi) = prevfloat(float(lo)) < x < nextfloat(float(hi))

"""
    _check_parameter_type(FT)

Throw an `ArgumentError` when `FT` is an integer or rational type.

Models and optical properties are computed in a floating-point type or in a dual number
of one. An integer type fails at the first fractional default, a rational type would
compute in exact rationals.
"""
function _check_parameter_type(::Type{FT}) where {FT}
    FT <: Union{Integer, Rational} &&
        throw(ArgumentError("element type $FT: models and optical properties need a " *
                            "floating-point type (Float64, Float32) or a dual number of " *
                            "one; use $(float(FT)) instead"))
    return nothing
end

"""
    _parameter_type(parameters::NamedTuple)

Return the element type of a model constructed without an explicit one.

It is `Float64` promoted with the type of every real-valued parameter, so that plain
numbers of any precision give a `Float64` model and a dual number (or a `BigFloat`)
gives a model that carries it.
"""
_parameter_type(parameters::NamedTuple) = _promote_reals(Float64, values(parameters)...)
_promote_reals(::Type{T}) where {T} = T
_promote_reals(::Type{T}, x::Real, rest...) where {T} =
    _promote_reals(promote_type(T, typeof(x)), rest...)
_promote_reals(::Type{T}, x, rest...) where {T} = _promote_reals(T, rest...)
