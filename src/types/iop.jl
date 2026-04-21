"""
$(TYPEDEF)

Bundle of inherent optical properties at one wavelength (or one wavelength
grid). This is the neutral currency of the package: every constituent
produces an `IOP`, constituents mix via `+`, and layers reduce to a single
`IOP` before being expanded into an [`OceanLayerOptics`](@ref).

An `IOP` is **wavelength-local**: the scalar-field version describes one λ;
the vector-field version describes a λ grid. Same struct, parametric fields.

# Fields
$(TYPEDFIELDS)

# Example

```julia
iop1 = IOP(0.1, 0.05, 0.018)       # pure water at 550 nm
iop2 = IOP(0.03, 0.30, 0.012)      # phytoplankton contribution
total = iop1 + iop2                # correctly energy-weights backscatter
```
"""
struct IOP{FT<:AbstractFloat, A<:Union{FT, AbstractArray{FT}}}
    "Absorption coefficient `a(λ)` `[1/m]`"
    a::A
    "Scattering coefficient `b(λ)` `[1/m]`"
    b::A
    "Backscattering coefficient `b_b(λ)` `[1/m]`. Optional: zero for pure absorbers."
    bb::A
end

# Lightweight constructors ----------------------------------------------------

"Scalar IOP for one wavelength."
IOP(a::FT, b::FT, bb::FT) where {FT<:AbstractFloat} = IOP{FT, FT}(a, b, bb)

"Per-wavelength-grid IOP backed by identically-typed vectors/matrices."
IOP(a::A, b::A, bb::A) where {FT<:AbstractFloat, A<:AbstractArray{FT}} =
    IOP{FT, A}(a, b, bb)

"Mixed-type convenience: promote to a common `AbstractFloat`."
function IOP(a::Real, b::Real, bb::Real)
    FT = promote_type(typeof(a), typeof(b), typeof(bb))
    FT <: AbstractFloat || (FT = float(FT))
    return IOP(FT(a), FT(b), FT(bb))
end

"Pure absorber: no scattering, no backscatter."
IOP_absorber(a) = IOP(a, zero(a), zero(a))

"Pure scatterer with backscatter fraction `B = bb/b`."
IOP_scatterer(b, B) = IOP(zero(b), b, b * B)

# Accessors -------------------------------------------------------------------

"Beam attenuation `c = a + b`."
attenuation(iop::IOP) = iop.a .+ iop.b

"Single-scattering albedo `ϖ = b / (a + b)`; zero when `a = b = 0`."
function single_scattering_albedo(iop::IOP{FT}) where {FT}
    c = attenuation(iop)
    # Broadcast-safe form; works for scalar or array IOPs.
    return @. ifelse(c > zero(FT), iop.b / c, zero(FT))
end

"Backscatter fraction `B = b_b / b`; zero when `b = 0`."
function backscatter_fraction(iop::IOP{FT}) where {FT}
    return @. ifelse(iop.b > zero(FT), iop.bb / iop.b, zero(FT))
end

# Algebra ---------------------------------------------------------------------
#
# IOPs compose additively — this is the physical mixing law for independent
# constituents, and matches vSmartMOM's `+(::CoreScatteringOpticalProperties)`.
# All three quantities (a, b, bb) are *extensive* (proportional to density),
# so they literally sum. Phase functions mix with scattering-weighted averages
# — that is handled one level up, in `OceanLayerOptics` assembly, because it
# is an operation on *functions*, not scalars.

Base.:+(x::IOP, y::IOP) = IOP(x.a .+ y.a, x.b .+ y.b, x.bb .+ y.bb)
Base.zero(iop::IOP) = IOP(zero(iop.a), zero(iop.b), zero(iop.bb))

"Multiply all IOPs by a scalar (e.g. to scale by layer thickness)."
Base.:*(s::Real, iop::IOP) = IOP(s .* iop.a, s .* iop.b, s .* iop.bb)
Base.:*(iop::IOP, s::Real) = s * iop

function Base.show(io::IO, iop::IOP{FT}) where {FT}
    a = iop.a isa AbstractArray ? "[$(length(iop.a)) λ]" : string(iop.a)
    b = iop.b isa AbstractArray ? "[$(length(iop.b)) λ]" : string(iop.b)
    bb = iop.bb isa AbstractArray ? "[$(length(iop.bb)) λ]" : string(iop.bb)
    print(io, "IOP{$FT}(a=$a, b=$b, bb=$bb)")
end
