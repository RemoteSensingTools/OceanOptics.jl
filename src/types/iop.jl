"""
$(TYPEDEF)

Absorption, scattering and backscattering coefficients `[1/m]` at one wavelength (number
fields) or on a wavelength grid (array fields).

`IOP`s of constituents in one volume add with `+` and scale with `*`. The type is a
standalone value: [`layer_optics`](@ref) does not use it and reports optical depths per
component in an [`OceanLayerOptics`](@ref).

# Fields
$(TYPEDFIELDS)

# Example

```julia
iop1 = IOP(0.1, 0.05, 0.018)       # pure water at 550 nm
iop2 = IOP(0.03, 0.30, 0.012)      # phytoplankton contribution
total = iop1 + iop2                # a, b and bb add; B = bb/b is scattering-weighted
```
"""
struct IOP{FT<:Real, A<:Union{FT, AbstractArray{FT}}}
    "Absorption coefficient `a(λ)` `[1/m]`"
    a::A
    "Scattering coefficient `b(λ)` `[1/m]`"
    b::A
    "Backscattering coefficient `b_b(λ)` `[1/m]`. Optional: zero for pure absorbers."
    bb::A
end

# Lightweight constructors ----------------------------------------------------

"""
    IOP(a, b, bb)

Construct an `IOP` from absorption `a`, scattering `b` and backscattering `bb` `[1/m]`.

Three numbers describe one wavelength and are promoted to a common floating-point type (dual numbers are kept); three
arrays of one type describe a wavelength grid.
"""
function IOP(a::FT, b::FT, bb::FT) where {FT<:Real}
    T = float(FT)
    return IOP{T, T}(T(a), T(b), T(bb))
end

IOP(a::A, b::A, bb::A) where {FT<:Real, A<:AbstractArray{FT}} =
    IOP{FT, A}(a, b, bb)
# Integer and rational arrays are converted to floating point, as the scalar method does.
IOP(a::A, b::A, bb::A) where {FT<:Union{Integer, Rational}, A<:AbstractArray{FT}} =
    IOP(float.(a), float.(b), float.(bb))

function IOP(a::Real, b::Real, bb::Real)
    FT = float(promote_type(typeof(a), typeof(b), typeof(bb)))
    return IOP(FT(a), FT(b), FT(bb))
end

"Pure absorber: no scattering, no backscatter."
IOP_absorber(a) = IOP(a, zero(a), zero(a))

"Pure scatterer with backscatter fraction `B = bb/b`."
IOP_scatterer(b, B) = IOP(zero(b), b, b * B)

# Accessors -------------------------------------------------------------------

"""
    attenuation(iop::IOP)

Return the beam attenuation `c = a + b` `[1/m]`.
"""
attenuation(iop::IOP) = iop.a .+ iop.b

"""
    single_scattering_albedo(iop::IOP)

Return the single-scattering albedo `ϖ = b / (a + b)`, zero when `a = b = 0`.
"""
function single_scattering_albedo(iop::IOP{FT}) where {FT}
    c = attenuation(iop)
    # Broadcast-safe form; works for scalar or array IOPs.
    return @. ifelse(c > zero(FT), iop.b / c, zero(FT))
end

"""
    backscatter_fraction(iop::IOP)

Return the backscatter fraction `B = b_b / b`, zero when `b = 0`.
"""
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

"""
    *(s::Real, iop::IOP)

Multiply the three coefficients of `iop` by the number `s`.
"""
Base.:*(s::Real, iop::IOP) = IOP(s .* iop.a, s .* iop.b, s .* iop.bb)
Base.:*(iop::IOP, s::Real) = s * iop

function Base.show(io::IO, iop::IOP{FT}) where {FT}
    a = iop.a isa AbstractArray ? "[$(length(iop.a)) λ]" : string(iop.a)
    b = iop.b isa AbstractArray ? "[$(length(iop.b)) λ]" : string(iop.b)
    bb = iop.bb isa AbstractArray ? "[$(length(iop.bb)) λ]" : string(iop.bb)
    print(io, "IOP{$FT}(a=$a, b=$b, bb=$bb)")
end
