# =============================================================================
# Refractive index of seawater at the air-water interface
# =============================================================================
#
# The interface index is the real index of seawater relative to air (what
# Snell's law and the Fresnel coefficients need at the sea surface). Two
# models:
#
#   FixedRefractiveIndex(n)          — benchmarks (e.g. n = 1.34)
#   QuanFry1995(; temperature, salinity)
#                                    — Quan & Fry (1995), Appl. Opt. 34,
#                                      3477: relative to air, fitted to
#                                      400–700 nm, 0–30 °C, 0–35 PSU
#
# Outside its fitted domain QuanFry1995 throws unless `extrapolate = true`
# is chosen explicitly.
# =============================================================================

"""
$(TYPEDEF)

Model for the real refractive index of seawater relative to air; see
[`refractive_index`](@ref).
"""
abstract type AbstractRefractiveIndexModel end

"""
$(TYPEDEF)

Wavelength-independent refractive index of seawater relative to air,
`n ≥ 1` (`n = 1` is an index-matched interface).

# Fields
$(TYPEDFIELDS)
"""
struct FixedRefractiveIndex{FT<:Real} <: AbstractRefractiveIndexModel
    "Index relative to air"
    n::FT

    function FixedRefractiveIndex{FT}(n::FT) where {FT<:Real}
        isfinite(n) && n ≥ 1 ||
            throw(ArgumentError("FixedRefractiveIndex: the index of seawater relative to air " *
                                "must be finite and ≥ 1, got n = $n; check that it is " *
                                "n_water / n_air and not its inverse"))
        new{FT}(n)
    end
end
FixedRefractiveIndex(n::Real) = FixedRefractiveIndex{float(typeof(n))}(float(n))

"""
$(TYPEDEF)

Quan & Fry (1995) empirical refractive index of seawater relative to air,

```
n = n₀ + (n₁ + n₂T + n₃T²) S + n₄T² + (n₅ + n₆S + n₇T)/λ + n₈/λ² + n₉/λ³
```

with `λ` in nm, `T` in °C and `S` in PSU, fitted to Austin & Halikas (1976)
data over 400–700 nm, 0–30 °C and 0–35 PSU.

# Fields
$(TYPEDFIELDS)

Queries outside the fitted domain throw unless `extrapolate = true`.
"""
struct QuanFry1995{FT<:Real} <: AbstractRefractiveIndexModel
    "Temperature `[°C]`"
    temperature::FT
    "Salinity `[PSU]`"
    salinity::FT
    "Allow use outside 400–700 nm, 0–30 °C, 0–35 PSU"
    extrapolate::Bool
end

function QuanFry1995(; temperature::Real, salinity::Real, extrapolate::Bool = false)
    FT = float(promote_type(typeof(temperature), typeof(salinity)))
    extrapolate || (_within(temperature, 0, 30) && _within(salinity, 0, 35)) ||
        throw(DomainError((temperature, salinity),
                          "QuanFry1995 is fitted to 0–30 °C and 0–35 PSU; pass extrapolate = true to use it outside"))
    return QuanFry1995{FT}(FT(temperature), FT(salinity), extrapolate)
end

"""
    refractive_index(model, λ) -> FT

Return the real refractive index of seawater relative to air at wavelength `λ` `[nm]`
(a number or a vector).
"""
refractive_index(m::FixedRefractiveIndex, λ::Real) = m.n
refractive_index(m::AbstractRefractiveIndexModel, λ::AbstractVector) =
    map(Base.Fix1(refractive_index, m), λ)

# The message is one literal: this method also runs per element in a device
# broadcast, where a string cannot be built.
function refractive_index(m::QuanFry1995{FT}, λ::Real) where {FT}
    m.extrapolate || _within(λ, 400, 700) ||
        throw(DomainError(λ, "QuanFry1995 is fitted to 400–700 nm; pass extrapolate = true to use it outside, or use FixedRefractiveIndex"))
    return convert(_result_type(FT, λ),
                   _quan_fry_index(_wide(λ), _wide(m.temperature), _wide(m.salinity)))
end

# Quan & Fry (1995) Eq. (4) coefficients; relative to air, λ in nm.
function _quan_fry_index(λ, T, S)
    n0 = 1.31405;   n1 = 1.779e-4;  n2 = -1.05e-6;  n3 = 1.6e-8;    n4 = -2.02e-6
    n5 = 15.868;    n6 = 0.01155;   n7 = -0.00423;  n8 = -4382.0;   n9 = 1.1455e6
    return n0 + (n1 + n2 * T + n3 * T^2) * S + n4 * T^2 +
           (n5 + n6 * S + n7 * T) / λ + n8 / λ^2 + n9 / λ^3
end
