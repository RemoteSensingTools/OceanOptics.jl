# =============================================================================
# Coloured Dissolved Organic Matter (CDOM)
# =============================================================================
#
# CDOM is the optically-active dissolved organic pool in natural waters
# (humic/fulvic acids, by-products of phytoplankton decay). Its absorption
# rises sharply toward the UV and decays exponentially into the visible:
#
#     a_g(λ) = a_g(λ_ref) · exp(-slope · (λ - λ_ref))        [1/m]
#
# Slope ranges globally from ≈ 0.011 nm⁻¹ in coastal waters (Babin et al.
# 2003, JGR 108, 3211) up to ≈ 0.020 nm⁻¹ in pigment-rich tropical waters,
# with an accepted open-ocean mean of 0.014 nm⁻¹ at λ_ref = 440 nm per
# Bricaud, Morel & Prieur (1981), Limnol. Oceanogr. 26, 43.
#
# CDOM does not scatter light in the plane-parallel RT sense: the organic
# molecules are small and dilute enough that scattering contributes
# negligibly to the VSF. Subtyping `AbstractAbsorber` inherits the
# package-wide zero fallbacks for `scattering` and `backscattering`.
# =============================================================================

"""
$(TYPEDEF)

Bricaud, Morel & Prieur (1981) exponential-decay CDOM absorption,
`a_g(λ) = a_g(λ_ref) · exp(-slope · (λ - λ_ref))`. Currently the only
parameterization shipped; future additions (e.g. separate slopes for the
UV and visible hinge regions) can subtype `AbstractCDOMModel` without
altering `CDOM` itself.
"""
struct BricaudMorelPrieur1981 <: AbstractCDOMModel end

"""
$(TYPEDEF)

Coloured Dissolved Organic Matter (CDOM) as a pure absorber.

# Fields
$(TYPEDFIELDS)

# Examples

```julia
# Open-ocean default: a_g(440 nm) = 0.03 m⁻¹, slope 0.014 nm⁻¹
cdom = CDOM{Float64}()

# Coastal / humic-rich case
cdom_coastal = CDOM{Float64}(a_ref = 0.5, λ_ref = 440.0, slope = 0.017)

# Evaluate
absorption(cdom, 440.0)   # ≈ 0.03
absorption(cdom, 500.0)   # exponentially smaller
```
"""
struct CDOM{FT, M<:AbstractCDOMModel} <: AbstractAbsorber{FT}
    "Parameterization-model tag (e.g. `BricaudMorelPrieur1981()`)"
    model::M
    "Reference absorption `a_g(λ_ref)` `[1/m]`"
    a_ref::FT
    "Reference wavelength `[nm]` at which `a_ref` is anchored"
    λ_ref::FT
    "Spectral slope `S` `[1/nm]`"
    slope::FT
end

"""
    CDOM{FT}(; model=BricaudMorelPrieur1981(), a_ref=0.03, λ_ref=440, slope=0.014)

Construct a CDOM absorber. Default values reproduce the Bricaud, Morel &
Prieur (1981) open-ocean average. `slope` outside the Babin et al. (2003)
global survey range of 0.005–0.030 nm⁻¹ triggers a non-fatal warning.
"""
function CDOM{FT}(;
        model::AbstractCDOMModel = BricaudMorelPrieur1981(),
        a_ref                    = FT(0.03),
        λ_ref                    = FT(440),
        slope                    = FT(0.014),
) where {FT<:AbstractFloat}
    a_ref ≥ zero(a_ref) ||
        throw(ArgumentError("CDOM: a_ref must be non-negative, got $a_ref"))
    λ_ref > zero(λ_ref) ||
        throw(ArgumentError("CDOM: λ_ref must be positive, got $λ_ref"))
    (FT(0.005) ≤ slope ≤ FT(0.030)) ||
        @warn "CDOM: slope $slope [1/nm] outside typical 0.005-0.030 range"
    return CDOM{FT, typeof(model)}(model, FT(a_ref), FT(λ_ref), FT(slope))
end

CDOM(; kwargs...) = CDOM{Float64}(; kwargs...)

# Absorption — analytic, branch-free, AD-transparent.
#
# The returned type is `promote_type(FT, typeof(λ))`: passing a plain
# `Float64` gives `FT`; passing a `ForwardDiff.Dual` gives `Dual` so the
# derivative w.r.t. λ is tracked through `exp`.
absorption(c::CDOM, λ::Real)           = c.a_ref * exp(-c.slope * (λ - c.λ_ref))
absorption(c::CDOM, λ::AbstractVector) = [absorption(c, λi) for λi in λ]

function Base.show(io::IO, c::CDOM{FT, M}) where {FT, M}
    print(io, "CDOM{$FT}(", nameof(M), ", a($(c.λ_ref))=$(c.a_ref), S=$(c.slope))")
end
