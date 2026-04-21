"""
$(TYPEDEF)

Complete bundle of optical properties for one ocean layer, ready to hand to
a radiative-transfer solver. This is the *output* of [`layer_optics`](@ref);
downstream code (vSmartMOM, HydroLight wrappers, Monte-Carlo) consumes it.

# Fields
$(TYPEDFIELDS)

# Relationship to `IOP`

An `IOP` carries the volumetric coefficients `(a, b, bb)`; an
`OceanLayerOptics` carries what a layer-by-layer RT solver actually wants:

* `τ = (a + b) · Δz`  — layer optical depth (extinction × thickness)
* `ϖ = b / (a + b)`   — single-scattering albedo
* `β_ℓ`              — Legendre-moment expansion of the (mixed) phase function

For polarized RT all six Greek coefficients `(α, β, γ, δ, ε, ζ)` are carried;
for scalar RT only `β` is meaningful and the others may be empty / absent.

# Fourier / Greek coefficients

Arrays are indexed `[ℓ, λ]` where `ℓ = 0:ℓ_max` and `λ` runs over the
wavelength grid. This layout matches vSmartMOM's `GreekCoefs` and, by
contracts elsewhere in the package, is what `to_core_scattering_properties`
hands off.
"""
struct OceanLayerOptics{FT, V<:AbstractVector{FT}, M<:AbstractMatrix{FT}}
    "Wavelength grid `[nm]`"
    λ::V
    "Layer thickness `[m]`"
    Δz::FT
    "Layer optical depth per wavelength, `τ(λ) = c(λ)·Δz`"
    τ::V
    "Single-scattering albedo per wavelength, `ϖ(λ)`"
    ϖ::V
    "Backscatter fraction per wavelength, `B(λ) = b_b(λ)/b(λ)`; diagnostic only"
    B::V

    # Greek coefficients β_ℓ(λ) — `[ℓ_max+1, nλ]`
    "Phase-function expansion coefficients β_ℓ(λ) (index ℓ = 0…ℓ_max; row 1 is ℓ=0)"
    β::M

    # Polarized Stokes elements — empty matrices for scalar-only layers
    "α_ℓ(λ), B[2,2] of the polarized scattering matrix (empty → scalar RT)"
    α::M
    "γ_ℓ(λ), B[1,2] = B[2,1]"
    γ::M
    "δ_ℓ(λ), B[4,4]"
    δ::M
    "ε_ℓ(λ), B[3,4] = -B[4,3]"
    ε::M
    "ζ_ℓ(λ), B[3,3]"
    ζ::M

    # Optional: source terms (fluorescence). `nothing` if no emitters.
    "Isotropic (Fourier m=0 only) fluorescence source `[W/m²/sr/nm]` per λ; `nothing` if absent"
    J_source::Union{Nothing, V}
end

# Convenience constructor: scalar-only layer (no polarization, no fluorescence).
# Julia's `where` clause only binds type variables from *positional* argument
# types, not from keyword-argument types — so we derive `FT` by hand from
# `eltype(typeof(λ))` rather than declaring it in the where. This also lets
# the caller pass `Δz` as any `Real` (it is converted on assignment).
function OceanLayerOptics(; λ::AbstractVector, Δz::Real,
                           τ::AbstractVector, ϖ::AbstractVector, B::AbstractVector,
                           β::AbstractMatrix,
                           α::AbstractMatrix = similar(β, 0, 0),
                           γ::AbstractMatrix = similar(β, 0, 0),
                           δ::AbstractMatrix = similar(β, 0, 0),
                           ε::AbstractMatrix = similar(β, 0, 0),
                           ζ::AbstractMatrix = similar(β, 0, 0),
                           J_source = nothing)
    length(λ) == length(τ) == length(ϖ) == length(B) ||
        throw(DimensionMismatch("OceanLayerOptics: λ/τ/ϖ/B length mismatch"))
    size(β, 2) == length(λ) ||
        throw(DimensionMismatch("OceanLayerOptics: β must have one column per λ"))
    FT = eltype(λ)
    V  = typeof(λ)
    M  = typeof(β)
    return OceanLayerOptics{FT, V, M}(λ, FT(Δz), τ, ϖ, B, β, α, γ, δ, ε, ζ, J_source)
end

"Number of Legendre moments stored (ℓ_max + 1)."
n_moments(o::OceanLayerOptics) = size(o.β, 1)

"Number of wavelengths."
n_wavelengths(o::OceanLayerOptics) = length(o.λ)

"Whether this layer carries polarized Greek coefficients."
is_polarized(o::OceanLayerOptics) = !isempty(o.α)

"Whether this layer carries a fluorescence source term."
has_source(o::OceanLayerOptics) = o.J_source !== nothing

function Base.show(io::IO, o::OceanLayerOptics{FT}) where {FT}
    println(io, "OceanLayerOptics{$FT}(Δz=$(o.Δz) m, $(n_wavelengths(o)) λ, ",
                "$(n_moments(o)-1) moments, ",
                is_polarized(o) ? "polarized" : "scalar",
                has_source(o) ? ", +fluorescence" : "", ")")
end
