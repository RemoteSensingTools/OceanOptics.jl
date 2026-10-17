"""
$(TYPEDEF)

Complete bundle of optical properties for one ocean layer, ready to hand to
a radiative-transfer solver. This is the *output* of [`layer_optics`](@ref);
downstream code (vSmartMOM, HydroLight wrappers, Monte-Carlo) consumes it.
The full consumer contract (units, conventions, zero-scattering behavior)
is documented on the *Consumer API* page.

# Fields
$(TYPEDFIELDS)

# Relationship to `IOP`

An `IOP` carries the volumetric coefficients `(a, b, bb)`; an
`OceanLayerOptics` carries what a layer-by-layer RT solver actually wants:

* `τ = (a + b) · Δz`  — layer optical depth (extinction × thickness)
* `ϖ = b / (a + b)`   — single-scattering albedo
* `β_ℓ`              — Legendre-moment expansion of the (mixed) phase function

For polarized RT all six Greek coefficients `(α, β, γ, δ, ϵ, ζ)` are carried;
for scalar RT only `β` is meaningful and the other five are empty `0×0`
matrices. A bundle is either scalar (all five empty) or polarized (all five
the same size as `β`); partial bundles are rejected.

# Fourier / Greek coefficients

Arrays are indexed `[ℓ+1, iλ]` where `ℓ = 0:ℓ_max` and `iλ` runs over the
wavelength grid. Coefficients follow `p(μ) = Σ_ℓ β_ℓ P_ℓ(μ)` with the
`(2ℓ+1)` factor already included, so `β₀ = 1`. A column `[:, iλ]` matches
the per-wavelength vector layout of vSmartMOM's `GreekCoefs`.

# Validation

Every constructor checks that `λ` is finite, positive and strictly
monotone; `Δz` is finite and positive; `τ ≥ 0`, `0 ≤ ϖ ≤ 1`, `0 ≤ B ≤ 1`;
all coefficients are finite; `β₀ = 1` to within `√eps(FT)`; and
`|β_ℓ| ≤ 2ℓ+1`, which every non-negative phase function satisfies. For a
polarized bundle, `α, γ, ϵ, ζ` must vanish at `ℓ < 2`, where their
generalized spherical functions are identically zero.

These are necessary conditions, not proof of a correct import: the
moment bound catches a doubled `(2ℓ+1)` factor only for strongly peaked
phase functions and never catches a missing one, so also compare
`β₁/3` with the published asymmetry parameter. Arrays are checked when
the bundle is constructed; mutating them afterwards is not re-checked.
"""
struct OceanLayerOptics{FT, V<:AbstractVector{FT}, M<:AbstractMatrix{FT}}
    "Wavelength grid `[nm]`, strictly monotone"
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
    "ϵ_ℓ(λ), B[3,4] = -B[4,3]"
    ϵ::M
    "ζ_ℓ(λ), B[3,3]"
    ζ::M

    # Optional: source terms (fluorescence). `nothing` if no emitters.
    "Isotropic (Fourier m=0 only) fluorescence source `[W/m²/sr/nm]` per λ; `nothing` if absent"
    J_source::Union{Nothing, V}

    function OceanLayerOptics{FT, V, M}(λ, Δz, τ, ϖ, B, β, α, γ, δ, ϵ, ζ,
                                        J_source) where {FT, V<:AbstractVector{FT},
                                                         M<:AbstractMatrix{FT}}
        o = new{FT, V, M}(λ, Δz, τ, ϖ, B, β, α, γ, δ, ϵ, ζ, J_source)
        _validate_layer_optics(o)
        return o
    end
end

# Positional form, as provided implicitly before the validating inner
# constructor existed.
OceanLayerOptics(λ::V, Δz::Real, τ::V, ϖ::V, B::V, β::M, α::M, γ::M, δ::M, ϵ::M, ζ::M,
                 J_source::Union{Nothing, V}) where {FT, V<:AbstractVector{FT},
                                                     M<:AbstractMatrix{FT}} =
    OceanLayerOptics{FT, V, M}(λ, FT(Δz), τ, ϖ, B, β, α, γ, δ, ϵ, ζ, J_source)

# Convenience constructor: scalar-only layer unless all five polarized
# Greek matrices are given. Julia's `where` clause only binds type
# variables from *positional* argument types, not from keyword-argument
# types — so we derive `FT` by hand from `eltype(λ)`. This also lets the
# caller pass `Δz` as any `Real` (it is converted on assignment).
function OceanLayerOptics(; λ::AbstractVector, Δz::Real,
                           τ::AbstractVector, ϖ::AbstractVector, B::AbstractVector,
                           β::AbstractMatrix,
                           α::AbstractMatrix = similar(β, 0, 0),
                           γ::AbstractMatrix = similar(β, 0, 0),
                           δ::AbstractMatrix = similar(β, 0, 0),
                           ϵ::AbstractMatrix = similar(β, 0, 0),
                           ζ::AbstractMatrix = similar(β, 0, 0),
                           J_source = nothing)
    FT = eltype(λ)
    V  = typeof(λ)
    M  = typeof(β)
    for (name, v) in ((:τ, τ), (:ϖ, ϖ), (:B, B))
        v isa V || throw(ArgumentError("OceanLayerOptics: $name must have the same " *
                                       "array type as λ ($V), got $(typeof(v))"))
    end
    for (name, m) in ((:α, α), (:γ, γ), (:δ, δ), (:ϵ, ϵ), (:ζ, ζ))
        m isa M || throw(ArgumentError("OceanLayerOptics: $name must have the same " *
                                       "array type as β ($M), got $(typeof(m))"))
    end
    J_source === nothing || J_source isa V ||
        throw(ArgumentError("OceanLayerOptics: J_source must be `nothing` or have " *
                            "the same array type as λ ($V), got $(typeof(J_source))"))
    return OceanLayerOptics{FT, V, M}(λ, FT(Δz), τ, ϖ, B, β, α, γ, δ, ϵ, ζ, J_source)
end

function _validate_layer_optics(o::OceanLayerOptics{FT}) where {FT}
    (; λ, Δz, τ, ϖ, B, β, J_source) = o
    nλ = length(λ)
    nλ ≥ 1 || throw(ArgumentError("OceanLayerOptics: empty wavelength grid"))
    length(τ) == length(ϖ) == length(B) == nλ ||
        throw(DimensionMismatch("OceanLayerOptics: λ/τ/ϖ/B length mismatch " *
                                "($(nλ)/$(length(τ))/$(length(ϖ))/$(length(B)))"))
    size(β, 1) ≥ 1 && size(β, 2) == nλ ||
        throw(DimensionMismatch("OceanLayerOptics: β must be (ℓ_max+1, nλ) = (≥1, $nλ), " *
                                "got $(size(β))"))

    all(x -> isfinite(x) && x > zero(FT), λ) ||
        throw(ArgumentError("OceanLayerOptics: wavelengths must be finite and positive"))
    _is_strictly_monotone(λ) ||
        throw(ArgumentError("OceanLayerOptics: wavelengths must be strictly increasing " *
                            "or strictly decreasing"))
    isfinite(Δz) && Δz > zero(FT) ||
        throw(ArgumentError("OceanLayerOptics: Δz must be finite and positive, got $Δz"))
    all(x -> isfinite(x) && x ≥ zero(FT), τ) ||
        throw(ArgumentError("OceanLayerOptics: τ must be finite and non-negative"))
    all(x -> isfinite(x) && zero(FT) ≤ x ≤ one(FT), ϖ) ||
        throw(ArgumentError("OceanLayerOptics: ϖ must lie in [0, 1]"))
    all(x -> isfinite(x) && zero(FT) ≤ x ≤ one(FT), B) ||
        throw(ArgumentError("OceanLayerOptics: B must lie in [0, 1]"))

    _validate_beta_moments(β, "OceanLayerOptics")

    greek = (α = o.α, γ = o.γ, δ = o.δ, ϵ = o.ϵ, ζ = o.ζ)
    empty_names = [name for (name, m) in pairs(greek) if isempty(m)]
    if length(empty_names) != length(greek)
        isempty(empty_names) ||
            throw(ArgumentError("OceanLayerOptics: polarized bundles need all five of " *
                                "α, γ, δ, ϵ, ζ; empty: " * join(empty_names, ", ")))
        for (name, m) in pairs(greek)
            size(m) == size(β) ||
                throw(DimensionMismatch("OceanLayerOptics: $name has size $(size(m)), " *
                                        "β has size $(size(β))"))
            all(isfinite, m) ||
                throw(ArgumentError("OceanLayerOptics: $name contains non-finite values"))
        end
        _validate_low_order_greek(greek, "OceanLayerOptics")
    end

    if J_source !== nothing
        length(J_source) == nλ ||
            throw(DimensionMismatch("OceanLayerOptics: J_source must have one entry per λ"))
        all(x -> isfinite(x) && x ≥ zero(FT), J_source) ||
            throw(ArgumentError("OceanLayerOptics: J_source must be finite and non-negative"))
    end
    return nothing
end

_is_strictly_increasing(x::AbstractVector) = all(i -> x[i] < x[i+1], firstindex(x):lastindex(x)-1)
_is_strictly_monotone(x::AbstractVector)   =
    _is_strictly_increasing(x) || all(i -> x[i] > x[i+1], firstindex(x):lastindex(x)-1)

# Checks shared by every β-moment container ([ℓ+1, column] layout):
# finite values, β₀ = 1 and the non-negative-phase-function bound
# |β_ℓ| ≤ 2ℓ+1 (|P_ℓ| ≤ 1 on [-1, 1]). Whole-array reductions; the
# offending entry is located only for the error message.
function _validate_beta_moments(β::AbstractMatrix{FT}, who::AbstractString) where {FT}
    all(isfinite, β) ||
        throw(ArgumentError("$who: β contains non-finite values"))
    tol = sqrt(eps(FT))
    β₀  = view(β, 1:1, :)
    if !all(x -> abs(x - one(FT)) ≤ tol, β₀)
        j = findfirst(x -> abs(x - one(FT)) > tol, vec(collect(β₀)))
        throw(ArgumentError("$who: β₀ = $(β₀[j]) in column $j; phase functions must be " *
                            "normalized so that β₀ = 1 (tolerance $tol)"))
    end
    bound = (2 .* (0:size(β, 1)-1) .+ 1) .* (one(FT) + tol)
    if !all(abs.(β) .≤ bound)
        I = findfirst(abs.(collect(β)) .> collect(bound))
        ℓ = I[1] - 1
        throw(ArgumentError("$who: |β_$ℓ| = $(abs(β[I])) > 2ℓ+1 = $(2ℓ+1) in column " *
                            "$(I[2]); no non-negative phase function allows this (was " *
                            "the (2ℓ+1) factor applied twice?)"))
    end
    return nothing
end

# α, γ, ϵ, ζ multiply generalized spherical functions with |m| or |n| = 2,
# which vanish identically for ℓ < 2; non-zero entries there signal a
# row-shifted or mis-indexed table.
function _validate_low_order_greek(greek::NamedTuple, who::AbstractString)
    for name in (:α, :γ, :ϵ, :ζ)
        m   = getfield(greek, name)
        tol = sqrt(eps(eltype(m)))
        all(x -> abs(x) ≤ tol, view(m, 1:min(2, size(m, 1)), :)) ||
            throw(ArgumentError("$who: $(name)_0 and $(name)_1 must be zero (their " *
                                "generalized spherical functions vanish for ℓ < 2); " *
                                "is the table shifted by a row?"))
    end
    return nothing
end

"Number of Legendre moments stored (ℓ_max + 1)."
n_moments(o::OceanLayerOptics) = size(o.β, 1)

"Number of wavelengths."
n_wavelengths(o::OceanLayerOptics) = length(o.λ)

"""
    is_polarized(o::OceanLayerOptics) -> Bool

Whether this layer carries all six polarized Greek coefficients. The
constructor guarantees that α, γ, δ, ϵ, ζ are either all empty or all
populated.
"""
is_polarized(o::OceanLayerOptics) = !isempty(o.α)

"Whether this layer carries a fluorescence source term."
has_source(o::OceanLayerOptics) = o.J_source !== nothing

function Base.show(io::IO, o::OceanLayerOptics{FT}) where {FT}
    println(io, "OceanLayerOptics{$FT}(Δz=$(o.Δz) m, $(n_wavelengths(o)) λ, ",
                "$(n_moments(o)-1) moments, ",
                is_polarized(o) ? "polarized" : "scalar",
                has_source(o) ? ", +fluorescence" : "", ")")
end
