# =============================================================================
# Tabulated Greek-coefficient phase function
# =============================================================================
#
# Phase function or phase matrix given directly by its expansion — from a
# published table, a Mie code or a fit — in this package's convention (see
# `GreekCoefficients`). Coefficients in other conventions must be mapped
# first: VLIDORT/Mishchenko-format β₁ and β₂ become γ = −β₁ and ϵ = −β₂, the
# flips vSmartMOM applies to such data. Angular scattering matrices belong
# in `TabulatedMuellerPhase`; per-wavelength tables in
# `WavelengthDependentPhase`.
# =============================================================================

"""
$(TYPEDEF)

Phase function or phase matrix given by its Greek expansion.

# Fields
$(TYPEDFIELDS)

# Construction

    TabulatedGreekPhase(g::GreekCoefficients; source = "")
    TabulatedGreekPhase(β; α, γ, δ, ϵ, ζ, source = "", normalize = false)

The keyword form builds the [`GreekCoefficients`](@ref) (pass none or all
of `α, γ, δ, ϵ, ζ`). `normalize = true` divides all six coefficients by
`β₀`, an explicit choice for published tables rounded to a few digits.
The expansion must also give a backscatter fraction in `[0, 1]`, which a
strongly peaked table truncated far before convergence can fail.

Requests for `ℓ_max` below the table length truncate it; requests above
pad with zeros, which is exact for the tabulated expansion.

# Example

```julia
# Rayleigh scattering without depolarization, as a phase matrix
pf = TabulatedGreekPhase([1.0, 0.0, 0.5]; α = [0.0, 0.0, 3.0],
                         γ = [0.0, 0.0, sqrt(1.5)], δ = [0.0, 1.5, 0.0],
                         ϵ = zeros(3), ζ = zeros(3), source = "Hansen & Travis (1974)")
greek_expansion(pf, 8; polarized = true).γ[3]    # √(3/2)
```
"""
struct TabulatedGreekPhase{FT} <: AbstractOceanPhaseFunction{FT}
    "The expansion"
    greek::GreekCoefficients{FT}
    "Provenance: citation, table, and any transformation applied"
    source::String

    function TabulatedGreekPhase{FT}(g::GreekCoefficients{FT}, source::String) where {FT<:AbstractFloat}
        B = backscatter_fraction(g)
        0 ≤ B ≤ 1 ||
            throw(ArgumentError("TabulatedGreekPhase: the expansion has backscatter fraction " *
                                "$B ∉ [0, 1] and is not a valid phase function (table " *
                                "truncated before convergence?)"))
        new{FT}(g, source)
    end
end

TabulatedGreekPhase(g::GreekCoefficients{FT}; source::AbstractString = "") where {FT} =
    TabulatedGreekPhase{FT}(g, String(source))

function TabulatedGreekPhase(β::AbstractVector; α = nothing, γ = nothing, δ = nothing,
                             ϵ = nothing, ζ = nothing, source::AbstractString = "",
                             normalize::Bool = false)
    if normalize
        β₀ = first(β)
        isfinite(β₀) && β₀ > 0 ||
            throw(ArgumentError("TabulatedGreekPhase: normalize = true needs a positive β₀"))
        sc = x -> x === nothing ? nothing : x ./ β₀
        β, α, γ, δ, ϵ, ζ = β ./ β₀, sc(α), sc(γ), sc(δ), sc(ϵ), sc(ζ)
    end
    return TabulatedGreekPhase(GreekCoefficients(β; α, γ, δ, ϵ, ζ); source)
end

is_polarizable(pf::TabulatedGreekPhase) = is_polarized(pf.greek)

phase_function_moments(pf::TabulatedGreekPhase, ℓ_max::Integer) =
    resize_moments(pf.greek, ℓ_max).β

_polarized_greek(pf::TabulatedGreekPhase, ℓ_max::Integer) = resize_moments(pf.greek, ℓ_max)

phase_function_value(pf::TabulatedGreekPhase, μ::Real) =
    _reconstruct_from_moments(pf.greek.β, μ)

backscatter_fraction(pf::TabulatedGreekPhase) = backscatter_fraction(pf.greek)
asymmetry_parameter(pf::TabulatedGreekPhase)  = asymmetry_parameter(pf.greek)

function Base.show(io::IO, pf::TabulatedGreekPhase{FT}) where {FT}
    print(io, "TabulatedGreekPhase{$FT}(ℓ ≤ ", n_moments(pf.greek) - 1, ", ",
          is_polarizable(pf) ? "phase matrix" : "phase function",
          isempty(pf.source) ? "" : ", \"$(pf.source)\"", ")")
end
