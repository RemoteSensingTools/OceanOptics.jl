# =============================================================================
# Wavelength-dependent phase functions
# =============================================================================
#
# A different phase function (or phase matrix) at each tabulated
# wavelength — e.g. Mie hydrosols, or the two Chowdhary et al. (2020) AOS-IV
# matrices that apply to different wavelengths. There is no interpolation
# between phase functions: lookup is exact, or — chosen explicitly — the
# nearest tabulated wavelength, so that a fine spectral grid within a band
# uses one band-representative phase matrix.
# =============================================================================

"""
$(TYPEDEF)

One phase function per tabulated wavelength.

# Fields
$(TYPEDFIELDS)

# Construction

    WavelengthDependentPhase(λ, phases; lookup = :exact)

`λ` `[nm]` strictly increasing, `phases` wavelength-independent phase
functions of one float type. With `lookup = :exact` every query must name
a tabulated wavelength (`phase_function_moments(pf, ℓ_max, λ)` etc.) and
other wavelengths throw; with `lookup = :nearest` the phase function of
the nearest tabulated wavelength is used. The phase function is
polarizable only if all entries are.
"""
struct WavelengthDependentPhase{FT, P<:AbstractOceanPhaseFunction{FT}} <: AbstractOceanPhaseFunction{FT}
    "Wavelengths `[nm]`, strictly increasing"
    λ::Vector{FT}
    "Phase function at each wavelength"
    phases::Vector{P}
    "`:exact` or `:nearest`"
    lookup::Symbol

    function WavelengthDependentPhase{FT, P}(λ, phases, lookup) where {FT, P}
        lookup in (:exact, :nearest) ||
            throw(ArgumentError("WavelengthDependentPhase: lookup must be :exact or :nearest"))
        _check_exact_λ_grid(λ, "WavelengthDependentPhase")
        length(phases) == length(λ) ||
            throw(DimensionMismatch("WavelengthDependentPhase: $(length(phases)) phase " *
                                    "functions for $(length(λ)) wavelengths"))
        any(is_wavelength_dependent, phases) &&
            throw(ArgumentError("WavelengthDependentPhase: entries must be wavelength independent"))
        new{FT, P}(λ, phases, lookup)
    end
end

function WavelengthDependentPhase(λ::AbstractVector, phases::AbstractVector{<:AbstractOceanPhaseFunction{FT}};
                                  lookup::Symbol = :exact) where {FT}
    P = typeintersect(eltype(phases), AbstractOceanPhaseFunction{FT})
    return WavelengthDependentPhase{FT, P}(collect(FT, λ), collect(P, phases), lookup)
end

is_wavelength_dependent(::WavelengthDependentPhase) = true
is_polarizable(pf::WavelengthDependentPhase) = all(is_polarizable, pf.phases)

function _phase_entry(pf::WavelengthDependentPhase, λ::Real)
    pf.lookup === :exact && return _exact_λ_index(pf.λ, λ, "WavelengthDependentPhase")
    i = clamp(searchsortedfirst(pf.λ, λ), 1, length(pf.λ))
    return (i > 1 && abs(pf.λ[i-1] - λ) ≤ abs(pf.λ[i] - λ)) ? i - 1 : i
end

_phase_at(pf::WavelengthDependentPhase, λ::Real) = pf.phases[_phase_entry(pf, λ)]
_phase_key(pf::WavelengthDependentPhase, λ::Real) = _phase_entry(pf, λ)

phase_function_moments(pf::WavelengthDependentPhase, ℓ_max::Integer, λ::Real) =
    phase_function_moments(_phase_at(pf, λ), ℓ_max)
greek_expansion(pf::WavelengthDependentPhase, ℓ_max::Integer, λ::Real; polarized::Bool = false) =
    greek_expansion(_phase_at(pf, λ), ℓ_max; polarized)
backscatter_fraction(pf::WavelengthDependentPhase, λ::Real) = backscatter_fraction(_phase_at(pf, λ))
asymmetry_parameter(pf::WavelengthDependentPhase, λ::Real)  = asymmetry_parameter(_phase_at(pf, λ))

for f in (:phase_function_moments, :greek_expansion)
    @eval $f(pf::WavelengthDependentPhase, ℓ_max::Integer; kwargs...) = _requires_wavelength(pf)
end
backscatter_fraction(pf::WavelengthDependentPhase) = _requires_wavelength(pf)
asymmetry_parameter(pf::WavelengthDependentPhase)  = _requires_wavelength(pf)
phase_function_value(pf::WavelengthDependentPhase, μ::Real) = _requires_wavelength(pf)

function Base.show(io::IO, pf::WavelengthDependentPhase{FT}) where {FT}
    print(io, "WavelengthDependentPhase{$FT}(", length(pf.λ), " λ, ",
          is_polarizable(pf) ? "phase matrices" : "phase functions", ")")
end
