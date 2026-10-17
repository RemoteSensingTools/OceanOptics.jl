# =============================================================================
# Layer assembly — take an OceanLayer and produce an OceanLayerOptics
# =============================================================================
#
# This is the main user-facing function of the elastic path. Given
#
#   • an OceanLayer (constituents + thickness)
#   • a wavelength grid λ_grid [nm]
#   • a truncation ℓ_max for the Legendre expansion
#
# produce a self-contained OceanLayerOptics bundle ready to feed into any
# plane-parallel RT solver using Matrix Operator / Adding-Doubling / Discrete
# Ordinates / etc.
#
# Algorithm
# ---------
# For each wavelength λ:
#   1. For each constituent c, compute (a_c, b_c, bb_c)
#   2. For each scattering constituent, obtain its phase function moments β_c
#   3. Sum extensive quantities:  a = Σ a_c,  b = Σ b_c,  bb = Σ bb_c
#   4. Mix phase-function moments scattering-weighted: β = Σ (b_c / b) β_c
#   5. Compute τ = (a + b) Δz,  ϖ = b / (a + b),  B = bb / b
#
# The whole operation is linear in the constituent count and wavelength count.
# =============================================================================

"""
    layer_optics(layer, λ_grid; ℓ_max=64, polarized=false,
                 scalar_phase_approximation=nothing, FT=Float64) -> OceanLayerOptics

Compute the complete optical bundle for one ocean layer at each wavelength
in `λ_grid`.

# Arguments
- `layer::OceanLayer` — the layer (constituents + thickness)
- `λ_grid` — strictly monotone vector of wavelengths `[nm]`
- `ℓ_max` — maximum Legendre order for the phase-function expansion
  (default 64, which is adequate for Gauss-Lobatto RT with N ≤ 32)
- `polarized` — also compute the polarized Greek coefficients
  `(α, γ, δ, ϵ, ζ)` via [`phase_matrix_moments`](@ref); populates the
  corresponding `OceanLayerOptics` fields. Default `false` for scalar RT.
- `scalar_phase_approximation` — what to do, when `polarized = true`, with
  scatterers whose phase function is scalar-only (`is_polarizable(pf) ==
  false`, e.g. Fournier-Forand particles). `nothing` (default) throws: a
  scalar phase function does not determine a phase matrix.
  `:depolarizing` explicitly selects the completely depolarizing
  approximation of [`depolarizing_phase_matrix_moments`](@ref) for those
  scatterers. The check applies to every scatterer in the layer, even one
  whose scattering coefficient happens to be zero. Ignored when
  `polarized = false`.
- `FT` — floating-point precision

# Returns
An [`OceanLayerOptics`](@ref) containing τ, ϖ, β_ℓ(λ) and the diagnostic
backscatter fraction B, all on the `λ_grid`. When `polarized = true` the
α/γ/δ/ϵ/ζ matrices are filled in the same `(ℓ_max+1, nλ)` layout.

# Mixing rules and checks
Absorption and scattering coefficients add; Greek coefficients mix with
weights `b_c / b`. A layer with no scattering at some wavelength gets
`ϖ = 0`, `B = 0`, `β = [1, 0, …]` and zero polarized coefficients there.
Phase functions with `is_wavelength_dependent(pf)` are evaluated at every
wavelength; all others once per layer.

If a constituent declares [`is_total_absorption`](@ref) (a
`TabulatedAbsorber` with `total = true`), it must be the only one, and
every other constituent must have zero absorption at every wavelength —
otherwise its absorbers would be counted twice.

# Notes
This function handles only the *elastic* IOPs. Inelastic source terms
(fluorescence, Raman) are computed separately by a downstream two-pass
RT solver from a scalar irradiance field `E°(τ, λ)` obtained from the
first pass (see Fell 1997 §2.2.9).

# Example

```julia
water = PureWater{Float64}()
chl   = Phytoplankton{Float64, Bricaud1998}(Chl = 0.5)
layer = OceanLayer(0.0, 5.0, [water, chl])
λ     = 400.0:5.0:700.0

opt = layer_optics(layer, collect(λ); ℓ_max = 64)
# opt.τ, opt.ϖ, opt.β all have one column/entry per λ.

# Polarized (Stokes-vector) RT. Fournier-Forand particles carry no phase
# matrix, so their treatment must be chosen explicitly:
opt_p = layer_optics(layer, collect(λ); ℓ_max = 64, polarized = true,
                     scalar_phase_approximation = :depolarizing)
# opt_p.α, opt_p.γ, opt_p.δ, opt_p.ϵ, opt_p.ζ are populated too.
```
"""
function layer_optics(layer::OceanLayer, λ_grid::AbstractVector{<:Real};
                      ℓ_max::Integer = 64,
                      polarized::Bool = false,
                      scalar_phase_approximation::Union{Nothing, Symbol} = nothing,
                      FT::Type = Float64)
    ℓ_max ≥ 0 || throw(ArgumentError("ℓ_max must be non-negative, got $ℓ_max"))
    scalar_phase_approximation in (nothing, :depolarizing) ||
        throw(ArgumentError("scalar_phase_approximation must be `nothing` or " *
                            ":depolarizing, got :$scalar_phase_approximation"))
    nλ = length(λ_grid)
    Δz = FT(thickness(layer))
    λ  = FT.(λ_grid)
    N  = ℓ_max + 1

    # Allocate outputs
    τ = zeros(FT, nλ)
    ϖ = zeros(FT, nλ)
    B = zeros(FT, nλ)
    β = zeros(FT, N, nλ)
    if polarized
        α = zeros(FT, N, nλ)
        γ = zeros(FT, N, nλ)
        δ = zeros(FT, N, nλ)
        ϵ = zeros(FT, N, nλ)
        ζ = zeros(FT, N, nλ)
    end

    constituents = layer.constituents
    i_total = findall(is_total_absorption, constituents)
    length(i_total) ≤ 1 ||
        throw(ArgumentError("layer_optics: $(length(i_total)) constituents declare total " *
                            "(bulk) absorption; at most one can"))
    i_bulk = isempty(i_total) ? 0 : only(i_total)

    # Split constituents by role — avoids re-dispatching per λ inside the loop
    # and lets us extract phase functions once.
    scatterers = filter(c -> c isa Union{AbstractScatterer, AbstractAbsorbingScatterer},
                        constituents)
    phase_fns  = [phase_function(c) for c in scatterers]
    if polarized && scalar_phase_approximation === nothing
        scalar_only = [nameof(typeof(c)) => nameof(typeof(pf))
                       for (c, pf) in zip(scatterers, phase_fns) if !is_polarizable(pf)]
        isempty(scalar_only) ||
            throw(ArgumentError("layer_optics(polarized = true): scalar-only phase " *
                                "functions for " * join(string.(scalar_only), ", ") *
                                ". Their phase matrices are unknown; supply a full phase " *
                                "matrix (e.g. TabulatedScatterer with a polarized " *
                                "TabulatedGreekPhase) or pass " *
                                "`scalar_phase_approximation = :depolarizing` to select " *
                                "the completely depolarizing approximation explicitly."))
    end

    # Greek coefficients per scatterer: computed once for wavelength-
    # independent phase functions (most ocean cases), per λ otherwise.
    # (λ-independent phase functions ignore the wavelength argument.)
    coefficients(pf, λ_i) = polarized ?
        (is_polarizable(pf) ? phase_matrix_moments(pf, ℓ_max, λ_i) :
                              depolarizing_phase_matrix_moments(pf, ℓ_max, λ_i)) :
        (β = phase_function_moments(pf, ℓ_max, λ_i),)
    fixed_coefficients = [is_wavelength_dependent(pf) || isempty(λ) ? nothing :
                          coefficients(pf, first(λ)) for pf in phase_fns]

    b_per_scatterer = zeros(FT, length(scatterers))
    @inbounds for (iλ, λ_i) in enumerate(λ)
        a_total  = zero(FT)
        b_total  = zero(FT)
        bb_total = zero(FT)

        for (i, c) in enumerate(constituents)
            a_c = FT(absorption(c, λ_i))
            if i_bulk != 0 && i != i_bulk && a_c != zero(FT)
                throw(ArgumentError("layer_optics: $(constituents[i_bulk]) is the " *
                                    "total absorption of the layer, but " *
                                    "$(nameof(typeof(c))) also absorbs (a = $a_c 1/m at " *
                                    "λ = $λ_i nm); this would count its absorbers twice"))
            end
            a_total += a_c
        end
        for (j, c) in enumerate(scatterers)
            b_j = FT(scattering(c, λ_i))
            b_per_scatterer[j] = b_j
            b_total  += b_j
            bb_total += FT(backscattering(c, λ_i))
        end

        c_total = a_total + b_total
        τ[iλ] = c_total * Δz
        ϖ[iλ] = c_total > zero(FT) ? b_total / c_total : zero(FT)
        B[iλ] = b_total > zero(FT) ? bb_total / b_total : zero(FT)

        coefs = [fixed === nothing ? coefficients(pf, λ_i) : fixed
                 for (fixed, pf) in zip(fixed_coefficients, phase_fns)]

        # Scattering-weighted phase-moment mixing. Pass the per-layer
        # moment count explicitly so empty-scatterer layers return a
        # correctly-sized (β₀ = 1, β_{ℓ≥1} = 0) vector rather than the
        # length-1 fallback, which would broadcast-fill all rows with 1.
        β[:, iλ] .= mix_phase_moments(b_per_scatterer, Vector{FT}[cf.β for cf in coefs];
                                      n_moments = N)

        if polarized
            # Mix each Greek coefficient with the same scattering weights.
            # The α/γ/δ/ϵ/ζ scalars (not β) have no `[1] = 1` normalization
            # invariant, so we use a direct weighted sum without the
            # fallback-to-unit-isotropic behavior of mix_phase_moments.
            α[:, iλ] .= _greek_mix(b_per_scatterer, coefs, :α, b_total, FT, N)
            γ[:, iλ] .= _greek_mix(b_per_scatterer, coefs, :γ, b_total, FT, N)
            δ[:, iλ] .= _greek_mix(b_per_scatterer, coefs, :δ, b_total, FT, N)
            ϵ[:, iλ] .= _greek_mix(b_per_scatterer, coefs, :ϵ, b_total, FT, N)
            ζ[:, iλ] .= _greek_mix(b_per_scatterer, coefs, :ζ, b_total, FT, N)
        end
    end

    return polarized ?
        OceanLayerOptics(; λ, Δz, τ, ϖ, B, β, α, γ, δ, ϵ = ϵ, ζ) :
        OceanLayerOptics(; λ, Δz, τ, ϖ, B, β)
end

# Scattering-weighted mixing of one Greek coefficient across constituents.
# Returns a zero vector when total scattering is zero.
function _greek_mix(b_per_scatterer::AbstractVector{FT},
                    phase_matrices, field::Symbol, b_total::FT,
                    ::Type{FT}, N::Integer) where {FT}
    out = zeros(FT, N)
    b_total == zero(FT) && return out
    @inbounds for (b, mats) in zip(b_per_scatterer, phase_matrices)
        coef = getfield(mats, field)
        out .+= (b / b_total) .* coef
    end
    return out
end
