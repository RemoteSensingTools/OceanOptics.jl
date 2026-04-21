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
    layer_optics(layer, λ_grid; ℓ_max=64, FT=Float64) -> OceanLayerOptics

Compute the complete optical bundle for one ocean layer at each wavelength
in `λ_grid`.

# Arguments
- `layer::OceanLayer` — the layer (constituents + thickness)
- `λ_grid` — vector of wavelengths [nm]; plain `AbstractVector{<:Real}` or
  `Unitful` with spectral units
- `ℓ_max` — maximum Legendre order for the phase-function expansion
  (default 64, which is adequate for Gauss-Lobatto RT with N ≤ 32)
- `FT` — floating-point precision

# Returns
An [`OceanLayerOptics`](@ref) containing τ, ϖ, β_ℓ(λ) and the diagnostic
backscatter fraction B, all on the `λ_grid`.

# Notes
This function handles only the *elastic* IOPs. Inelastic source terms
(fluorescence, Raman) are computed separately by
[`inelastic_sources`](@ref) from a scalar irradiance field `E°(τ, λ)`
obtained by a first-pass elastic solve (see Fell 1997 §2.2.9).

# Example

```julia
water = PureWater{Float64}()
chl   = Phytoplankton{Float64, Bricaud1998}(Chl = 0.5)
layer = OceanLayer(0.0, 5.0, [water, chl])
λ     = 400.0:5.0:700.0

opt = layer_optics(layer, collect(λ); ℓ_max = 64)
# opt.τ, opt.ϖ, opt.β all have one column/entry per λ.
```
"""
function layer_optics(layer::OceanLayer, λ_grid::AbstractVector{<:Real};
                      ℓ_max::Integer = 64,
                      FT::Type = Float64)
    @assert ℓ_max ≥ 0 "ℓ_max must be non-negative"
    nλ = length(λ_grid)
    Δz = FT(thickness(layer))
    λ  = FT.(λ_grid)

    # Allocate outputs
    τ = zeros(FT, nλ)
    ϖ = zeros(FT, nλ)
    B = zeros(FT, nλ)
    β = zeros(FT, ℓ_max + 1, nλ)

    # Split constituents by role — avoids re-dispatching per λ inside the loop
    # and lets us extract phase functions once.
    scatterers    = filter(c -> c isa Union{AbstractScatterer, AbstractAbsorbingScatterer},
                           layer.constituents)
    phase_fns     = [phase_function(c) for c in scatterers]
    # Pre-compute Legendre moments once per scattering constituent — they
    # are wavelength-independent in most practical ocean cases (phase
    # functions depend on λ only through the backscatter fraction, which
    # itself varies weakly; future wavelength-dependent phase functions
    # would shift this computation inside the λ loop).
    phase_moments = [phase_function_moments(pf, ℓ_max) for pf in phase_fns]

    @inbounds for (iλ, λ_i) in enumerate(λ)
        # Absorption + scattering coefficients per constituent
        a_total  = zero(FT)
        b_total  = zero(FT)
        bb_total = zero(FT)
        # Scattering coefficient per scattering constituent — needed for
        # phase-function mixing below
        b_per_scatterer = zeros(FT, length(scatterers))

        for c in layer.constituents
            a_total += FT(absorption(c, λ_i))
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

        # Scattering-weighted phase-moment mixing. Pass the per-layer
        # moment count explicitly so empty-scatterer layers return a
        # correctly-sized (β₀ = 1, β_{ℓ≥1} = 0) vector rather than the
        # length-1 fallback, which would broadcast-fill all rows with 1.
        β_iλ = mix_phase_moments(b_per_scatterer, phase_moments;
                                 n_moments = ℓ_max + 1)
        β[:, iλ] .= β_iλ
    end

    return OceanLayerOptics(; λ, Δz, τ, ϖ, B, β)
end
