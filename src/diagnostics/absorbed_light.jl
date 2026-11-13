# =============================================================================
# Light absorbed per component, and the emission it drives
# =============================================================================
#
# A radiative-transfer solver supplies, for each layer and wavelength, the
# layer-mean spectral scalar irradiance
#
#     Ē₀(λ) = (1/Δz) ∫_layer E₀(z, λ) dz        [W m⁻² nm⁻¹]
#
# (E₀: total intensity integrated over all directions, direct and diffuse).
# In a homogeneous layer every component absorbs
#
#     A_c(λ) = a_c(λ) ∫ E₀ dz = τ_c(λ) Ē₀(λ)     [W m⁻² nm⁻¹],
#
# so the absorption components of `OceanLayerOptics` give the light absorbed
# by each constituent or pigment group. These helpers turn that into photon
# units, absorbed PAR, and the spectral power emitted by an inelastic
# process. Spectral integrals over a window [lo, hi] integrate the linear
# interpolant of the integrand on the layer's wavelength grid exactly (the
# trapezoidal rule, with the end intervals clipped at lo and hi); a solver
# with its own bin weights may integrate itself.
#
# A solver's irradiance can be negative at round-off level where the field
# has decayed. Entries no more negative than `roundoff_rtol` (default √eps of
# the irradiance's float type) times the largest |Ē₀| of the layer are set to
# zero; anything more negative is an error.
# =============================================================================

# h c N_A in J m mol⁻¹ (SI 2019 exact constants)
const _PLANCK_SI   = 6.62607015e-34    # J s
const _LIGHT_SI    = 299_792_458.0     # m s⁻¹
const _AVOGADRO_SI = 6.02214076e23     # mol⁻¹

"Photon flux per unit power at wavelength `λ` `[nm]`: µmol photons J⁻¹."
_µmol_per_joule(λ) =
    convert(typeof(float(λ)), _wide(λ) * 1e-9 / (_PLANCK_SI * _LIGHT_SI * _AVOGADRO_SI) * 1e6)

# Default relative size of negative irradiance entries treated as round-off: √eps of
# the irradiance's float type.
_roundoff_rtol(::AbstractVector{T}) where {T<:Real} = sqrt(eps(float(T)))

# The irradiance with round-off negatives set to zero; throws on anything else
# that is not a finite, non-negative spectrum on the layer grid.
function _checked_irradiance(o::OceanLayerOptics, E₀::AbstractVector{T}, who,
                             roundoff_rtol) where {T<:Real}
    roundoff_rtol ≥ 0 ||
        throw(ArgumentError("$who: roundoff_rtol must be non-negative, got $roundoff_rtol"))
    length(E₀) == n_wavelengths(o) ||
        throw(DimensionMismatch("$who: $(length(E₀)) irradiance values for " *
                                "$(n_wavelengths(o)) wavelengths"))
    all(isfinite, E₀) || throw(ArgumentError("$who: scalar irradiance must be finite"))
    allowance = roundoff_rtol * maximum(abs, E₀)
    lowest    = minimum(E₀)
    lowest ≥ -allowance ||
        throw(ArgumentError("$who: scalar irradiance must be non-negative; its lowest entry, " *
                            "$lowest, is below the round-off allowance −$allowance " *
                            "(roundoff_rtol = $roundoff_rtol times the largest entry). " *
                            "Check the sign and units of the field"))
    return [max(x, zero(x)) for x in E₀]
end

"""
    absorbed_power(o::OceanLayerOptics, E₀; roundoff_rtol = √eps) -> Matrix

Return the spectral power absorbed by each absorption component of layer `o`,
`A[c, k] = τ_c(λ_k) · E₀[k]`.

`E₀[k]` is the layer-mean spectral scalar irradiance at `o.λ[k]` (`(1/Δz) ∫ E₀ dz` over
the layer, total intensity over all directions, direct and diffuse). Rows follow
`o.absorbers`; units are those of `E₀` (`W m⁻² nm⁻¹` ⇒ absorbed `W m⁻² nm⁻¹` per layer).

`E₀` must be finite and non-negative. Negative entries at round-off level,
`E₀[k] ≥ −roundoff_rtol · maximum(abs, E₀)`, are treated as zero; more negative
entries throw an `ArgumentError`. The default `roundoff_rtol` is `√eps(T)` with `T` the
float type of `E₀` (`1.5e-8` in `Float64`); `roundoff_rtol = 0` rejects every negative
entry. The same rule and keyword apply in
[`absorbed_photons`](@ref), [`absorbed_par`](@ref) and
[`inelastic_emission`](@ref).
"""
function absorbed_power(o::OceanLayerOptics, E₀::AbstractVector{<:Real};
                        roundoff_rtol::Real = _roundoff_rtol(E₀))
    E = _checked_irradiance(o, E₀, "absorbed_power", roundoff_rtol)
    return [c.τ[k] * E[k] for c in o.absorbers, k in eachindex(o.λ)]
end

"""
    absorbed_photons(o::OceanLayerOptics, E₀) -> Matrix

Return [`absorbed_power`](@ref) in photon units: µmol photons m⁻² s⁻¹ nm⁻¹ per
absorption component, for `E₀` in `W m⁻² nm⁻¹` (photon energy `hc/λ`).
"""
function absorbed_photons(o::OceanLayerOptics, E₀::AbstractVector{<:Real};
                          roundoff_rtol::Real = _roundoff_rtol(E₀))
    A = absorbed_power(o, E₀; roundoff_rtol)
    return A .* _µmol_per_joule.(permutedims(o.λ))
end

"""
    absorbed_par(o::OceanLayerOptics, E₀; range = (400, 700), roundoff_rtol = √eps) -> Vector

Return the photosynthetically available radiation absorbed by each absorption component
of layer `o`: [`absorbed_photons`](@ref) integrated over `range` `[nm]`, in
µmol photons m⁻² s⁻¹.

The integral is that of the linear interpolant on the layer grid (trapezoidal rule, end
intervals clipped at the range limits), so the grid must span `range` but need not
contain its limits. Entries follow `o.absorbers`; see [`absorber_index`](@ref).

For phytoplankton with a pigment partition this gives the PAR absorbed by
each pigment group, e.g. by chlorophyll.
"""
function absorbed_par(o::OceanLayerOptics, E₀::AbstractVector{<:Real};
                      range::Tuple{<:Real, <:Real} = _PAR_RANGE_NM,
                      roundoff_rtol::Real = _roundoff_rtol(E₀))
    w = _window_weights(o.λ, range, "absorbed_par")
    N = absorbed_photons(o, E₀; roundoff_rtol)
    return N * w
end

# Weights w (on the points of λ, in either order) such that Σ w_k y_k is the
# integral over [lo, hi] of the linear interpolant of y: the trapezoidal
# rule, with the intervals that straddle lo or hi clipped there.
function _window_weights(λ::AbstractVector{FT}, range, who) where {FT}
    lo, hi = range
    lo < hi || throw(ArgumentError("$who: range must be (lo, hi) with lo < hi, got $range"))
    λmin, λmax = extrema(λ)
    λmin ≤ lo && λmax ≥ hi ||
        throw(ArgumentError("$who: the wavelength grid ($(λmin)–$(λmax) nm) does not span " *
                            "the integration range $(lo)–$(hi) nm"))
    order = sortperm(λ)
    w = zeros(FT, length(λ))
    for (a, b) in zip(@view(order[1:end-1]), @view(order[2:end]))
        xa, xb = λ[a], λ[b]
        u, v   = max(xa, lo), min(xb, hi)
        v > u || continue
        h      = xb - xa
        tu, tv = (u - xa) / h, (v - xa) / h
        wb     = h * (tv^2 - tu^2) / 2           # ∫ t dx over the clipped interval
        w[a]  += h * (tv - tu) - wb              # ∫ (1 − t) dx
        w[b]  += wb
    end
    return w
end

"""
    excitation_optical_depth(o::OceanLayerOptics, j) -> Vector

Return the excitation optical depth `τᴵ(λ) = aᴵ(λ) Δz` of the layer's `j`-th inelastic
component on the layer grid.

For a process with excitation pathways (fluorescence) it is formed from the stored
absorption components, `quantum_yield(p) · Σ_k efficiency[k] · absorbers[excitation[k]].τ`,
which equals `Δz · excitation_absorption(p, λ)`; otherwise (water Raman) it is
`Δz · excitation_absorption(p, λ)`.
"""
excitation_optical_depth(o::OceanLayerOptics{FT}, j::Integer) where {FT} =
    _excitation_optical_depth(o, o.inelastic[j])::Vector{FT}

# `o.inelastic` holds components of different process types, so the callers
# dispatch here at run time and assert the result type.
function _excitation_optical_depth(o::OceanLayerOptics{FT}, c::InelasticComponent) where {FT}
    isempty(c.excitation) &&
        return [FT(thickness(o) * excitation_absorption(c.process, λ)) for λ in o.λ]
    τ = zeros(FT, n_wavelengths(o))
    for (i, η) in zip(c.excitation, c.efficiency)
        τ .+= η .* o.absorbers[i].τ
    end
    return FT(quantum_yield(c.process)) .* τ
end

"""
    inelastic_emission(o::OceanLayerOptics, j, E₀, λ_em;
                       range = excitation_range(p), roundoff_rtol = √eps) -> Vector

Return the spectral power `[W m⁻² nm⁻¹]` emitted in all directions by layer `o` through
its `j`-th inelastic component (process `p`), at emission wavelengths `λ_em` `[nm]`:

```
P(λ) = Δz ∫ radiance_source_kernel(p, λ′, λ) Ē₀(λ′) dλ′
     = ∫ τᴵ(λ′) Ē₀(λ′) f(λ′, λ) λ′/λ dλ′      (photon conversion: fluorescence)
     = ∫ τᴵ(λ′) Ē₀(λ′) f(λ′, λ) dλ′           (energy conversion: default water Raman)
```

according to [`conversion_units`](@ref)`(p)`,
with `τᴵ` from [`excitation_optical_depth`](@ref) (the layer's stored
absorption components) and `E₀` the layer-mean scalar irradiance on the
layer grid (as in [`absorbed_power`](@ref)), integrated over `range` as
in [`absorbed_par`](@ref). For an isotropic process the volume source of
radiance is `P / (4π Δz)` per unit emission wavelength.

For fluorescence, the photons emitted equal the quantum yield times the
photons the excitation components deliver,
`∫ P λ/(hc) dλ = φ Σ_k η_k ∫ A_k λ′/(hc) dλ′` (to the extent `λ_em`
covers the emission band); reabsorption of the emission is left to the
radiative transfer.
"""
function inelastic_emission(o::OceanLayerOptics{FT}, j::Integer, E₀::AbstractVector{<:Real},
                            λ_em::AbstractVector{<:Real};
                            range::Tuple{<:Real, <:Real} = excitation_range(o.inelastic[j].process),
                            roundoff_rtol::Real = _roundoff_rtol(E₀)) where {FT}
    E = _checked_irradiance(o, E₀, "inelastic_emission", roundoff_rtol)
    T = promote_type(FT, eltype(E), float(eltype(λ_em)))
    return _inelastic_emission(T, o, o.inelastic[j], E, λ_em, range)::Vector{T}
end

function _inelastic_emission(::Type{T}, o::OceanLayerOptics, c::InelasticComponent, E, λ_em,
                             range) where {T}
    p  = c.process
    w  = _window_weights(o.λ, range, "inelastic_emission")
    x  = w .* _excitation_optical_depth(o, c) .* E        # excitation per grid point
    ks = [k for k in eachindex(o.λ) if x[k] != 0]
    u  = conversion_units(p)
    return T[sum(_energy_kernel(u, x[k] * emission(p, o.λ[k], λ), o.λ[k], λ) for k in ks;
                 init = zero(T))
             for λ in λ_em]
end
