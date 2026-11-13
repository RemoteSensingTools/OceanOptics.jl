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
# =============================================================================

# h c N_A in J m mol⁻¹ (SI 2019 exact constants)
const _PLANCK_SI   = 6.62607015e-34    # J s
const _LIGHT_SI    = 299_792_458.0     # m s⁻¹
const _AVOGADRO_SI = 6.02214076e23     # mol⁻¹

"Photon flux per unit power at wavelength `λ` `[nm]`: µmol photons J⁻¹."
_µmol_per_joule(λ) = λ * 1e-9 / (_PLANCK_SI * _LIGHT_SI * _AVOGADRO_SI) * 1e6

function _check_irradiance(o::OceanLayerOptics, E₀::AbstractVector, who)
    length(E₀) == n_wavelengths(o) ||
        throw(DimensionMismatch("$who: $(length(E₀)) irradiance values for " *
                                "$(n_wavelengths(o)) wavelengths"))
    all(x -> isfinite(x) && x ≥ 0, E₀) ||
        throw(ArgumentError("$who: scalar irradiance must be finite and non-negative"))
    return nothing
end

"""
    absorbed_power(o::OceanLayerOptics, E₀) -> Matrix

Spectral power absorbed by each absorption component of layer `o`,
`A[c, k] = τ_c(λ_k) · E₀[k]`, given the layer-mean spectral scalar
irradiance `E₀[k]` at `o.λ[k]` (`(1/Δz) ∫ E₀ dz` over the layer, total
intensity over all directions, direct and diffuse). Rows follow
`o.absorbers`; units are those of `E₀` (`W m⁻² nm⁻¹` ⇒ absorbed
`W m⁻² nm⁻¹` per layer).
"""
function absorbed_power(o::OceanLayerOptics, E₀::AbstractVector{<:Real})
    _check_irradiance(o, E₀, "absorbed_power")
    return [c.τ[k] * E₀[k] for c in o.absorbers, k in eachindex(o.λ)]
end

"""
    absorbed_photons(o::OceanLayerOptics, E₀) -> Matrix

[`absorbed_power`](@ref) in photon units: µmol photons m⁻² s⁻¹ nm⁻¹ per
absorption component, for `E₀` in `W m⁻² nm⁻¹` (photon energy `hc/λ`).
"""
function absorbed_photons(o::OceanLayerOptics, E₀::AbstractVector{<:Real})
    A = absorbed_power(o, E₀)
    return A .* _µmol_per_joule.(permutedims(o.λ))
end

"""
    absorbed_par(o::OceanLayerOptics, E₀; range = (400, 700)) -> Vector

Photosynthetically available radiation absorbed by each absorption
component of layer `o`: [`absorbed_photons`](@ref) integrated over
`range` `[nm]`, in µmol photons m⁻² s⁻¹. The integral is that of the
linear interpolant on the layer grid (trapezoidal rule, end intervals
clipped at the range limits), so the grid must span `range` but need not
contain its limits. Entries follow `o.absorbers`; see
[`absorber_index`](@ref).

For phytoplankton with a pigment partition this gives the PAR absorbed by
each pigment group, e.g. by chlorophyll.
"""
function absorbed_par(o::OceanLayerOptics, E₀::AbstractVector{<:Real};
                      range::Tuple{<:Real, <:Real} = (400, 700))
    w = _window_weights(o.λ, range, "absorbed_par")
    N = absorbed_photons(o, E₀)
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

Excitation optical depth `τᴵ(λ) = aᴵ(λ) Δz` of the layer's `j`-th
inelastic component on the layer grid. For a process with excitation
pathways (fluorescence) it is formed from the stored absorption
components, `quantum_yield(p) · Σ_k efficiency[k] · absorbers[excitation[k]].τ`,
which equals `Δz · excitation_absorption(p, λ)`; otherwise (water Raman)
it is `Δz · excitation_absorption(p, λ)`.
"""
function excitation_optical_depth(o::OceanLayerOptics{FT}, j::Integer) where {FT}
    c = o.inelastic[j]
    isempty(c.excitation) &&
        return [FT(thickness(o) * excitation_absorption(c.process, λ)) for λ in o.λ]
    τ = zeros(FT, n_wavelengths(o))
    for (i, η) in zip(c.excitation, c.efficiency)
        τ .+= η .* o.absorbers[i].τ
    end
    return FT(quantum_yield(c.process)) .* τ
end

"""
    inelastic_emission(o::OceanLayerOptics, j, E₀, λ_em; range = excitation_range(p)) -> Vector

Spectral power `[W m⁻² nm⁻¹]` emitted in all directions by layer `o`
through its `j`-th inelastic component (process `p`), at emission
wavelengths `λ_em` `[nm]`:

`P(λ) = ∫ τᴵ(λ′) Ē₀(λ′) f(λ′, λ) λ′/λ dλ′ = Δz ∫ radiance_source_kernel(p, λ′, λ) Ē₀(λ′) dλ′`,

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
                            range::Tuple{<:Real, <:Real} = excitation_range(o.inelastic[j].process)) where {FT}
    _check_irradiance(o, E₀, "inelastic_emission")
    p  = o.inelastic[j].process
    w  = _window_weights(o.λ, range, "inelastic_emission")
    x  = w .* excitation_optical_depth(o, j) .* E₀        # excitation per grid point
    ks = [k for k in eachindex(o.λ) if x[k] != 0]
    return [sum(x[k] * emission(p, o.λ[k], λ) * o.λ[k] / λ for k in ks; init = zero(FT) * zero(float(λ)))
            for λ in λ_em]
end
