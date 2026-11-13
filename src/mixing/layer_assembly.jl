# =============================================================================
# Layer and column assembly
# =============================================================================
#
# For one layer and wavelength grid:
#   1. one AbsorptionComponent per absorber and absorption part:
#      τ = a_part Δz (with the total-absorption exclusivity rule);
#   2. one ScatteringComponent per scattering constituent: τ_c = b_c Δz, the
#      phase function's expansion (once, or per wavelength if it depends on
#      wavelength), its backscatter fraction and a polarization label;
#   3. one InelasticComponent per inelastic process, linked to the
#      absorption components of its excitation pathways.
# =============================================================================

"""
    layer_optics(layer, λ; ℓ_max = 64, polarized = false,
                 scalar_phase_approximation = nothing, FT = Float64) -> OceanLayerOptics

Optical properties of one [`OceanLayer`](@ref) at wavelengths `λ` `[nm]`
(strictly monotone, increasing or decreasing).

# Keywords
- `ℓ_max` — highest expansion order. The expansions are exact through
  `ℓ_max` (to the quadrature tolerance) and not truncated further: δ-M or
  similar truncation is the solver's choice and needs `ℓ_max` at least as
  large as its fitting order.
- `polarized` — also compute phase-matrix coefficients `α, γ, δ, ϵ, ζ`.
- `scalar_phase_approximation` — with `polarized = true`, how to treat
  scatterers whose phase function is scalar-only (Fournier-Forand,
  Henyey-Greenstein, Petzold): `nothing` throws; `:depolarizing` uses the
  completely depolarizing approximation, labelled `:depolarizing` on the
  component. Applies to every such scatterer, even one with zero
  scattering. Ignored when `polarized = false`.
- `FT` — floating-point type of the result.

# Result
An [`OceanLayerOptics`](@ref) with, in constituent order,
- one [`AbsorptionComponent`](@ref) per absorbing constituent and
  absorption part (see [`absorption_parts`](@ref): one part, `:total`,
  unless e.g. phytoplankton carries a pigment partition);
- one [`ScatteringComponent`](@ref) per scattering constituent;
- one [`InelasticComponent`](@ref) per process in `layer.inelastic`, with
  the indices of the absorption components that excite it.
Mixed quantities follow from [`absorption_optical_depth`](@ref),
[`optical_depth`](@ref), [`single_scattering_albedo`](@ref) and
[`greek_expansion`](@ref)`(o, iλ)`.

If a constituent declares [`is_total_absorption`](@ref) (a
`TabulatedAbsorber` with `total = true`), it must be the only one, and
every other constituent must have zero absorption at every wavelength.

Every absorber that excites a fluorescence process must itself be a
constituent of the layer (the same object, or an equal one): otherwise its
absorption would drive emission without attenuating the light.

# Example
```julia
water = PureWater{Float64}()
phyto = Phytoplankton{Float64, Bricaud1998}(Chl = 0.5)
layer = OceanLayer(0.0, 5.0, AbstractOceanConstituent{Float64}[water, phyto])
opt   = layer_optics(layer, collect(400.0:5.0:700.0); ℓ_max = 64)
optical_depth(opt), single_scattering_albedo(opt)
opt.scatterers[2].label               # "Phytoplankton"
```
"""
function layer_optics(layer::OceanLayer, λ_grid::AbstractVector{<:Real};
                      ℓ_max::Integer = 64,
                      polarized::Bool = false,
                      scalar_phase_approximation::Union{Nothing, Symbol} = nothing,
                      FT::Type{<:AbstractFloat} = Float64)
    return _layer_optics(layer, λ_grid, ℓ_max, polarized, scalar_phase_approximation, FT,
                         Dict{Any, Any}())
end

# `cache` maps (phase_id, ℓ_max, FT, phase key) to (expansion, backscatter
# fraction), so that each distinct expansion is computed once — per layer,
# or per column when `column_optics` shares it between layers (which then
# share the same GreekCoefficients objects).
function _layer_optics(layer::OceanLayer, λ_grid, ℓ_max, polarized,
                       scalar_phase_approximation, FT, cache)
    ℓ_max ≥ 0 || throw(ArgumentError("ℓ_max must be non-negative, got $ℓ_max"))
    scalar_phase_approximation in (nothing, :depolarizing) ||
        throw(ArgumentError("scalar_phase_approximation must be `nothing` or " *
                            ":depolarizing, got :$scalar_phase_approximation"))
    λ = collect(FT, λ_grid)
    !isempty(λ) && all(x -> isfinite(x) && x > 0, λ) && _is_strictly_monotone(λ) ||
        throw(ArgumentError("layer_optics: wavelengths must be finite, positive and " *
                            "strictly monotone"))
    Δz = FT(thickness(layer))
    constituents = layer.constituents

    # --- absorption components -----------------------------------------------
    absorbing = [c for c in constituents if c isa AbstractAbsorber]
    i_total   = findall(is_total_absorption, absorbing)
    length(i_total) ≤ 1 ||
        throw(ArgumentError("layer_optics: $(length(i_total)) constituents declare total " *
                            "(bulk) absorption; at most one can"))
    i_bulk    = isempty(i_total) ? 0 : only(i_total)
    absorbers = AbsorptionComponent{FT}[]
    first_part = Int[]            # index of each absorber's first component
    for (i, c) in enumerate(absorbing)
        parts = absorption_parts(c)
        τ     = zeros(FT, length(parts), length(λ))
        for (k, λk) in enumerate(λ)
            a = absorption_by_part(c, λk)
            if i_bulk != 0 && i != i_bulk && any(!iszero, a)
                throw(ArgumentError("layer_optics: $(absorbing[i_bulk]) is the total " *
                                    "absorption of the layer, but $(nameof(typeof(c))) also " *
                                    "absorbs (a = $(sum(a)) 1/m at λ = $λk nm); this would " *
                                    "count its absorbers twice"))
            end
            for j in eachindex(parts)
                τ[j, k] = FT(a[j]) * Δz
            end
        end
        push!(first_part, length(absorbers) + 1)
        for (j, part) in enumerate(parts)
            push!(absorbers, AbsorptionComponent{FT}(_component_label(c), part, τ[j, :]))
        end
    end

    # --- inelastic components ---------------------------------------------
    inelastic = map(layer.inelastic) do p
        pathways = excitation_pathways(p)
        index    = map(pathways) do x
            i = _matching_absorber(absorbing, x.absorber, p, layer)
            first_part[i] + findfirst(==(x.part), absorption_parts(absorbing[i])) - 1
        end
        InelasticComponent{FT}(string(nameof(typeof(p))), p, collect(Int, index),
                               FT[x.efficiency for x in pathways])
    end

    # --- scattering components -------------------------------------------
    scatterers = filter(c -> c isa Union{AbstractScatterer, AbstractAbsorbingScatterer},
                        constituents)
    if polarized && scalar_phase_approximation === nothing
        scalar_only = [string(nameof(typeof(c)), " (", nameof(typeof(scattering_phase(c))), ")")
                       for c in scatterers if !is_polarizable(scattering_phase(c))]
        isempty(scalar_only) ||
            throw(ArgumentError("layer_optics(polarized = true): scalar-only phase functions " *
                                "for " * join(scalar_only, ", ") * ". Their phase matrices " *
                                "are unknown; supply a phase matrix (TabulatedScatterer with " *
                                "a TabulatedMuellerPhase or polarized TabulatedGreekPhase) or " *
                                "pass `scalar_phase_approximation = :depolarizing`."))
    end
    components = [_scattering_component(c, λ, Δz, ℓ_max, polarized, FT, cache)
                  for c in scatterers]
    return OceanLayerOptics{FT}(λ, FT(layer.depth_top), FT(layer.depth_bottom), absorbers,
                                components, collect(InelasticComponent{FT}, inelastic),
                                ℓ_max + 1, polarized)
end

function _scattering_component(c, λ::Vector{FT}, Δz::FT, ℓ_max, polarized, ::Type{FT},
                               cache) where {FT}
    pf    = scattering_phase(c)
    label = !polarized ? :scalar : is_polarizable(pf) ? :full : :depolarizing
    pid   = hash((objectid(pf), label))
    expand(args...) = label === :scalar ? greek_expansion(pf, ℓ_max, args...) :
                      label === :full   ? greek_expansion(pf, ℓ_max, args...; polarized = true) :
                                          depolarizing_greek_expansion(pf, ℓ_max, args...)
    entry(key, args...) = get!(cache, (pid, ℓ_max, FT, key)) do
        (_convert_greek(FT, expand(args...)), FT(backscatter_fraction(pf, args...)))
    end
    if is_wavelength_dependent(pf)
        keys  = [_phase_key(pf, λk) for λk in λ]
        slot  = Dict{Any, Int}()
        first_λ = Int[]
        for (k, key) in enumerate(keys)
            haskey(slot, key) || (slot[key] = length(slot) + 1; push!(first_λ, k))
        end
        entries = [entry(keys[k], λ[k]) for k in first_λ]
        index   = length(entries) == 1 ? Int[] : [slot[key] for key in keys]
    else
        entries = [entry(nothing)]
        index   = Int[]
    end
    greek = [e[1] for e in entries]
    B     = [e[2] for e in entries]
    τ     = [FT(scattering(c, λk)) * Δz for λk in λ]
    return ScatteringComponent{FT}(_component_label(c), pid, τ, greek, B, index, label)
end

_component_label(c::AbstractOceanConstituent) = string(nameof(typeof(c)))

# The constituent of the layer that is the exciting absorber: the same object,
# or else an equal one; it must be unique, so that excitation is attributed
# to one absorption component.
function _matching_absorber(absorbing, absorber, p, layer)
    where_ = "$(nameof(typeof(p))) in the layer at $(layer.depth_top)–$(layer.depth_bottom) m"
    i = findall(c -> c === absorber, absorbing)
    isempty(i) && (i = findall(c -> isequal(c, absorber), absorbing))
    isempty(i) &&
        throw(ArgumentError("layer_optics: $where_ is excited by $absorber, which is not a " *
                            "constituent of that layer. Add it to the constituents, so that " *
                            "the light it absorbs is also removed from the radiation field."))
    length(i) == 1 ||
        throw(ArgumentError("layer_optics: $where_ is excited by $absorber, which appears " *
                            "$(length(i)) times among the constituents; merge the equal " *
                            "constituents into one"))
    return only(i)
end

_convert_greek(::Type{FT}, g::GreekCoefficients{FT}) where {FT} = g
_convert_greek(::Type{FT}, g::GreekCoefficients) where {FT} =
    GreekCoefficients{FT}((FT.(getfield(g, f)) for f in fieldnames(GreekCoefficients))...)

"""
    column_optics(column, λ; refractive_index, bottom, ℓ_max = 64, polarized = false,
                  scalar_phase_approximation = nothing, FT = Float64) -> OceanColumnOptics

Optical description of an [`OceanColumn`](@ref) for a coupled
atmosphere–ocean solver: [`layer_optics`](@ref) for every layer on the
same grid and settings, the interface refractive index, and the bottom.

Components with the same phase function and polarization treatment share
one `phase_id` and the same expansion objects across layers.

`refractive_index` (an [`AbstractRefractiveIndexModel`](@ref), e.g.
`FixedRefractiveIndex(1.34)` or `QuanFry1995(; temperature, salinity)`)
and `bottom` are required: neither has a defensible default. The bottom
is `BlackBottom()` (deep water), `LambertianBottom(albedo)` (one albedo,
or one per wavelength of `λ`) or a [`TabulatedLambertianBottom`](@ref),
which is interpolated onto `λ`.

# Example
```julia
water  = PureWater{Float64}()
column = uniform_column(AbstractOceanConstituent{Float64}[water]; n_layers = 4, depth = 20.0)
col    = column_optics(column, [450.0, 550.0]; refractive_index = FixedRefractiveIndex(1.34),
                       bottom = BlackBottom(), polarized = true)
```
"""
function column_optics(column::OceanColumn, λ_grid::AbstractVector{<:Real};
                       refractive_index::AbstractRefractiveIndexModel,
                       bottom::AbstractOceanBottom,
                       ℓ_max::Integer = 64,
                       polarized::Bool = false,
                       scalar_phase_approximation::Union{Nothing, Symbol} = nothing,
                       FT::Type{<:AbstractFloat} = Float64)
    cache  = Dict{Any, Any}()
    layers = [_layer_optics(l, λ_grid, ℓ_max, polarized, scalar_phase_approximation, FT, cache)
              for l in column]
    λ = first(layers).λ
    n = FT.(OceanOptics.refractive_index(refractive_index, λ))
    return OceanColumnOptics(layers, n, _resolve_bottom(bottom, λ, FT))
end

_resolve_bottom(b::AbstractOceanBottom, λ, FT) = b
_resolve_bottom(b::LambertianBottom, λ, FT) = LambertianBottom(FT.(b.albedo))
function _resolve_bottom(b::TabulatedLambertianBottom, λ, FT)
    lo, hi = first(b.λ), last(b.λ)
    all(x -> lo ≤ x ≤ hi, λ) ||
        throw(DomainError(extrema(λ), "TabulatedLambertianBottom" *
                          (isempty(b.source) ? "" : " ($(b.source))") * " covers $(lo)–$(hi) nm, " *
                          "the column grid $(minimum(λ))–$(maximum(λ)) nm"))
    tbl = SpectralTable(FT.(b.λ), FT.(b.albedo))
    return LambertianBottom(FT[linterp(tbl, x) for x in λ])
end
