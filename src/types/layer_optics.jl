# =============================================================================
# Layer and column optics — the solver-facing output
# =============================================================================
#
# A layer is described the way plane-parallel solvers combine optical
# properties: one absorption optical depth plus one scattering component
# per scattering constituent, each with its own scattering optical depth
# and phase-matrix expansion. Mixed quantities (total τ, ϖ, the mixed
# expansion) are derived on demand. Keeping components separate lets a
# solver build each component's phase operators once per band and truncate
# them individually, as it does for aerosol types, and it keeps the
# provenance of every phase matrix visible.
# =============================================================================

_is_strictly_increasing(x::AbstractVector) = all(i -> x[i] < x[i+1], firstindex(x):lastindex(x)-1)
_is_strictly_monotone(x::AbstractVector) =
    _is_strictly_increasing(x) || all(i -> x[i] > x[i+1], firstindex(x):lastindex(x)-1)

const _POLARIZATION_LABELS = (:scalar, :full, :depolarizing)

"""
$(TYPEDEF)

Scattering by one constituent of a layer.

# Fields
$(TYPEDFIELDS)

`greek` and `backscatter_fraction` hold the component's *distinct*
expansions: one when the phase function does not depend on wavelength;
otherwise `phase_index[iλ]` names the entry used at each wavelength of the
layer. See [`greek_expansion`](@ref)`(c, iλ)` and
[`is_wavelength_dependent`](@ref)`(c)`.

`phase_id` identifies the phase function and polarization treatment:
components with equal `phase_id` (in any layer of a column) have identical
expansions — [`column_optics`](@ref) shares the same objects — so a solver
can build their phase operators once.

`polarization` records how the phase-matrix coefficients were obtained:
`:scalar` (phase function only, `β`), `:full` (the constituent's phase
matrix), or `:depolarizing` (its phase function with the completely
depolarizing approximation, selected explicitly).
"""
struct ScatteringComponent{FT}
    "Constituent the component came from (type name, or the source of tabulated data)"
    label::String
    "Identity of the phase function and polarization treatment (equal ⇒ identical expansions)"
    phase_id::UInt64
    "Scattering optical depth `b·Δz` per wavelength"
    τ::Vector{FT}
    "Distinct phase-matrix expansions"
    greek::Vector{GreekCoefficients{FT}}
    "Backscatter fraction of the phase function (not of its truncated expansion), per entry of `greek`"
    backscatter_fraction::Vector{FT}
    "Entry of `greek` used at each wavelength; empty when there is a single entry"
    phase_index::Vector{Int}
    "`:scalar`, `:full` or `:depolarizing`"
    polarization::Symbol

    function ScatteringComponent{FT}(label, phase_id, τ, greek, backscatter_fraction,
                                     phase_index, polarization) where {FT}
        c = new{FT}(label, phase_id, τ, greek, backscatter_fraction, phase_index, polarization)
        _validate_component(c)
        return c
    end
end

function ScatteringComponent(label::AbstractString, τ::AbstractVector{FT},
                             greek::AbstractVector{GreekCoefficients{FT}}, B::AbstractVector,
                             polarization::Symbol; phase_index::AbstractVector{<:Integer} = Int[],
                             phase_id::Integer = hash((label, polarization))) where {FT}
    return ScatteringComponent{FT}(String(label), UInt64(phase_id), collect(τ), collect(greek),
                                   collect(FT, B), collect(Int, phase_index), polarization)
end

function _validate_component(c::ScatteringComponent)
    who = "ScatteringComponent \"$(c.label)\""
    all(x -> isfinite(x) && x ≥ 0, c.τ) ||
        throw(ArgumentError("$who: scattering optical depths must be finite and non-negative"))
    isempty(c.greek) && throw(ArgumentError("$who: no phase-matrix expansion"))
    length(c.backscatter_fraction) == length(c.greek) ||
        throw(DimensionMismatch("$who: $(length(c.backscatter_fraction)) backscatter " *
                                "fractions for $(length(c.greek)) expansions"))
    all(x -> isfinite(x) && 0 ≤ x ≤ 1, c.backscatter_fraction) ||
        throw(ArgumentError("$who: backscatter fractions must lie in [0, 1]"))
    c.polarization in _POLARIZATION_LABELS ||
        throw(ArgumentError("$who: polarization must be one of $(_POLARIZATION_LABELS)"))
    if isempty(c.phase_index)
        length(c.greek) == 1 ||
            throw(ArgumentError("$who: $(length(c.greek)) expansions need a phase_index"))
    else
        length(c.phase_index) == length(c.τ) ||
            throw(DimensionMismatch("$who: phase_index has $(length(c.phase_index)) entries " *
                                    "for $(length(c.τ)) wavelengths"))
        all(i -> 1 ≤ i ≤ length(c.greek), c.phase_index) ||
            throw(ArgumentError("$who: phase_index refers to missing expansions"))
    end
    g1 = first(c.greek)
    all(g -> n_moments(g) == n_moments(g1), c.greek) ||
        throw(ArgumentError("$who: expansions differ in length"))
    pol = c.polarization != :scalar
    all(g -> is_polarized(g) == pol, c.greek) ||
        throw(ArgumentError("$who: polarization label :$(c.polarization) does not match " *
                            "the expansions"))
    return nothing
end

"Number of expansion terms (`ℓ_max + 1`)."
n_moments(c::ScatteringComponent) = n_moments(first(c.greek))
"Whether the component's phase matrix differs between wavelengths."
is_wavelength_dependent(c::ScatteringComponent) = !isempty(c.phase_index)
"Whether the component carries phase-matrix coefficients."
is_polarized(c::ScatteringComponent) = c.polarization != :scalar

"""
    greek_expansion(c::ScatteringComponent, iλ) -> GreekCoefficients

Expansion of component `c` at the `iλ`-th wavelength of its layer.
"""
greek_expansion(c::ScatteringComponent, iλ::Integer) =
    is_wavelength_dependent(c) ? c.greek[c.phase_index[iλ]] : first(c.greek)

backscatter_fraction(c::ScatteringComponent, iλ::Integer) =
    is_wavelength_dependent(c) ? c.backscatter_fraction[c.phase_index[iλ]] :
                                 first(c.backscatter_fraction)

# =============================================================================
# Layer
# =============================================================================

"""
$(TYPEDEF)

Optical properties of one homogeneous ocean layer on a wavelength grid,
produced by [`layer_optics`](@ref).

# Fields
$(TYPEDFIELDS)

# Derived quantities

[`optical_depth`](@ref), [`scattering_optical_depth`](@ref),
[`single_scattering_albedo`](@ref), [`backscatter_fraction`](@ref) and
the mixed expansion [`greek_expansion`](@ref)`(o, iλ)` — the
scattering-weighted mixture of the components.

# Validation

`λ` finite, positive and strictly monotone (increasing or decreasing);
`0 ≤ depth_top < depth_bottom`; optical depths finite and non-negative;
every component on the same wavelength grid, with `n_moments` terms and
the layer's polarization. Arrays are checked at construction; mutating
them afterwards is not re-checked.
"""
struct OceanLayerOptics{FT}
    "Wavelengths `[nm]`, strictly monotone"
    λ::Vector{FT}
    "Depth of the layer top below the surface `[m]`"
    depth_top::FT
    "Depth of the layer bottom below the surface `[m]`"
    depth_bottom::FT
    "Absorption optical depth `a·Δz` per wavelength (all absorbers)"
    τ_abs::Vector{FT}
    "Scattering components, one per scattering constituent"
    components::Vector{ScatteringComponent{FT}}
    "Number of expansion terms, `ℓ_max + 1`"
    n_moments::Int
    "Whether every component carries phase-matrix coefficients"
    polarized::Bool

    function OceanLayerOptics{FT}(λ, depth_top, depth_bottom, τ_abs, components,
                                  n_moments, polarized) where {FT}
        o = new{FT}(λ, depth_top, depth_bottom, τ_abs, components, n_moments, polarized)
        _validate_layer(o)
        return o
    end
end

function _validate_layer(o::OceanLayerOptics)
    who = "OceanLayerOptics"
    nλ  = length(o.λ)
    nλ ≥ 1 || throw(ArgumentError("$who: empty wavelength grid"))
    all(x -> isfinite(x) && x > 0, o.λ) ||
        throw(ArgumentError("$who: wavelengths must be finite and positive"))
    _is_strictly_monotone(o.λ) ||
        throw(ArgumentError("$who: wavelengths must be strictly increasing or strictly decreasing"))
    isfinite(o.depth_top) && isfinite(o.depth_bottom) && 0 ≤ o.depth_top < o.depth_bottom ||
        throw(ArgumentError("$who: need 0 ≤ depth_top < depth_bottom, got " *
                            "$(o.depth_top)–$(o.depth_bottom) m"))
    length(o.τ_abs) == nλ ||
        throw(DimensionMismatch("$who: $(length(o.τ_abs)) absorption optical depths for $nλ wavelengths"))
    all(x -> isfinite(x) && x ≥ 0, o.τ_abs) ||
        throw(ArgumentError("$who: absorption optical depths must be finite and non-negative"))
    o.n_moments ≥ 1 || throw(ArgumentError("$who: n_moments must be positive"))
    for c in o.components
        length(c.τ) == nλ ||
            throw(DimensionMismatch("$who: component \"$(c.label)\" has $(length(c.τ)) " *
                                    "optical depths for $nλ wavelengths"))
        n_moments(c) == o.n_moments ||
            throw(ArgumentError("$who: component \"$(c.label)\" has $(n_moments(c)) terms, " *
                                "layer has $(o.n_moments)"))
        is_polarized(c) == o.polarized ||
            throw(ArgumentError("$who: component \"$(c.label)\" is :$(c.polarization) in a " *
                                (o.polarized ? "polarized" : "scalar") * " layer"))
    end
    return nothing
end

"Geometric thickness `[m]`."
thickness(o::OceanLayerOptics) = o.depth_bottom - o.depth_top
"Number of wavelengths."
n_wavelengths(o::OceanLayerOptics) = length(o.λ)
"Number of expansion terms, `ℓ_max + 1`."
n_moments(o::OceanLayerOptics) = o.n_moments

"""
    is_polarized(x) -> Bool

Whether `x` (layer, column, component or [`GreekCoefficients`](@ref))
carries all six phase-matrix coefficients.
"""
is_polarized(o::OceanLayerOptics) = o.polarized

"Absorption optical depth per wavelength."
absorption_optical_depth(o::OceanLayerOptics) = o.τ_abs

"Scattering optical depth per wavelength, summed over components."
function scattering_optical_depth(o::OceanLayerOptics{FT}) where {FT}
    τ = zeros(FT, n_wavelengths(o))
    for c in o.components
        τ .+= c.τ
    end
    return τ
end

"Total (extinction) optical depth `τ_abs + τ_scat` per wavelength."
optical_depth(o::OceanLayerOptics) = o.τ_abs .+ scattering_optical_depth(o)

"""
    single_scattering_albedo(o::OceanLayerOptics) -> Vector

`ϖ = τ_scat / (τ_abs + τ_scat)` per wavelength; zero where `τ = 0`.
"""
function single_scattering_albedo(o::OceanLayerOptics{FT}) where {FT}
    τs = scattering_optical_depth(o)
    τ  = o.τ_abs .+ τs
    return [τ[i] > 0 ? τs[i] / τ[i] : zero(FT) for i in eachindex(τ)]
end

"""
    backscatter_fraction(o::OceanLayerOptics) -> Vector

Scattering-weighted backscatter fraction of the components' phase
functions per wavelength; zero where the layer does not scatter.
"""
function backscatter_fraction(o::OceanLayerOptics{FT}) where {FT}
    τs = scattering_optical_depth(o)
    B  = zeros(FT, n_wavelengths(o))
    for c in o.components, i in eachindex(B)
        B[i] += c.τ[i] * backscatter_fraction(c, i)
    end
    return [τs[i] > 0 ? B[i] / τs[i] : zero(FT) for i in eachindex(B)]
end

"""
    greek_expansion(o::OceanLayerOptics, iλ) -> GreekCoefficients

Expansion of the layer's mixed phase matrix at its `iλ`-th wavelength: the
components weighted by their scattering optical depths. Where the layer
does not scatter, the isotropic expansion `β = [1, 0, …]` (zero
polarization terms).
"""
function greek_expansion(o::OceanLayerOptics{FT}, iλ::Integer) where {FT}
    w = [c.τ[iλ] for c in o.components]
    if isempty(w) || sum(w) == 0
        β = zeros(FT, o.n_moments); β[1] = one(FT)
        g = GreekCoefficients(β)
        return o.polarized ? depolarized(g) : g
    end
    return mix_greek(w, [greek_expansion(c, iλ) for c in o.components])
end

function Base.show(io::IO, o::OceanLayerOptics{FT}) where {FT}
    print(io, "OceanLayerOptics{$FT}($(o.depth_top)–$(o.depth_bottom) m, ",
          n_wavelengths(o), " λ, ", length(o.components), " scattering components, ",
          "ℓ ≤ $(o.n_moments - 1), ", o.polarized ? "polarized" : "scalar", ")")
end

# =============================================================================
# Bottom boundary descriptors
# =============================================================================

"Lower boundary of an ocean column (for the solver to apply)."
abstract type AbstractOceanBottom end

"Perfectly absorbing bottom."
struct BlackBottom <: AbstractOceanBottom end

"""
$(TYPEDEF)

Lambertian (isotropically reflecting, depolarizing) bottom with albedo per
wavelength (or one value for all wavelengths), in `[0, 1]`.

# Fields
$(TYPEDFIELDS)
"""
struct LambertianBottom{FT} <: AbstractOceanBottom
    "Albedo: one value, or one per wavelength of the column"
    albedo::Vector{FT}

    function LambertianBottom{FT}(albedo::Vector{FT}) where {FT}
        isempty(albedo) && throw(ArgumentError("LambertianBottom: no albedo"))
        all(a -> isfinite(a) && 0 ≤ a ≤ 1, albedo) ||
            throw(ArgumentError("LambertianBottom: albedo must lie in [0, 1]"))
        new{FT}(albedo)
    end
end

LambertianBottom(albedo::Real) = LambertianBottom{float(typeof(albedo))}([float(albedo)])
LambertianBottom(albedo::AbstractVector{<:Real}) =
    LambertianBottom{float(eltype(albedo))}(collect(float(eltype(albedo)), albedo))

# =============================================================================
# Column
# =============================================================================

"""
$(TYPEDEF)

Optical description of a whole water column, produced by
[`column_optics`](@ref): contiguous layers from the surface down on one
wavelength grid, the refractive index of the water at the surface relative
to air, and the bottom.

# Fields
$(TYPEDFIELDS)

Layers share `λ`, `n_moments` and polarization; the first starts at the
surface and each starts where the previous one ends. `refractive_index`
is real and given per wavelength.
"""
struct OceanColumnOptics{FT, B<:AbstractOceanBottom}
    "Layers, surface first"
    layers::Vector{OceanLayerOptics{FT}}
    "Refractive index of the surface water relative to air, per wavelength"
    refractive_index::Vector{FT}
    "Lower boundary"
    bottom::B

    function OceanColumnOptics{FT, B}(layers, refractive_index, bottom) where {FT, B}
        col = new{FT, B}(layers, refractive_index, bottom)
        _validate_column(col)
        return col
    end
end

OceanColumnOptics(layers::Vector{OceanLayerOptics{FT}}, n::AbstractVector,
                  bottom::B) where {FT, B<:AbstractOceanBottom} =
    OceanColumnOptics{FT, B}(layers, collect(FT, n), bottom)

function _validate_column(col::OceanColumnOptics)
    who = "OceanColumnOptics"
    isempty(col.layers) && throw(ArgumentError("$who: no layers"))
    l1 = first(col.layers)
    l1.depth_top == 0 || throw(ArgumentError("$who: the first layer must start at the surface"))
    for (i, l) in enumerate(col.layers)
        l.λ == l1.λ || throw(ArgumentError("$who: layer $i has a different wavelength grid"))
        l.n_moments == l1.n_moments ||
            throw(ArgumentError("$who: layer $i has $(l.n_moments) expansion terms, layer 1 $(l1.n_moments)"))
        l.polarized == l1.polarized ||
            throw(ArgumentError("$who: layers differ in polarization"))
        i == 1 || l.depth_top == col.layers[i-1].depth_bottom ||
            throw(ArgumentError("$who: layers $(i-1) and $i are not contiguous"))
    end
    nλ = length(l1.λ)
    length(col.refractive_index) == nλ ||
        throw(DimensionMismatch("$who: $(length(col.refractive_index)) refractive indices for $nλ wavelengths"))
    all(n -> isfinite(n) && n > 0, col.refractive_index) ||
        throw(ArgumentError("$who: refractive indices must be finite and positive"))
    if col.bottom isa LambertianBottom
        length(col.bottom.albedo) in (1, nλ) ||
            throw(DimensionMismatch("$who: bottom albedo needs 1 or $nλ values"))
    end
    return nothing
end

"Number of layers."
n_layers(col::OceanColumnOptics) = length(col.layers)
"Depths of the layer interfaces `[m]`, from the surface (0) to the bottom."
interface_depths(col::OceanColumnOptics) =
    [first(col.layers).depth_top; [l.depth_bottom for l in col.layers]]
"Depth of the column bottom `[m]`."
bottom_depth(col::OceanColumnOptics) = last(col.layers).depth_bottom
"Wavelength grid `[nm]`."
wavelengths(col::OceanColumnOptics) = first(col.layers).λ
n_wavelengths(col::OceanColumnOptics) = length(wavelengths(col))
n_moments(col::OceanColumnOptics) = first(col.layers).n_moments
is_polarized(col::OceanColumnOptics) = first(col.layers).polarized

"""
    bottom_albedo(col::OceanColumnOptics) -> Vector

Bottom albedo per wavelength: zeros for a [`BlackBottom`](@ref).
"""
bottom_albedo(col::OceanColumnOptics{FT, BlackBottom}) where {FT} = zeros(FT, n_wavelengths(col))
function bottom_albedo(col::OceanColumnOptics{FT, <:LambertianBottom}) where {FT}
    a = col.bottom.albedo
    return length(a) == 1 ? fill(FT(only(a)), n_wavelengths(col)) : FT.(a)
end

Base.length(col::OceanColumnOptics) = n_layers(col)
Base.getindex(col::OceanColumnOptics, i) = col.layers[i]
Base.iterate(col::OceanColumnOptics, st...) = iterate(col.layers, st...)
Base.firstindex(col::OceanColumnOptics) = 1
Base.lastindex(col::OceanColumnOptics) = n_layers(col)

function Base.show(io::IO, col::OceanColumnOptics{FT}) where {FT}
    print(io, "OceanColumnOptics{$FT}($(n_layers(col)) layers, 0–$(bottom_depth(col)) m, ",
          n_wavelengths(col), " λ, ", nameof(typeof(col.bottom)), ")")
end
