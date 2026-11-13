# =============================================================================
# Layer and column optics — the solver-facing output
# =============================================================================
#
# A layer is described the way plane-parallel solvers combine optical
# properties: absorption components (one per absorber, or per absorption
# part such as a pigment group), scattering components (one per scatterer,
# with its phase-matrix expansion) and the inelastic processes active in
# the layer, each linked to the absorption components that excite it.
# Mixed quantities (total τ, ϖ, the mixed expansion) are derived on demand.
# Keeping components separate lets a solver build each component's phase
# operators once per band and truncate them individually, as it does for
# aerosol types; lets it report the light absorbed by each constituent
# (absorbed PAR, fluorescence excitation); and keeps the provenance of
# every phase matrix visible.
#
# The numeric arrays of every result struct have the array type `V` (and `I`
# for integer indices) as a type parameter, and `Adapt.adapt(to, x)` moves
# them to a device. Labels, symbols, content keys, the lists of components
# and layers, and the inelastic components stay on the host: the
# structs are containers, a kernel receives their arrays. `adapt` keeps the
# element type, keeps expansions that components share shared, and builds the
# components without validating them again. Validation and
# the per-wavelength derived quantities are reductions and broadcasts, so
# they run where the arrays live; functions that take a wavelength index
# `iλ` read single elements and are host operations.
# =============================================================================

# Neighbour comparisons as one broadcast over two views and a reduction.
_is_strictly_increasing(x::AbstractVector) =
    length(x) < 2 || all(@views x[begin+1:end] .> x[begin:end-1])
_is_strictly_monotone(x::AbstractVector) =
    _is_strictly_increasing(x) || all(@views x[begin+1:end] .< x[begin:end-1])

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

`phase_id[j]` is the content key of `greek[j]` (see [`phase_id`](@ref)):
entries with equal keys — in any component, layer, column or session —
hold equal coefficients, and entries that differ in any coefficient, in
`ℓ_max` or in element type have different keys, so a solver can build
phase operators once per distinct key.

`polarization` records how the phase-matrix coefficients were obtained:
`:scalar` (phase function only, `β`), `:full` (the constituent's phase
matrix), or `:depolarizing` (its phase function with the completely
depolarizing approximation, selected explicitly).
"""
struct ScatteringComponent{FT<:Real, V<:AbstractVector{FT}, I<:AbstractVector{Int}}
    "Constituent the component came from (type name, or the source of tabulated data)"
    label::String
    "Content key of each entry of `greek` (equal ⇔ equal coefficients, `ℓ_max` and element type)"
    phase_id::Vector{UInt64}
    "Scattering optical depth `b·Δz` per wavelength"
    τ::V
    "Distinct phase-matrix expansions"
    greek::Vector{GreekCoefficients{FT, V}}
    "Backscatter fraction of the phase function (not of its truncated expansion), per entry of `greek`"
    backscatter_fraction::V
    "Entry of `greek` used at each wavelength; empty when there is a single entry"
    phase_index::I
    "`:scalar`, `:full` or `:depolarizing`"
    polarization::Symbol

    function ScatteringComponent{FT, V, I}(label, τ, greek, backscatter_fraction, phase_index,
                                           polarization) where {FT, V<:AbstractVector{FT},
                                                                I<:AbstractVector{Int}}
        c = new{FT, V, I}(label, map(phase_id, greek), τ, greek, backscatter_fraction,
                          phase_index, polarization)
        _validate_component(c)
        return c
    end

    ScatteringComponent{FT, V, I}(::_Unchecked, label, ids, τ, greek, backscatter_fraction,
                                  phase_index,
                                  polarization) where {FT, V<:AbstractVector{FT},
                                                       I<:AbstractVector{Int}} =
        new{FT, V, I}(label, ids, τ, greek, backscatter_fraction, phase_index, polarization)
end

# The array types are those of `τ` and `phase_index`.
function ScatteringComponent{FT}(label, τ, greek, backscatter_fraction, phase_index,
                                 polarization) where {FT}
    t, index = _dense(FT, τ), _dense(Int, phase_index)
    return ScatteringComponent{FT, typeof(t), typeof(index)}(label, t, greek, backscatter_fraction,
                                                             index, polarization)
end

function ScatteringComponent(label::AbstractString, τ::AbstractVector{FT},
                             greek::AbstractVector{<:GreekCoefficients{FT}}, B::AbstractVector,
                             polarization::Symbol;
                             phase_index::AbstractVector{<:Integer} = Int[]) where {FT}
    return ScatteringComponent{FT}(String(label), _dense_copy(FT, τ), collect(greek),
                                   _dense_copy(FT, B), _dense_copy(Int, phase_index), polarization)
end

Adapt.adapt_structure(to, c::ScatteringComponent) = _adapt_scatterer(to, c, IdDict{Any, Any}())

# `shared` maps each host expansion to its adapted copy, so that expansions
# shared between components and layers are moved once and stay shared. The
# `phase_id` keys are those of the host arrays, which hold the same values.
function _adapt_scatterer(to, c::ScatteringComponent{FT}, shared::IdDict) where {FT}
    τ, index = adapt(to, c.τ), adapt(to, c.phase_index)
    _check_adapt_eltype(to, FT, τ)
    _check_adapt_eltype(to, Int, index)
    V     = typeof(τ)
    greek = GreekCoefficients{FT, V}[get!(() -> adapt(to, g), shared, g) for g in c.greek]
    return ScatteringComponent{FT, V, typeof(index)}(
        _Unchecked(), c.label, c.phase_id, τ, greek, adapt(to, c.backscatter_fraction), index,
        c.polarization)
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
        n = length(c.greek)        # the reduction captures a number, not the component
        all(i -> 1 ≤ i ≤ n, c.phase_index) ||
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

n_moments(c::ScatteringComponent) = n_moments(first(c.greek))
"""
    is_wavelength_dependent(c::ScatteringComponent)

Return whether the phase matrix of component `c` differs between wavelengths.
"""
is_wavelength_dependent(c::ScatteringComponent) = !isempty(c.phase_index)
is_polarized(c::ScatteringComponent) = c.polarization != :scalar

"""
    greek_expansion(c::ScatteringComponent, iλ) -> GreekCoefficients

Return the expansion of component `c` at the `iλ`-th wavelength of its layer.
"""
greek_expansion(c::ScatteringComponent, iλ::Integer) =
    is_wavelength_dependent(c) ? c.greek[c.phase_index[iλ]] : first(c.greek)

backscatter_fraction(c::ScatteringComponent, iλ::Integer) =
    is_wavelength_dependent(c) ? c.backscatter_fraction[c.phase_index[iλ]] :
                                 first(c.backscatter_fraction)

# Backscatter fraction at every wavelength of the layer: gathered through
# `phase_index`, or the single entry, which broadcasts against the grid.
_backscatter_on_grid(c::ScatteringComponent) =
    is_wavelength_dependent(c) ? c.backscatter_fraction[c.phase_index] : c.backscatter_fraction

"""
    phase_id(c::ScatteringComponent, iλ) -> key

Return the content key of the expansion of component `c` at the `iλ`-th wavelength of
its layer, equal to [`phase_id`](@ref)`(greek_expansion(c, iλ))`.
"""
phase_id(c::ScatteringComponent, iλ::Integer) =
    is_wavelength_dependent(c) ? c.phase_id[c.phase_index[iλ]] : first(c.phase_id)

# =============================================================================
# Absorption and inelastic components
# =============================================================================

"""
$(TYPEDEF)

Absorption by one constituent of a layer, or by one part of it (see
[`absorption_parts`](@ref)): `τ = a·Δz` per wavelength.

# Fields
$(TYPEDFIELDS)
"""
struct AbsorptionComponent{FT<:Real, V<:AbstractVector{FT}}
    "Constituent the component came from (type name, or the source of tabulated data)"
    label::String
    "Absorption part of the constituent (`:total` unless it splits its absorption)"
    part::Symbol
    "Absorption optical depth `a·Δz` per wavelength"
    τ::V

    function AbsorptionComponent{FT, V}(label, part, τ) where {FT, V<:AbstractVector{FT}}
        all(x -> isfinite(x) && x ≥ 0, τ) ||
            throw(ArgumentError("AbsorptionComponent \"$label\" (:$part): absorption optical " *
                                "depths must be finite and non-negative"))
        new{FT, V}(label, part, τ)
    end

    AbsorptionComponent{FT, V}(::_Unchecked, label, part, τ) where {FT, V<:AbstractVector{FT}} =
        new{FT, V}(label, part, τ)
end

function AbsorptionComponent{FT}(label, part, τ) where {FT}
    t = _dense(FT, τ)
    return AbsorptionComponent{FT, typeof(t)}(label, part, t)
end

AbsorptionComponent(label::AbstractString, part::Symbol, τ::AbstractVector{FT}) where {FT} =
    AbsorptionComponent{FT}(String(label), part, _dense_copy(FT, τ))

function Adapt.adapt_structure(to, c::AbsorptionComponent{FT}) where {FT}
    τ = adapt(to, c.τ)
    _check_adapt_eltype(to, FT, τ)
    return AbsorptionComponent{FT, typeof(τ)}(_Unchecked(), c.label, c.part, τ)
end

"""
$(TYPEDEF)

An inelastic process (fluorescence, water Raman) active in a layer, with
the absorption components of the layer that excite it.

# Fields
$(TYPEDFIELDS)

For fluorescence, `excitation[k]` indexes the layer's `absorbers` and
`efficiency[k]` is the energy-transfer efficiency η of that pathway, so
that the excitation coefficient on the layer grid is
`quantum_yield(process) · Σ_k efficiency[k] · absorbers[excitation[k]].τ / Δz`
— the same as [`excitation_absorption`](@ref)`(process, λ)`. Water Raman
has no excitation components.

The component is host data: `excitation` and `efficiency` have one entry per
excitation pathway and are read together, and `Adapt.adapt` returns the component
unchanged.
"""
struct InelasticComponent{FT<:Real, P<:AbstractOceanInelasticProcess}
    "Process type name"
    label::String
    "The process, with its kernels (excitation, emission, quantum yield)"
    process::P
    "Indices of the layer's absorption components that excite the process"
    excitation::Vector{Int}
    "Energy-transfer efficiency η of each excitation component"
    efficiency::Vector{FT}

    function InelasticComponent{FT, P}(label, process, excitation, efficiency) where {FT, P}
        length(excitation) == length(efficiency) ||
            throw(DimensionMismatch("InelasticComponent \"$label\": $(length(efficiency)) " *
                                    "efficiencies for $(length(excitation)) excitation components"))
        all(η -> isfinite(η) && 0 ≤ η ≤ 1, efficiency) ||
            throw(ArgumentError("InelasticComponent \"$label\": efficiencies must lie in [0, 1]"))
        new{FT, P}(label, process, excitation, efficiency)
    end
end

InelasticComponent{FT}(label, process::P, excitation,
                       efficiency) where {FT, P<:AbstractOceanInelasticProcess} =
    InelasticComponent{FT, P}(label, process, excitation, efficiency)

# The process is a model description; it stays on the host with the indices
# and efficiencies of its excitation pathways.
Adapt.adapt_structure(to, c::InelasticComponent) = c

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

[`absorption_optical_depth`](@ref) (all absorbers),
[`optical_depth`](@ref), [`scattering_optical_depth`](@ref),
[`single_scattering_albedo`](@ref), [`backscatter_fraction`](@ref) and
the mixed expansion [`greek_expansion`](@ref)`(o, iλ)` — the
scattering-weighted mixture of the scattering components. Light absorbed
per component, given the scalar irradiance: [`absorbed_power`](@ref),
[`absorbed_photons`](@ref), [`absorbed_par`](@ref).

# Validation

`λ` finite, positive and strictly monotone (increasing or decreasing);
`0 ≤ depth_top < depth_bottom`; optical depths finite and non-negative;
every component on the same wavelength grid; scattering components with
`n_moments` terms and the layer's polarization; inelastic components
referring to existing absorption components. Arrays are checked at
construction; mutating them afterwards is not re-checked.

# Devices

`λ` and the arrays of the components have the array type `V` (integer
indices `I`). `Adapt.adapt(to, o)` returns the layer with these arrays on
the device of `to`, for example `adapt(CuArray, o)`; labels, parts,
`phase_id` keys, the component lists and the inelastic components
stay on the host. The element type is kept: an adaptor that changes it
(`CuArray{Float32}` for a `Float64` layer) is refused; build the optics with
`layer_optics(Float32, …)` instead. Expansions shared by components stay
shared. On a device,
[`absorption_optical_depth`](@ref), [`scattering_optical_depth`](@ref),
[`optical_depth`](@ref), [`single_scattering_albedo`](@ref) and
[`backscatter_fraction`](@ref)`(o)` return device arrays. Functions of a
wavelength index (`greek_expansion(o, iλ)`) and the absorbed-light and
inelastic diagnostics are host operations: apply them before `adapt`.
"""
struct OceanLayerOptics{FT<:Real, V<:AbstractVector{FT}, I<:AbstractVector{Int}}
    "Wavelengths `[nm]`, strictly monotone"
    λ::V
    "Depth of the layer top below the surface `[m]`"
    depth_top::FT
    "Depth of the layer bottom below the surface `[m]`"
    depth_bottom::FT
    "Absorption components, one per absorber or absorption part"
    absorbers::Vector{AbsorptionComponent{FT, V}}
    "Scattering components, one per scattering constituent"
    scatterers::Vector{ScatteringComponent{FT, V, I}}
    "Inelastic processes, with the absorption components that excite them"
    inelastic::Vector{InelasticComponent{FT}}
    "Number of expansion terms, `ℓ_max + 1`"
    n_moments::Int
    "Whether every scattering component carries phase-matrix coefficients"
    polarized::Bool

    function OceanLayerOptics{FT, V, I}(λ, depth_top, depth_bottom, absorbers, scatterers,
                                        inelastic, n_moments,
                                        polarized) where {FT, V<:AbstractVector{FT},
                                                          I<:AbstractVector{Int}}
        o = new{FT, V, I}(λ, depth_top, depth_bottom, absorbers, scatterers, inelastic,
                          n_moments, polarized)
        _validate_layer(o)
        return o
    end
end

# The array type is that of `λ`; the index type that of the scattering
# components, `Vector{Int}` if their vector does not fix it.
function OceanLayerOptics{FT}(λ, depth_top, depth_bottom, absorbers, scatterers, inelastic,
                              n_moments, polarized) where {FT}
    x = _dense(FT, λ)
    return OceanLayerOptics{FT, typeof(x), _index_type(scatterers)}(
        x, depth_top, depth_bottom, absorbers, scatterers, inelastic, n_moments, polarized)
end

_index_type(::AbstractVector{ScatteringComponent{FT, V, I}}) where {FT, V, I} = I
_index_type(::AbstractVector) = Vector{Int}

Adapt.adapt_structure(to, o::OceanLayerOptics) = _adapt_layer(to, o, IdDict{Any, Any}())

function _adapt_layer(to, o::OceanLayerOptics{FT}, shared::IdDict) where {FT}
    λ     = adapt(to, o.λ)
    index = adapt(to, _phase_index_prototype(o))
    _check_adapt_eltype(to, FT, λ)
    _check_adapt_eltype(to, Int, index)
    V, I  = typeof(λ), typeof(index)
    return OceanLayerOptics{FT, V, I}(
        λ, o.depth_top, o.depth_bottom,
        AbsorptionComponent{FT, V}[adapt(to, c) for c in o.absorbers],
        ScatteringComponent{FT, V, I}[_adapt_scatterer(to, c, shared) for c in o.scatterers],
        o.inelastic, o.n_moments, o.polarized)
end

# An index array of the layer's index type, to find its type after `adapt`
# also when the layer has no scattering component.
_phase_index_prototype(o::OceanLayerOptics) =
    isempty(o.scatterers) ? Int[] : first(o.scatterers).phase_index

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
    for c in o.absorbers
        length(c.τ) == nλ ||
            throw(DimensionMismatch("$who: absorption component \"$(c.label)\" (:$(c.part)) has " *
                                    "$(length(c.τ)) optical depths for $nλ wavelengths"))
    end
    for c in o.inelastic
        all(i -> 1 ≤ i ≤ length(o.absorbers), c.excitation) ||
            throw(ArgumentError("$who: inelastic component \"$(c.label)\" refers to a missing " *
                                "absorption component"))
    end
    o.n_moments ≥ 1 || throw(ArgumentError("$who: n_moments must be positive"))
    for c in o.scatterers
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

thickness(o::OceanLayerOptics) = o.depth_bottom - o.depth_top
"""
    n_wavelengths(x)

Return the number of wavelengths of the grid of a layer's or a column's optics.
"""
n_wavelengths(o::OceanLayerOptics) = length(o.λ)
n_moments(o::OceanLayerOptics) = o.n_moments

"""
    is_polarized(x) -> Bool

Return whether `x` (layer, column, component or [`GreekCoefficients`](@ref)) carries
phase-matrix coefficients in addition to the phase function.
"""
is_polarized(o::OceanLayerOptics) = o.polarized

"""
    absorber_index(o::OceanLayerOptics, label, part = :total) -> Int

Return the index in `o.absorbers` of the absorption component with `label` and `part`,
e.g. `absorber_index(o, "Phytoplankton", :chlorophyll)`.

The same index names the rows of [`absorbed_power`](@ref) and the entries of
[`absorbed_par`](@ref). Throws if there is no such component or more than one.
"""
function absorber_index(o::OceanLayerOptics, label::AbstractString, part::Symbol = :total)
    i = findall(c -> c.label == label && c.part === part, o.absorbers)
    length(i) == 1 && return only(i)
    throw(ArgumentError(isempty(i) ?
        "absorber_index: no absorption component \"$label\" (:$part); the layer has " *
        join(("\"$(c.label)\" (:$(c.part))" for c in o.absorbers), ", ") :
        "absorber_index: $(length(i)) absorption components \"$label\" (:$part); " *
        "give equal constituents distinct sources or index o.absorbers directly"))
end

# Sum of the components' optical depths, in an array like `λ`.
function _summed_optical_depth(λ::AbstractVector, components)
    τ = zero(λ)
    for c in components
        τ .+= c.τ
    end
    return τ
end

"""
    absorption_optical_depth(o::OceanLayerOptics)

Return the absorption optical depth of the layer per wavelength, summed over its
absorption components.
"""
absorption_optical_depth(o::OceanLayerOptics) = _summed_optical_depth(o.λ, o.absorbers)

"""
    scattering_optical_depth(o::OceanLayerOptics)

Return the scattering optical depth of the layer per wavelength, summed over its
scattering components.
"""
scattering_optical_depth(o::OceanLayerOptics) = _summed_optical_depth(o.λ, o.scatterers)

"""
    optical_depth(o::OceanLayerOptics)

Return the extinction optical depth `τ_abs + τ_scat` of the layer per wavelength.
"""
optical_depth(o::OceanLayerOptics) = absorption_optical_depth(o) .+ scattering_optical_depth(o)

"""
    single_scattering_albedo(o::OceanLayerOptics) -> AbstractVector

Return the single-scattering albedo `ϖ = τ_scat / (τ_abs + τ_scat)` of the layer per
wavelength, zero where `τ = 0`.
"""
function single_scattering_albedo(o::OceanLayerOptics{FT}) where {FT}
    τs = scattering_optical_depth(o)
    τ  = absorption_optical_depth(o) .+ τs
    return ifelse.(τ .> 0, τs ./ τ, zero(FT))
end

"""
    backscatter_fraction(o::OceanLayerOptics) -> AbstractVector

Return the scattering-weighted backscatter fraction of the components' phase functions
per wavelength, zero where the layer does not scatter.
"""
function backscatter_fraction(o::OceanLayerOptics{FT}) where {FT}
    τs = scattering_optical_depth(o)
    B  = zero(o.λ)
    for c in o.scatterers
        B .+= c.τ .* _backscatter_on_grid(c)
    end
    return ifelse.(τs .> 0, B ./ τs, zero(FT))
end

"""
    greek_expansion(o::OceanLayerOptics, iλ) -> GreekCoefficients

Return the expansion of the layer's mixed phase matrix at its `iλ`-th wavelength: the
components weighted by their scattering optical depths.

Where the layer does not scatter, the result is the isotropic expansion `β = [1, 0, …]`
with zero polarization terms.
"""
function greek_expansion(o::OceanLayerOptics{FT}, iλ::Integer) where {FT}
    w = [c.τ[iλ] for c in o.scatterers]
    if isempty(w) || sum(w) == 0
        # `fill`, not `zeros(FT, n)`: for the abstract `FT` of a not fully typed layer
        # inference (and JET) cannot exclude the `zeros` method that returns a matrix.
        β = fill(zero(FT), o.n_moments); β[1] = one(FT)
        g = GreekCoefficients(β)
        return o.polarized ? depolarized(g) : g
    end
    return mix_greek(w, [greek_expansion(c, iλ) for c in o.scatterers])
end

function Base.show(io::IO, o::OceanLayerOptics{FT}) where {FT}
    print(io, "OceanLayerOptics{$FT}($(o.depth_top)–$(o.depth_bottom) m, ",
          n_wavelengths(o), " λ, ", length(o.absorbers), " absorption, ",
          length(o.scatterers), " scattering",
          isempty(o.inelastic) ? "" : ", $(length(o.inelastic)) inelastic", " components, ",
          "ℓ ≤ $(o.n_moments - 1), ", o.polarized ? "polarized" : "scalar", ")")
end

# =============================================================================
# Bottom boundary descriptors
# =============================================================================

"""
$(TYPEDEF)

Lower boundary of an ocean column, for the solver to apply.
"""
abstract type AbstractOceanBottom end

"""
$(TYPEDEF)

Perfectly absorbing bottom.
"""
struct BlackBottom <: AbstractOceanBottom end

"""
$(TYPEDEF)

Lambertian (isotropically reflecting, depolarizing) bottom with albedo per
wavelength (or one value for all wavelengths), in `[0, 1]`.

# Fields
$(TYPEDFIELDS)
"""
struct LambertianBottom{FT<:Real, V<:AbstractVector{FT}} <: AbstractOceanBottom
    "Albedo: one value, or one per wavelength of the column"
    albedo::V

    function LambertianBottom{FT, V}(albedo) where {FT, V<:AbstractVector{FT}}
        isempty(albedo) && throw(ArgumentError("LambertianBottom: no albedo"))
        all(a -> isfinite(a) && 0 ≤ a ≤ 1, albedo) ||
            throw(ArgumentError("LambertianBottom: albedo must lie in [0, 1]"))
        new{FT, V}(albedo)
    end

    LambertianBottom{FT, V}(::_Unchecked, albedo) where {FT, V<:AbstractVector{FT}} =
        new{FT, V}(albedo)
end

function LambertianBottom{FT}(albedo) where {FT}
    a = _dense(FT, albedo)
    return LambertianBottom{FT, typeof(a)}(a)
end

LambertianBottom(albedo::Real) = LambertianBottom{float(typeof(albedo))}([float(albedo)])
LambertianBottom(albedo::AbstractVector{<:Real}) =
    LambertianBottom{float(eltype(albedo))}(_dense_copy(float(eltype(albedo)), albedo))

function Adapt.adapt_structure(to, b::LambertianBottom{FT}) where {FT}
    albedo = adapt(to, b.albedo)
    _check_adapt_eltype(to, FT, albedo)
    return LambertianBottom{FT, typeof(albedo)}(_Unchecked(), albedo)
end

"""
$(TYPEDEF)

Lambertian bottom whose albedo spectrum is tabulated on its own wavelength
grid, e.g. a measured sediment, seagrass or coral reflectance.
[`column_optics`](@ref) interpolates it linearly onto the column's
wavelengths (which must lie within the table) and stores the result as a
[`LambertianBottom`](@ref).

# Fields
$(TYPEDFIELDS)

# Construction

    TabulatedLambertianBottom(λ, albedo; source = "")
"""
struct TabulatedLambertianBottom{FT<:Real} <: AbstractOceanBottom
    "Wavelengths `[nm]`, strictly increasing"
    λ::Vector{FT}
    "Albedo at each `λ`, in `[0, 1]`"
    albedo::Vector{FT}
    "Provenance"
    source::String

    function TabulatedLambertianBottom{FT}(λ, albedo, source) where {FT}
        who = "TabulatedLambertianBottom"
        length(λ) ≥ 2 || throw(ArgumentError("$who: need at least two wavelengths"))
        all(x -> isfinite(x) && x > 0, λ) && _is_strictly_increasing(λ) ||
            throw(ArgumentError("$who: wavelengths must be finite, positive and strictly increasing"))
        length(albedo) == length(λ) ||
            throw(DimensionMismatch("$who: $(length(albedo)) albedos for $(length(λ)) wavelengths"))
        all(a -> isfinite(a) && 0 ≤ a ≤ 1, albedo) ||
            throw(ArgumentError("$who: albedo must lie in [0, 1]"))
        new{FT}(λ, albedo, source)
    end
end

function TabulatedLambertianBottom(λ::AbstractVector{<:Real}, albedo::AbstractVector{<:Real};
                                   source::AbstractString = "")
    FT = float(promote_type(eltype(λ), eltype(albedo)))
    return TabulatedLambertianBottom{FT}(collect(FT, λ), collect(FT, albedo), String(source))
end

# =============================================================================
# Column
# =============================================================================

"""
$(TYPEDEF)

Optical description of a whole water column, produced by
[`column_optics`](@ref): contiguous layers from the surface down on one
wavelength grid, the refractive index of the water at the surface relative
to air, and the bottom: `BlackBottom`, or `LambertianBottom` with one
albedo per wavelength or one for all (a [`TabulatedLambertianBottom`](@ref)
given to [`column_optics`](@ref) is resolved into the latter). Other
`AbstractOceanBottom` subtypes pass through for solvers that know them.

# Fields
$(TYPEDFIELDS)

Layers share `λ`, `n_moments` and polarization; the first starts at the
surface and each starts where the previous one ends. `refractive_index`
is real and given per wavelength.

# Devices

`Adapt.adapt(to, col)` returns the column with every numeric array (the
layers' arrays, `refractive_index`, the albedo of a `LambertianBottom`)
on the device of `to`, in the same element type; see
[`OceanLayerOptics`](@ref) for what stays on the host. Expansions that
layers share on the host are moved once and stay shared on the device.
"""
struct OceanColumnOptics{FT<:Real, B<:AbstractOceanBottom, V<:AbstractVector{FT},
                         I<:AbstractVector{Int}}
    "Layers, surface first"
    layers::Vector{OceanLayerOptics{FT, V, I}}
    "Refractive index of the surface water relative to air, per wavelength"
    refractive_index::V
    "Lower boundary"
    bottom::B

    function OceanColumnOptics{FT, B, V, I}(layers, refractive_index,
                                            bottom) where {FT, B, V<:AbstractVector{FT},
                                                           I<:AbstractVector{Int}}
        col = new{FT, B, V, I}(layers, refractive_index, bottom)
        _validate_column(col)
        return col
    end
end

# The array types are those of the layers.
OceanColumnOptics{FT, B}(layers::AbstractVector{OceanLayerOptics{FT, V, I}}, refractive_index,
                         bottom) where {FT, B, V, I} =
    OceanColumnOptics{FT, B, V, I}(layers, refractive_index, bottom)

# The element type of `layers` need not fix the array types: the first layer does.
function OceanColumnOptics(layers::AbstractVector{<:OceanLayerOptics}, n::AbstractVector,
                           bottom::AbstractOceanBottom)
    isempty(layers) && throw(ArgumentError("OceanColumnOptics: no layers"))
    return _column_optics(convert(Vector{typeof(first(layers))}, layers), n, bottom)
end

_column_optics(layers::Vector{OceanLayerOptics{FT, V, I}}, n, bottom::B) where {FT, V, I, B} =
    OceanColumnOptics{FT, B, V, I}(layers, _dense_copy(FT, n), bottom)

function Adapt.adapt_structure(to, col::OceanColumnOptics{FT}) where {FT}
    n = adapt(to, col.refractive_index)
    _check_adapt_eltype(to, FT, n)
    shared = IdDict{Any, Any}()
    layers = map(l -> _adapt_layer(to, l, shared), col.layers)
    bottom = adapt(to, col.bottom)
    return OceanColumnOptics{FT, typeof(bottom)}(layers, n, bottom)
end

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
    all(n -> isfinite(n) && n ≥ 1, col.refractive_index) ||
        throw(ArgumentError("$who: the refractive index of water relative to air must be " *
                            "finite and ≥ 1; check that it is n_water / n_air and not its " *
                            "inverse"))
    _check_bottom(col.bottom, nλ)
    return nothing
end

# What a column requires of its bottom on a grid of `nλ` wavelengths.
_check_bottom(::AbstractOceanBottom, nλ) = nothing
_check_bottom(b::LambertianBottom, nλ) =
    length(b.albedo) in (1, nλ) ||
        throw(DimensionMismatch("OceanColumnOptics: bottom albedo needs 1 or $nλ values"))
_check_bottom(::TabulatedLambertianBottom, nλ) =
    throw(ArgumentError("OceanColumnOptics: a TabulatedLambertianBottom must be interpolated " *
                        "onto the column grid; pass it to column_optics"))

n_layers(col::OceanColumnOptics) = length(col.layers)
"""
    interface_depths(col::OceanColumnOptics)

Return the depths `[m]` of the layer interfaces, from the surface (0) to the bottom.
"""
interface_depths(col::OceanColumnOptics) =
    [first(col.layers).depth_top; [l.depth_bottom for l in col.layers]]
bottom_depth(col::OceanColumnOptics) = last(col.layers).depth_bottom
"""
    wavelengths(col::OceanColumnOptics)

Return the wavelength grid `[nm]` of the column.
"""
wavelengths(col::OceanColumnOptics) = first(col.layers).λ
n_wavelengths(col::OceanColumnOptics) = length(wavelengths(col))
n_moments(col::OceanColumnOptics) = first(col.layers).n_moments
is_polarized(col::OceanColumnOptics) = first(col.layers).polarized

"""
    bottom_albedo(col::OceanColumnOptics) -> AbstractVector

Return the bottom albedo per wavelength: zeros for a [`BlackBottom`](@ref), the stored
albedo for a [`LambertianBottom`](@ref).

Any other bottom type has no Lambertian albedo and throws an `ArgumentError`.
"""
bottom_albedo(col::OceanColumnOptics{FT, BlackBottom}) where {FT<:Real} = zero(wavelengths(col))
bottom_albedo(col::OceanColumnOptics) =
    throw(ArgumentError("bottom_albedo: a $(typeof(col.bottom)) has no Lambertian albedo; " *
                        "only BlackBottom and LambertianBottom do. Check `col.bottom` and " *
                        "treat this bottom type in the solver, or define " *
                        "`bottom_albedo(::OceanColumnOptics{FT, <:$(nameof(typeof(col.bottom)))})`"))
# One albedo broadcasts against the grid; one per wavelength is copied.
bottom_albedo(col::OceanColumnOptics{FT, <:LambertianBottom}) where {FT<:Real} =
    FT.(col.bottom.albedo) .+ zero(wavelengths(col))

Base.length(col::OceanColumnOptics) = n_layers(col)
Base.getindex(col::OceanColumnOptics, i) = col.layers[i]
Base.iterate(col::OceanColumnOptics, st...) = iterate(col.layers, st...)
Base.firstindex(col::OceanColumnOptics) = 1
Base.lastindex(col::OceanColumnOptics) = n_layers(col)

function Base.show(io::IO, col::OceanColumnOptics{FT}) where {FT}
    print(io, "OceanColumnOptics{$FT}($(n_layers(col)) layers, 0–$(bottom_depth(col)) m, ",
          n_wavelengths(col), " λ, ", nameof(typeof(col.bottom)), ")")
end
