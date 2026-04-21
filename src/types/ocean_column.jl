# =============================================================================
# OceanColumn — stack of contiguous OceanLayers, surface-first
# =============================================================================
#
# An `OceanColumn` is an ordered, contiguous, surface-anchored stack of
# `OceanLayer`s. "Contiguous" here means no gaps and no overlaps: layer
# `i`'s bottom equals layer `i+1`'s top, exactly. The surface anchor
# (`layers[1].depth_top == 0`) reflects the physical convention that
# vertical RT inside the water column is referenced to the air-water
# interface.
#
# Construction helpers
# --------------------
#
#   uniform_column(constituents; n_layers, depth)
#       N equal-thickness layers from the surface to `depth`.
#
#   fell_column(constituents; depth)
#       Fell (1997) §5.2.1 layering: 1 m in 0–10 m, 2 m in 10–20 m,
#       5 m below. Mirrors the grid used in Fell's North-Sea coupled
#       atmosphere-ocean tests.
#
# Both helpers take a `fluorophores` keyword that is broadcast to every
# layer (typical pattern when chlorophyll fluorescence or Raman is
# active throughout the water column).
# =============================================================================

"""
$(TYPEDEF)

Ordered stack of `OceanLayer{FT}` spanning from the sea surface
(`layers[1].depth_top == 0`) downward. Contiguity (no gaps, no
overlaps) is enforced at construction time.

# Fields
$(TYPEDFIELDS)
"""
struct OceanColumn{FT}
    "Surface-first, contiguous list of layers"
    layers::Vector{OceanLayer{FT}}

    function OceanColumn{FT}(layers::Vector{OceanLayer{FT}}) where {FT}
        _validate_column_contiguity(layers)
        new{FT}(layers)
    end
end

OceanColumn(layers::Vector{OceanLayer{FT}}) where {FT} = OceanColumn{FT}(layers)

function _validate_column_contiguity(layers::Vector{<:OceanLayer})
    isempty(layers) &&
        throw(ArgumentError("OceanColumn: need at least one layer"))
    layers[1].depth_top == 0 ||
        throw(ArgumentError("OceanColumn: first layer must start at z = 0 " *
                            "(got $(layers[1].depth_top))"))
    @inbounds for i in 1:length(layers)-1
        layers[i].depth_bottom == layers[i+1].depth_top ||
            throw(ArgumentError("OceanColumn: layers $i and $(i+1) are not " *
                                "contiguous: depth_bottom=$(layers[i].depth_bottom) " *
                                "≠ depth_top=$(layers[i+1].depth_top)"))
    end
    return nothing
end

# --- Basic accessors ---------------------------------------------------------

"Number of layers in the column."
n_layers(col::OceanColumn) = length(col.layers)

"Total geometric thickness `[m]` — surface to column bottom."
thickness(col::OceanColumn) = col.layers[end].depth_bottom - col.layers[1].depth_top

"Depth of the column bottom `[m]`."
bottom_depth(col::OceanColumn) = col.layers[end].depth_bottom

"Per-layer midpoint depths `[m]`."
midpoint_depths(col::OceanColumn) = [midpoint_depth(l) for l in col.layers]

"Per-layer thicknesses `[m]`."
layer_thicknesses(col::OceanColumn) = [thickness(l) for l in col.layers]

# --- Collection interface ---------------------------------------------------

Base.length(col::OceanColumn)          = n_layers(col)
Base.getindex(col::OceanColumn, i)     = col.layers[i]
Base.firstindex(col::OceanColumn)      = firstindex(col.layers)
Base.lastindex(col::OceanColumn)       = lastindex(col.layers)
Base.iterate(col::OceanColumn)         = iterate(col.layers)
Base.iterate(col::OceanColumn, st)     = iterate(col.layers, st)
Base.eltype(::Type{OceanColumn{FT}}) where {FT} = OceanLayer{FT}

function Base.show(io::IO, col::OceanColumn{FT}) where {FT}
    print(io, "OceanColumn{$FT}(", n_layers(col), " layers, 0–",
              bottom_depth(col), " m)")
end

# =============================================================================
# Constructor helpers
# =============================================================================

"""
    uniform_column(constituents; n_layers=10, depth=50.0, fluorophores=[])

Build an `OceanColumn` with `n_layers` equal-thickness layers spanning
0 → `depth` metres. Each layer carries the same `constituents` and
`fluorophores` lists.
"""
function uniform_column(constituents::AbstractVector;
                        n_layers::Integer = 10,
                        depth              = 50.0,
                        fluorophores::AbstractVector = AbstractOceanInelasticProcess{Float64}[])
    n_layers > 0 ||
        throw(ArgumentError("uniform_column: n_layers must be positive"))
    depth > 0 ||
        throw(ArgumentError("uniform_column: depth must be positive"))

    FT = _column_float_type(constituents, fluorophores)
    dz = FT(depth) / n_layers

    layers = Vector{OceanLayer{FT}}(undef, n_layers)
    for i in 1:n_layers
        z_top = (i - 1) * dz
        z_bot = i * dz
        layers[i] = OceanLayer(z_top, z_bot, constituents;
                               fluorophores = fluorophores)
    end
    return OceanColumn{FT}(layers)
end

"""
    fell_column(constituents; depth=60.0, fluorophores=[])

Build an `OceanColumn` following Fell (1997) §5.2.1: 1 m layers in
0–10 m, 2 m layers in 10–20 m, 5 m layers below. Truncates gracefully
if `depth < 20`; the deepest layer is trimmed so `depth_bottom ==
depth` exactly.
"""
function fell_column(constituents::AbstractVector;
                     depth::Real       = 60.0,
                     fluorophores::AbstractVector = AbstractOceanInelasticProcess{Float64}[])
    depth > 0 ||
        throw(ArgumentError("fell_column: depth must be positive"))

    FT = _column_float_type(constituents, fluorophores)
    D  = FT(depth)

    # Boundary list: 0 → 10 at Δ=1, 10 → 20 at Δ=2, 20 → D at Δ=5
    boundaries = FT[0]
    _append_staircase!(boundaries, D, FT(1),  FT(10))
    _append_staircase!(boundaries, D, FT(2),  FT(20))
    _append_staircase!(boundaries, D, FT(5),  D)

    # Ensure the final boundary is exactly D (if the staircase under-shot)
    last(boundaries) == D || push!(boundaries, D)

    layers = [OceanLayer(boundaries[i], boundaries[i+1], constituents;
                         fluorophores = fluorophores)
              for i in 1:length(boundaries)-1]
    return OceanColumn{FT}(layers)
end

"""
Append grid points `(last_boundary + Δ, last_boundary + 2Δ, …)` up to
`upper` (capped at the column bottom `D`). Caller invariant: the
previous stage's last boundary is already in `boundaries`.
"""
function _append_staircase!(boundaries::Vector{FT}, D::FT, Δ::FT, upper::FT) where {FT}
    last(boundaries) ≥ min(upper, D) && return boundaries
    next = last(boundaries) + Δ
    while next ≤ min(upper, D)
        push!(boundaries, next)
        next += Δ
    end
    return boundaries
end

# Infer FT from the constituents or fluorophores (first non-empty one wins).
function _column_float_type(constituents::AbstractVector,
                            fluorophores::AbstractVector)
    if !isempty(constituents)
        return _constituent_float_type(first(constituents))
    elseif !isempty(fluorophores)
        return _fluorophore_float_type(first(fluorophores))
    else
        return Float64
    end
end
