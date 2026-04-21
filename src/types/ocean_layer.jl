# =============================================================================
# OceanLayer — one homogeneous vertical slab
# =============================================================================
#
# Minimal Phase-1 implementation. An OceanColumn (stack with contiguity
# validation + Fell-grid helpers) will be added in Phase 2 when vSmartMOM
# needs to drive a stack through adding-doubling.
# =============================================================================

"""
$(TYPEDEF)

One vertically-homogeneous ocean layer.

# Fields
$(TYPEDFIELDS)
"""
struct OceanLayer{FT}
    "Depth of layer top [m], 0 at surface"
    depth_top::FT
    "Depth of layer bottom [m]"
    depth_bottom::FT
    "Elastic optical constituents active in this layer"
    constituents::Vector{<:AbstractOceanConstituent{FT}}
    "Inelastic processes (fluorescence, Raman) active in this layer"
    fluorophores::Vector{<:AbstractOceanInelasticProcess{FT}}

    function OceanLayer{FT}(depth_top, depth_bottom, constituents,
                            fluorophores) where {FT}
        dt, db = FT(depth_top), FT(depth_bottom)
        db > dt || throw(ArgumentError("depth_bottom ($db) must be > depth_top ($dt)"))
        dt ≥ zero(FT) || throw(ArgumentError("depth_top must be non-negative"))
        new{FT}(dt, db, constituents, fluorophores)
    end
end

# User-friendly constructor — pick FT from the first constituent (or first
# fluorophore, if there are no elastic constituents).
function OceanLayer(depth_top, depth_bottom, constituents::AbstractVector;
                    fluorophores::AbstractVector = AbstractOceanInelasticProcess{Float64}[])
    FT = if !isempty(constituents)
            _constituent_float_type(first(constituents))
         elseif !isempty(fluorophores)
            _fluorophore_float_type(first(fluorophores))
         else
            Float64
         end
    cvec = convert(Vector{AbstractOceanConstituent{FT}}, constituents)
    fvec = convert(Vector{AbstractOceanInelasticProcess{FT}}, fluorophores)
    return OceanLayer{FT}(FT(depth_top), FT(depth_bottom), cvec, fvec)
end

_constituent_float_type(::AbstractOceanConstituent{FT}) where {FT} = FT
_fluorophore_float_type(::AbstractOceanInelasticProcess{FT}) where {FT} = FT

"Geometric thickness [m]"
thickness(layer::OceanLayer) = layer.depth_bottom - layer.depth_top

"Midpoint depth [m]"
midpoint_depth(layer::OceanLayer) = (layer.depth_top + layer.depth_bottom) / 2

function Base.show(io::IO, layer::OceanLayer{FT}) where {FT}
    nf = length(layer.fluorophores)
    print(io, "OceanLayer{$FT}($(layer.depth_top)-$(layer.depth_bottom) m, ",
              "$(length(layer.constituents)) constituents",
              nf > 0 ? ", $nf fluorophores" : "",
              ")")
end
