# =============================================================================
# OceanLayer — one homogeneous vertical slab
# =============================================================================
#
# One homogeneous slab. Stacks of contiguous layers are `OceanColumn`s
# (ocean_column.jl).
# =============================================================================

"""
$(TYPEDEF)

One vertically-homogeneous ocean layer.

# Fields
$(TYPEDFIELDS)
"""
struct OceanLayer{FT<:Real}
    "Depth of layer top [m], 0 at surface"
    depth_top::FT
    "Depth of layer bottom [m]"
    depth_bottom::FT
    "Elastic optical constituents active in this layer"
    constituents::Vector{AbstractOceanConstituent{FT}}
    "Inelastic processes (fluorescence, Raman) active in this layer"
    inelastic::Vector{AbstractOceanInelasticProcess{FT}}

    function OceanLayer{FT}(depth_top, depth_bottom, constituents,
                            inelastic) where {FT}
        dt, db = FT(depth_top), FT(depth_bottom)
        db > dt || throw(ArgumentError("depth_bottom ($db) must be > depth_top ($dt)"))
        dt ≥ zero(FT) || throw(ArgumentError("depth_top must be non-negative"))
        new{FT}(dt, db, constituents, inelastic)
    end
end

# User-friendly constructor — pick FT from the first constituent (or first
# fluorophore, if there are no elastic constituents).
function OceanLayer(depth_top, depth_bottom, constituents::AbstractVector;
                    inelastic::AbstractVector = AbstractOceanInelasticProcess{Float64}[])
    FT = if !isempty(constituents)
            _constituent_float_type(first(constituents))
         elseif !isempty(inelastic)
            _inelastic_float_type(first(inelastic))
         else
            Float64
         end
    cvec = convert(Vector{AbstractOceanConstituent{FT}}, constituents)
    fvec = convert(Vector{AbstractOceanInelasticProcess{FT}}, inelastic)
    return OceanLayer{FT}(FT(depth_top), FT(depth_bottom), cvec, fvec)
end

_constituent_float_type(::AbstractOceanConstituent{FT}) where {FT} = FT
_inelastic_float_type(::AbstractOceanInelasticProcess{FT}) where {FT} = FT

"""
    thickness(x)

Return the geometric thickness `[m]` of a layer, a column or the optics of a layer.
"""
thickness(layer::OceanLayer) = layer.depth_bottom - layer.depth_top

"""
    midpoint_depth(layer::OceanLayer)

Return the depth `[m]` of the middle of `layer`.
"""
midpoint_depth(layer::OceanLayer) = (layer.depth_top + layer.depth_bottom) / 2

function Base.show(io::IO, layer::OceanLayer{FT}) where {FT}
    nf = length(layer.inelastic)
    print(io, "OceanLayer{$FT}($(layer.depth_top)-$(layer.depth_bottom) m, ",
              "$(length(layer.constituents)) constituents",
              nf > 0 ? ", $nf inelastic processes" : "",
              ")")
end
