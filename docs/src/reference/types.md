# Types

Foundational type spine — the abstract hierarchy, the volumetric IOP
bundle, and the layer / column containers. Phase-function and
inelastic-process abstract types are cross-referenced from
[Phase functions](@ref) and [Inelastic](@ref); material-model tags
from [Materials](@ref).

## Abstract constituent hierarchy

```@docs
AbstractOceanConstituent
AbstractAbsorber
AbstractScatterer
AbstractAbsorbingScatterer
```

## Volumetric IOPs

```@docs
IOP
single_scattering_albedo
```

## Layers and columns

```@docs
OceanLayer
OceanColumn
thickness
midpoint_depth
n_layers
bottom_depth
midpoint_depths
layer_thicknesses
uniform_column
fell_column
```

## RT-solver output

```@docs
OceanLayerOptics
n_moments
n_wavelengths
layer_optics
```
