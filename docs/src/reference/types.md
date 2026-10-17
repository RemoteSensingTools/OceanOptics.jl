# Types

Foundational type spine — the abstract hierarchy, the volumetric IOP
bundle, layers and columns, and the solver-facing optics. Phase-function
and inelastic-process types are documented on [Phase functions](@ref)
and [Inelastic processes](@ref); material models on [Materials](@ref).

## Module

```@docs
OceanOptics
```

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

## Solver output

```@docs
layer_optics
column_optics
OceanLayerOptics
ScatteringComponent
OceanColumnOptics
AbstractOceanBottom
BlackBottom
LambertianBottom
bottom_albedo
interface_depths
wavelengths
n_moments
n_wavelengths
is_polarized
optical_depth
absorption_optical_depth
scattering_optical_depth
single_scattering_albedo
```
