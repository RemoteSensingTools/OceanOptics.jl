# Phase functions

## Concrete phase functions and matrices

```@docs
RayleighWaterPhase
FournierForandPhase
HenyeyGreensteinPhase
PetzoldPhase
TabulatedGreekPhase
TabulatedMuellerPhase
WavelengthDependentPhase
```

## Expansions

```@docs
GreekCoefficients
resize_moments
depolarized
mix_greek
evaluate_scattering_matrix
```

## Shared interface

```@docs
AbstractOceanPhaseFunction
phase_function_value
phase_function_moments
greek_expansion
depolarizing_greek_expansion
backscatter_fraction
asymmetry_parameter
is_polarizable
is_wavelength_dependent
```
