# Inelastic processes

Fluorescence, water Raman, and CDOM fluorescence packaged as kernels
consumed by a two-pass RT solver (Fell 1997 §2.2.9).

## Spectral redistribution

```@docs
AbstractSpectralRedistribution
GaussianEmission
WavenumberShiftEmission
redistribution
fwhm_to_sigma
```

## Concrete processes

```@docs
IsotropicFluorescence
CDOMFluorescence
WaterRaman
chlorophyll_fluorescence
raman_peak_wavelength
```

## Trait interface

```@docs
AbstractOceanInelasticProcess
excitation_absorption
emission
excitation_range
is_isotropic
inelastic_coefficient
```
