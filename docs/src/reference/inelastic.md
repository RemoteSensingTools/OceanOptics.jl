# Inelastic processes

Fluorescence, water Raman, and CDOM fluorescence packaged as kernels
consumed by a two-pass RT solver (Fell 1997 §2.2.9), and the light
absorbed per layer and component that drives them.

## Spectral redistribution

```@docs
AbstractSpectralRedistribution
TabulatedEmission
chloroplast_emission
GaussianEmission
WavenumberShiftEmission
redistribution
fwhm_to_sigma
```

## Concrete processes

```@docs
ExcitationPathway
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
excitation_pathways
is_isotropic
inelastic_coefficient
radiance_source_kernel
```

## Light absorbed per layer

```@docs
absorbed_power
absorbed_photons
absorbed_par
inelastic_emission
excitation_optical_depth
absorber_index
```
