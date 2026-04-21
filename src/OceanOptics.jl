"""
    OceanOptics

Compute Inherent Optical Properties (IOPs) and phase-function expansions
of natural waters from biogeochemical state. Produces an
`OceanLayerOptics` bundle that is neutral to any particular RT solver —
the primary consumer is `vSmartMOM.jl`'s `OceanSurface`, but any
plane-parallel code (HydroLight-style, Monte Carlo, SOS) can consume
the same output.

# Phase 1 scope (this release)

This first release covers the **elastic** optical path end-to-end:

- `PureWater` with choice of absorption model (Pope-Fry, Smith-Baker,
  Mason-Cone-Fry, or a combined piecewise) and scattering model
  (Morel or Zhang-Hu)
- `FournierForandPhase` for ocean particles (Mobley 2002
  parameterization), `RayleighWaterPhase` for molecular water
- Generic Legendre-moment machinery (`phase_function_moments`) that
  works for any `AbstractOceanPhaseFunction`
- `layer_optics(layer, λ)` — the main user entry point producing a
  complete `OceanLayerOptics` bundle
"""
module OceanOptics

using DocStringExtensions
using FastGaussQuadrature: gausslegendre
using Interpolations: linear_interpolation, Flat
using LinearAlgebra
using QuadGK: quadgk

# -----------------------------------------------------------------------------
# Type spine
# -----------------------------------------------------------------------------
include("types/abstract_types.jl")
include("types/iop.jl")
include("types/layer_optics.jl")
include("types/ocean_layer.jl")
include("types/ocean_column.jl")

# -----------------------------------------------------------------------------
# I/O: reference-data loaders and interpolation
# -----------------------------------------------------------------------------
include("io/data_loaders.jl")

# -----------------------------------------------------------------------------
# Phase functions
# -----------------------------------------------------------------------------
include("phase/abstract_phase.jl")
include("phase/moments.jl")
include("phase/fournier_forand.jl")
include("phase/rayleigh_water.jl")
include("phase/henyey_greenstein.jl")
include("phase/petzold.jl")

# -----------------------------------------------------------------------------
# Materials
# -----------------------------------------------------------------------------
include("materials/pure_water.jl")
include("materials/cdom.jl")
include("materials/phytoplankton.jl")
include("materials/nap.jl")

# -----------------------------------------------------------------------------
# Inelastic processes (fluorescence, Raman)
# -----------------------------------------------------------------------------
include("inelastic/redistribution.jl")
include("inelastic/processes.jl")
include("inelastic/fluorescence.jl")
include("inelastic/raman.jl")

# -----------------------------------------------------------------------------
# Mixing & assembly
# -----------------------------------------------------------------------------
include("mixing/layer_assembly.jl")

# =============================================================================
# Public API
# =============================================================================

# Abstract hierarchy
export AbstractOceanConstituent, AbstractAbsorber, AbstractScatterer, AbstractAbsorbingScatterer
export AbstractOceanPhaseFunction
export AbstractOceanInelasticProcess, AbstractSpectralRedistribution

# Concrete materials
export PureWater
export AbstractPureWaterAbsorptionModel, SmithBaker1981, PopeFry1997, MasonConeFry2016, CombinedPopeFryMason
export AbstractPureWaterScatteringModel, Morel1974, ZhangHu2009
export CDOM, AbstractCDOMModel, BricaudMorelPrieur1981
export Phytoplankton, AbstractPhytoplanktonModel, Bricaud1998, Gordon1992
export NonAlgalParticles, AbstractNAPModel, BabinBricaud2003, FellMineral

# Phase functions
export FournierForandPhase, RayleighWaterPhase, HenyeyGreensteinPhase, PetzoldPhase

# Bundles
export IOP, OceanLayer, OceanColumn, OceanLayerOptics
export thickness, midpoint_depth, n_moments, n_wavelengths
export n_layers, bottom_depth, midpoint_depths, layer_thicknesses
export uniform_column, fell_column

# Core methods
export absorption, scattering, backscattering, attenuation, phase_function
export phase_function_value, phase_function_moments
export backscatter_fraction, asymmetry_parameter
export single_scattering_albedo
export layer_optics, mix_phase_moments

# Data loading (power-user API)
export SpectralTable, linterp, data_path, load_spectral_csv

# Inelastic processes (Phase 2)
export GaussianEmission, WavenumberShiftEmission, fwhm_to_sigma
export IsotropicFluorescence, CDOMFluorescence, WaterRaman
export chlorophyll_fluorescence
export excitation_absorption, excitation_range, emission, is_isotropic
export quantum_yield, inelastic_coefficient, redistribution
export raman_peak_wavelength

end # module
