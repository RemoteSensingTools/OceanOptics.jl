"""
    OceanOptics

Inherent optical properties (IOPs) and phase-matrix expansions of natural
waters, for plane-parallel radiative-transfer solvers — primarily the
coupled atmosphere–ocean solver in vSmartMOM.jl, without depending on it.

- Constituents: `PureWater`, `CDOM`, `Phytoplankton`, `NonAlgalParticles`,
  and exact published values (`TabulatedAbsorber`, `TabulatedScatterer`).
- Phase functions and matrices: `RayleighWaterPhase`, `FournierForandPhase`,
  `HenyeyGreensteinPhase`, `PetzoldPhase`, `TabulatedGreekPhase`,
  `TabulatedMuellerPhase`, `WavelengthDependentPhase`, all expanded as
  `GreekCoefficients` in vSmartMOM's convention.
- `layer_optics` → `OceanLayerOptics`: absorption optical depth plus one
  scattering component per constituent; `column_optics` →
  `OceanColumnOptics` with the interface refractive index and the bottom.
- Inelastic process kernels (fluorescence, water Raman) for the solver.

Units, conventions and validation rules are on the *Consumer API* docs page.
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
include("types/greek.jl")
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
include("phase/angular.jl")
include("phase/fournier_forand.jl")
include("phase/rayleigh_water.jl")
include("phase/henyey_greenstein.jl")
include("phase/petzold.jl")
include("phase/tabulated_greek.jl")
include("phase/tabulated_mueller.jl")
include("phase/wavelength_phase.jl")

# -----------------------------------------------------------------------------
# Materials
# -----------------------------------------------------------------------------
include("materials/refractive_index.jl")
include("materials/pure_water.jl")
include("materials/cdom.jl")
include("materials/phytoplankton.jl")
include("materials/nap.jl")
include("materials/tabulated.jl")

# -----------------------------------------------------------------------------
# Inelastic processes (fluorescence, Raman)
# -----------------------------------------------------------------------------
include("inelastic/redistribution.jl")
include("inelastic/processes.jl")
include("inelastic/fluorescence.jl")
include("inelastic/raman.jl")

# -----------------------------------------------------------------------------
# Layer and column assembly
# -----------------------------------------------------------------------------
include("mixing/layer_assembly.jl")

# =============================================================================
# Public API
# =============================================================================

# Abstract hierarchy
export AbstractOceanConstituent, AbstractAbsorber, AbstractScatterer, AbstractAbsorbingScatterer
export AbstractOceanPhaseFunction, AbstractRefractiveIndexModel, AbstractOceanBottom
export AbstractOceanInelasticProcess, AbstractSpectralRedistribution

# Constituents
export PureWater
export AbstractPureWaterAbsorptionModel, IOCCG2018, PopeFry1997, Kou1993, SmithBaker1981
export AbstractPureWaterScatteringModel, Morel1974, ZhangHu2009
export CDOM, AbstractCDOMModel, BricaudMorelPrieur1981
export Phytoplankton, AbstractPhytoplanktonModel, Bricaud1998, Gordon1992
export NonAlgalParticles, AbstractNAPModel, BabinBricaud2003, FellMineral
export TabulatedAbsorber, TabulatedScatterer
export absorption, scattering, backscattering, attenuation, scattering_phase, is_total_absorption

# Refractive index of seawater
export FixedRefractiveIndex, QuanFry1995, refractive_index

# Phase functions and expansions
export FournierForandPhase, RayleighWaterPhase, HenyeyGreensteinPhase, PetzoldPhase
export TabulatedGreekPhase, TabulatedMuellerPhase, WavelengthDependentPhase
export GreekCoefficients, resize_moments, depolarized, mix_greek, evaluate_scattering_matrix
export phase_function_value, phase_function_moments, greek_expansion
export depolarizing_greek_expansion
export backscatter_fraction, asymmetry_parameter, is_polarizable, is_wavelength_dependent

# Layers, columns and their optics
export IOP, OceanLayer, OceanColumn
export thickness, midpoint_depth, n_layers, bottom_depth, midpoint_depths, layer_thicknesses
export uniform_column, fell_column
export ScatteringComponent, OceanLayerOptics, OceanColumnOptics, BlackBottom, LambertianBottom
export layer_optics, column_optics
export n_moments, n_wavelengths, wavelengths, is_polarized
export optical_depth, absorption_optical_depth, scattering_optical_depth
export single_scattering_albedo, bottom_albedo, interface_depths

# Data loading (power-user API)
export SpectralTable, linterp, data_path, load_spectral_csv

# Inelastic processes
export GaussianEmission, WavenumberShiftEmission, fwhm_to_sigma
export IsotropicFluorescence, CDOMFluorescence, WaterRaman
export chlorophyll_fluorescence
export excitation_absorption, excitation_range, emission, is_isotropic
export quantum_yield, inelastic_coefficient, radiance_source_kernel, redistribution
export raman_peak_wavelength

end # module
