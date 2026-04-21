# Consumer API

OceanOptics.jl publishes a stable, RT-agnostic interface. Downstream
solvers (vSmartMOM, HydroLight-style ODE codes, Monte-Carlo, SOS, …)
consume it *without* OceanOptics having any knowledge of the solver.
`Project.toml` therefore lists no RT-solver package under `[deps]` or
`[weakdeps]`; adapters live on the solver side.

## The stable contract

A downstream solver should consume only the symbols below. Anything
else in the package is internal and can change in any minor version.

### Output bundle

See [`OceanLayerOptics`](@ref) and [`layer_optics`](@ref).

### Constituent hierarchy

[`AbstractOceanConstituent`](@ref), [`AbstractAbsorber`](@ref),
[`AbstractScatterer`](@ref), [`AbstractAbsorbingScatterer`](@ref), and
the trait functions [`absorption`](@ref), [`scattering`](@ref),
[`backscattering`](@ref), [`attenuation`](@ref),
[`phase_function`](@ref). Canonical docstrings in [Materials](@ref).

### Inelastic-process hierarchy

[`AbstractOceanInelasticProcess`](@ref) and trait functions
[`excitation_absorption`](@ref), [`emission`](@ref),
[`excitation_range`](@ref), [`is_isotropic`](@ref),
[`inelastic_coefficient`](@ref). Canonical docstrings in
[Inelastic](@ref).

### Layer / column

[`OceanLayer`](@ref) and [`OceanColumn`](@ref). Canonical docstrings
in [Types](@ref).

## The two-pass solve convention (Fell 1997 §2.2.9)

Inelastic sources require two passes:

1. **Pass 1** — solve the elastic RT problem using `OceanLayerOptics` to
   get the scalar irradiance field `E°(τ, λ')`.

2. **Pass 2** — for each inelastic process `p` in
   `layer.fluorophores`, evaluate the source kernel

   ```math
   J^I(\tau, \lambda) = \int_{\lambda'} b^I(\lambda', \lambda) \,
                                       \langle L(\tau, \lambda') \rangle_\Omega \,
                                       p^I(\cos\theta; \lambda', \lambda) \, d\lambda'
   ```

   with `b^I(λ', λ) = excitation_absorption(p, λ') × emission(p, λ', λ)`.
   For isotropic processes (`is_isotropic(p) == true`) the angular factor
   `p^I` reduces to a scalar and the source contributes only to the
   `m = 0` Fourier mode.

The adapter code that packs these kernels into a solver's internal
format (e.g. vSmartMOM's `OceanRS <: AbstractRamanType`) lives in the
solver package.
