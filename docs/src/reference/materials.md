# Materials

## Pure (sea)water

```@docs
PureWater
AbstractPureWaterAbsorptionModel
AbstractPureWaterScatteringModel
IOCCG2018
PopeFry1997
Kou1993
SmithBaker1981
Morel1974
ZhangHu2009
```

## Refractive index of seawater

```@docs
AbstractRefractiveIndexModel
FixedRefractiveIndex
QuanFry1995
refractive_index
```

## CDOM

```@docs
CDOM
AbstractCDOMModel
BricaudMorelPrieur1981
```

## Phytoplankton

```@docs
Phytoplankton
AbstractPhytoplanktonModel
Bricaud1998
Gordon1992
```

## Non-algal particles

```@docs
NonAlgalParticles
AbstractNAPModel
BabinBricaud2003
FellMineral
```

## Tabulated constituents

Exact published coefficients, without an empirical model; see
[Consumer API](@ref) for the absorption-scope rule.

```@docs
TabulatedAbsorber
TabulatedScatterer
is_total_absorption
```

## Shared traits

```@docs
absorption
scattering
backscattering
attenuation
scattering_phase
```
