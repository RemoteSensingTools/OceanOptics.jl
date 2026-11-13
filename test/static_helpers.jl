# Fixtures shared by the JET, inference and Float32 tests.

iop_constituents(::Type{FT}) where {FT} =
    (PureWater{FT}(),
     PureWater{FT}(absorption_model = PopeFry1997(), scattering_model = Morel1974()),
     CDOM{FT}(), Phytoplankton{FT}(Chl = FT(0.5)),
     Phytoplankton{FT}(Chl = FT(0.5), pigments = ProspectPigments(carotenoid_ratio = FT(0.3))),
     Phytoplankton{FT}(Chl = FT(0.5), pigments = Bidigare1990Pigments(chl_a = FT(0.5), psc = FT(0.3))),
     NonAlgalParticles{FT}(), NonAlgalParticles{FT, FellMineral}())

scalar_phases(::Type{FT}) where {FT} =
    (FournierForandPhase{FT}(backscatter_fraction = FT(0.018)),
     HenyeyGreensteinPhase{FT}(g = FT(0.7)), RayleighWaterPhase{FT}(), PetzoldPhase{FT}())

# Water, CDOM and phytoplankton with a pigment partition, chlorophyll
# fluorescence and water Raman, 0–5 m.
function inelastic_layer(::Type{FT}) where {FT}
    phyto = Phytoplankton{FT}(Chl = FT(0.5), pigments = ProspectPigments(carotenoid_ratio = FT(0.3)))
    fluo  = chlorophyll_fluorescence(phyto; quantum_yield = FT(0.02), carotenoid_transfer = FT(0.5))
    return OceanLayer(FT(0), FT(5), AbstractOceanConstituent{FT}[PureWater{FT}(), CDOM{FT}(), phyto];
                      inelastic = AbstractOceanInelasticProcess{FT}[fluo, WaterRaman{FT}()])
end
