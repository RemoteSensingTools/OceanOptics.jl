using Test
using OceanOptics
using ForwardDiff

@testset "OceanOptics.jl" begin
    include("test_iop.jl")
    include("test_spectral_table.jl")
    include("test_pure_water.jl")
    include("test_refractive_index.jl")
    include("test_cdom.jl")
    include("test_greek.jl")
    include("test_phase.jl")
    include("test_phytoplankton.jl")
    include("test_pigments.jl")
    include("test_nap.jl")
    include("test_redistribution.jl")
    include("test_fluorescence.jl")
    include("test_raman.jl")
    include("test_tabulated.jl")
    include("test_hydrolight_reference.jl")
    include("test_layer_optics.jl")
    include("test_absorbed_light.jl")
    include("test_column.jl")
    include("test_solver_contract.jl")
    include("test_parameters.jl")
    include("test_ad.jl")
    include("static_helpers.jl")
    include("test_inference.jl")
    include("test_float32.jl")
    include("test_device.jl")
    include("test_aqua.jl")
    include("test_explicit_imports.jl")
    include("test_jet.jl")
end
