# The device tests of the package on a CUDA device. CUDA is a dependency of
# this environment only; see README.md in this directory.

using CUDA: CUDA, CuArray
using OceanOptics
using Test

if CUDA.functional()
    CUDA.allowscalar(false)
    const DEVICE_ARRAY = CuArray
    @info "OceanOptics device tests on $(CUDA.name(CUDA.device()))"
    @testset "OceanOptics.jl on CUDA" begin
        include(joinpath(@__DIR__, "..", "test_device.jl"))
    end
else
    @info "CUDA is not functional on this host; the CUDA device tests are skipped"
end
