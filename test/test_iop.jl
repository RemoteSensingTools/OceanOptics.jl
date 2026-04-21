@testset "IOP" begin
    # Scalar construction and accessors
    iop = IOP(0.1, 0.05, 0.018)
    @test iop.a  == 0.1
    @test iop.b  == 0.05
    @test iop.bb == 0.018

    # Derived quantities
    @test attenuation(iop) ≈ 0.15
    @test single_scattering_albedo(iop) ≈ 0.05 / 0.15
    @test backscatter_fraction(iop)    ≈ 0.018 / 0.05

    # Degenerate case: no absorption, no scattering
    zero_iop = IOP(0.0, 0.0, 0.0)
    @test single_scattering_albedo(zero_iop) == 0.0
    @test backscatter_fraction(zero_iop)    == 0.0

    # Additive mixing of independent constituents
    water  = IOP(0.01, 0.002, 0.001)
    cdom   = IOP(0.1,  0.0,   0.0)
    phyto  = IOP(0.02, 0.05,  0.0005)
    total  = water + cdom + phyto
    @test total.a  ≈ 0.13
    @test total.b  ≈ 0.052
    @test total.bb ≈ 0.0015

    # Scalar multiplication (e.g. for scaling by layer thickness)
    @test (2.0 * iop).a ≈ 0.2
    @test (iop * 3.0).b ≈ 0.15

    # zero() on an instance
    z = zero(iop)
    @test z.a == 0.0 && z.b == 0.0 && z.bb == 0.0

    # Vector IOP for a wavelength grid
    λ_iop = IOP([0.01, 0.02], [0.005, 0.006], [0.002, 0.0025])
    @test λ_iop.a isa Vector{Float64}
    @test (λ_iop + λ_iop).a == [0.02, 0.04]
    @test single_scattering_albedo(λ_iop) ≈ [0.005/0.015, 0.006/0.026]

    # Mixed-type promotion
    promoted = IOP(0.1, 0.05f0, 1//10)
    @test promoted.a isa Float64   # promoted to Float64 since 0.1 is Float64
end
