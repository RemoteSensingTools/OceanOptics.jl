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

@testset "Integer and rational element types" begin
    # Arrays of integers or rationals become floating-point IOPs, as scalars do.
    iop = IOP([1, 2], [3, 4], [0, 1])
    @test iop isa IOP{Float64} && iop.a == [1.0, 2.0] && iop.bb == [0.0, 1.0]
    @test IOP([1//2], [1//4], [0//1]) isa IOP{Float64}
    @test IOP(1, 2, 0) isa IOP{Float64}
    # Models and optics are not built in integer or rational arithmetic.
    for T in (Int, Rational{Int})
        @test_throws "floating-point type" CDOM{T}()
        @test_throws "floating-point type" PureWater{T}()
        @test_throws "floating-point type" Phytoplankton{T}(Chl = 1)
        @test_throws "floating-point type" NonAlgalParticles{T}()
        @test_throws "floating-point type" NonAlgalParticles{T, FellMineral}()
        @test_throws "floating-point type" GordonMorel1983{T}()
        @test_throws "floating-point type" HenyeyGreensteinPhase{T}()
        @test_throws "floating-point type" RayleighWaterPhase{T}()
        @test_throws "floating-point type" FournierForandPhase{T}(backscatter_fraction = 0.01)
    end
    layer = OceanLayer(0.0, 1.0, [PureWater(), CDOM()])
    @test_throws "floating-point type" layer_optics(Int, layer, [500.0])
    @test_throws "floating-point type" layer_optics(Rational{Int}, layer, [500.0])
    col = OceanColumn([layer])
    @test_throws "floating-point type" column_optics(Rational{Int}, col, [500.0];
                                                     refractive_index = FixedRefractiveIndex(1.34),
                                                     bottom = BlackBottom())
end
