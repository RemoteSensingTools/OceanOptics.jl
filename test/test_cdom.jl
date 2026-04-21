@testset "CDOM" begin
    cdom = CDOM{Float64}(a_ref = 0.03, λ_ref = 440.0, slope = 0.014)

    # Anchor point: at λ_ref, a_g equals a_ref exactly
    @test absorption(cdom, 440.0) == 0.03

    # Exponential decay: a(λ+1/S) = a(λ)/e
    @test absorption(cdom, 440.0 + 1/0.014) ≈ 0.03 / ℯ rtol=1e-10

    # Monotonic decay into the red
    @test absorption(cdom, 500.0) < absorption(cdom, 440.0)
    @test absorption(cdom, 400.0) > absorption(cdom, 440.0)

    # CDOM is a pure absorber: scattering/backscattering default to 0
    @test scattering(cdom, 500.0)     == 0.0
    @test backscattering(cdom, 500.0) == 0.0
    @test attenuation(cdom, 500.0)    == absorption(cdom, 500.0)

    # Non-negativity
    @test absorption(cdom, 400.0) > 0
    @test absorption(cdom, 700.0) > 0

    # Vectorized
    λ = [400.0, 440.0, 500.0, 600.0, 700.0]
    a = absorption(cdom, λ)
    @test length(a) == 5
    @test a[2] == 0.03
    @test issorted(a, rev=true)

    # Constructor edge cases
    @test_throws ArgumentError CDOM{Float64}(a_ref = -0.01)
    @test_throws ArgumentError CDOM{Float64}(λ_ref = -1.0)

    # Float32 parametric
    cdom32 = CDOM{Float32}(a_ref = 0.03f0, slope = 0.014f0)
    @test absorption(cdom32, 440.0f0) isa Float32
end
