@testset "NonAlgalParticles" begin

    @testset "Babin-Bricaud (2003) detritus + minerals" begin
        nap = NonAlgalParticles{Float64, BabinBricaud2003}(
            a_ref = 0.05, b_ref = 0.45, λ_ref = 443.0, slope = 0.011)

        @test absorption(nap, 443.0) ≈ 0.05         # anchor
        @test absorption(nap, 443.0 + 1/0.011) ≈ 0.05/ℯ rtol=1e-10
        @test absorption(nap, 500.0) < absorption(nap, 443.0)

        # Scattering: b(λ) = b_ref * (λ_ref / λ)
        @test scattering(nap, 443.0) ≈ 0.45
        @test scattering(nap, 886.0) ≈ 0.225 rtol=1e-12

        # Backscatter scales with the configured B
        bf = backscatter_fraction(phase_function(nap))
        @test backscattering(nap, 443.0) ≈ scattering(nap, 443.0) * bf
    end

    @testset "FellMineral beam-attenuation split" begin
        nap = NonAlgalParticles{Float64, FellMineral}(
            c_ref = 0.40, ω_nap = 0.95, λ_ref = 440.0)

        # c(440) = 0.40, split: a = (1-ω)·c = 0.02; b = ω·c = 0.38
        @test absorption(nap, 440.0) ≈ 0.02 rtol=1e-12
        @test scattering(nap, 440.0) ≈ 0.38 rtol=1e-12
        @test attenuation(nap, 440.0) ≈ 0.40 rtol=1e-12

        # Spectral form: c ∝ 1/λ
        @test attenuation(nap, 880.0) ≈ 0.20 rtol=1e-12
    end

    @testset "Construction validation" begin
        @test_throws ArgumentError NonAlgalParticles{Float64, FellMineral}(ω_nap = 1.1)
        @test_throws ArgumentError NonAlgalParticles{Float64, FellMineral}(c_ref = -1.0)
        @test_throws ArgumentError NonAlgalParticles{Float64, BabinBricaud2003}(a_ref = -0.1)
    end
end
