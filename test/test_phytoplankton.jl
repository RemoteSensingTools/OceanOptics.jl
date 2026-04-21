@testset "Phytoplankton" begin

    @testset "Morel & Maritorena (2001) scattering" begin
        # b_p(λ) = 0.30 · Chl^0.62 · (550/λ)
        p = Phytoplankton{Float64, Bricaud1995}(Chl = 1.0)
        @test scattering(p, 550.0) ≈ 0.30        rtol=1e-12
        @test scattering(p, 440.0) ≈ 0.30 * 550/440  rtol=1e-12

        # Chlorophyll scaling
        p2 = Phytoplankton{Float64, Bricaud1995}(Chl = 2.0)
        @test scattering(p2, 550.0) ≈ 0.30 * 2.0^0.62  rtol=1e-12

        # Backscattering derives from the Fournier-Forand phase function
        bf = backscatter_fraction(phase_function(p))
        @test backscattering(p, 550.0) ≈ scattering(p, 550.0) * bf
    end

    @testset "Bricaud (1995) absorption errors helpfully without data" begin
        p = Phytoplankton{Float64, Bricaud1995}(Chl = 1.0)
        @test_throws ErrorException absorption(p, 440.0)
    end

    @testset "Gordon (1992) absorption is a reserved stub" begin
        p = Phytoplankton{Float64, Gordon1992}(Chl = 1.0)
        @test_throws ErrorException absorption(p, 440.0)
    end

    @testset "Construction sanity" begin
        @test_throws ArgumentError Phytoplankton{Float64, Bricaud1995}(Chl = -0.1)

        # Custom backscatter fraction routes into FF renormalization
        p = Phytoplankton{Float64, Bricaud1995}(Chl = 0.5, B = 0.015)
        @test backscatter_fraction(phase_function(p)) ≈ 0.015 rtol=0.05
    end
end
