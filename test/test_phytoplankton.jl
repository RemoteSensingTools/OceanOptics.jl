@testset "Phytoplankton" begin

    @testset "Morel & Maritorena (2001) scattering" begin
        # b_p(λ) = 0.30 · Chl^0.62 · (550/λ)
        p = Phytoplankton{Float64, Bricaud1998}(Chl = 1.0)
        @test scattering(p, 550.0) ≈ 0.30        rtol=1e-12
        @test scattering(p, 440.0) ≈ 0.30 * 550/440  rtol=1e-12

        # Chlorophyll scaling
        p2 = Phytoplankton{Float64, Bricaud1998}(Chl = 2.0)
        @test scattering(p2, 550.0) ≈ 0.30 * 2.0^0.62  rtol=1e-12

        # Backscattering derives from the Fournier-Forand phase function
        bf = backscatter_fraction(phase_function(p))
        @test backscattering(p, 550.0) ≈ scattering(p, 550.0) * bf
    end

    @testset "Bricaud (1998 + UV-extended) absorption" begin
        # a_φ(λ, Chl=1) = A(λ)·1^{1-E(λ)} = A(λ) exactly; tabulated A(440) is
        # around 0.0394 m²/mg (Bricaud 1998 + Morrison-Nelson + Vasilkov update).
        p = Phytoplankton{Float64, Bricaud1998}(Chl = 1.0)
        @test 0.02 < absorption(p, 440.0) < 0.08

        # Red chlorophyll peak (676 nm) stronger than the green minimum (600 nm)
        @test absorption(p, 676.0) > absorption(p, 600.0)

        # a_φ ∝ Chl^{1-E} with E < 1 ⇒ quadrupling Chl gives less than 4×
        p2 = Phytoplankton{Float64, Bricaud1998}(Chl = 4.0)
        @test absorption(p2, 440.0) > absorption(p, 440.0)
        @test absorption(p2, 440.0) < 4 * absorption(p, 440.0)
    end

    @testset "Gordon (1992) absorption is a reserved stub" begin
        p = Phytoplankton{Float64, Gordon1992}(Chl = 1.0)
        @test_throws ErrorException absorption(p, 440.0)
    end

    @testset "Construction sanity" begin
        @test_throws ArgumentError Phytoplankton{Float64, Bricaud1998}(Chl = -0.1)

        # Custom backscatter fraction routes into FF renormalization
        p = Phytoplankton{Float64, Bricaud1998}(Chl = 0.5, B = 0.015)
        @test backscatter_fraction(phase_function(p)) ≈ 0.015 rtol=0.05
    end
end
