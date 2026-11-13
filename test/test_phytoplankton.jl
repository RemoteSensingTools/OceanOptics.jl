@testset "Phytoplankton" begin

    @testset "Gordon & Morel (1983) particulate scattering" begin
        # b_p(λ) = 0.30 · Chl^0.62 · (550/λ), Gordon & Morel (1983) as cited by
        # Fell (1997), Eq. 5.5. Morel & Maritorena (2001), Eq. 14, would give
        # 0.416 at Chl = 1 and 550 nm.
        p = Phytoplankton{Float64, Bricaud1998}(Chl = 1.0)
        @test scattering(p, 550.0) ≈ 0.30        rtol=1e-12
        @test scattering(p, 440.0) ≈ 0.30 * 550/440  rtol=1e-12

        # Chlorophyll scaling
        p2 = Phytoplankton{Float64, Bricaud1998}(Chl = 2.0)
        @test scattering(p2, 550.0) ≈ 0.30 * 2.0^0.62  rtol=1e-12

        # Backscattering derives from the Fournier-Forand phase function
        bf = backscatter_fraction(scattering_phase(p))
        @test backscattering(p, 550.0) ≈ scattering(p, 550.0) * bf
    end

    @testset "Bricaud et al. (1998) particulate absorption a_p = A·Chl^E" begin
        # Reference: the coefficients of HydroLight's AEmidUVabs.txt, whose header
        # reads "ap(wave) = A(wave)*[chl]^E(wave)": A(440) = 0.0520, E(440) = 0.6350
        # (Bricaud et al. 1998: a_p(440) = 0.052 Chl^0.635) and A(676) = 0.0203,
        # E(676) = 0.8223. The specific-absorption form A·Chl^(1−E) differs from
        # these by a factor 2.6 (440 nm) and 9.6 (676 nm) at Chl = 0.03.
        for Chl in (0.03, 0.1, 1.0, 10.0)
            p = Phytoplankton{Float64, Bricaud1998}(; Chl)
            @test absorption(p, 440.0) ≈ 0.0520 * Chl^0.6350 rtol = 1e-12
            @test absorption(p, 676.0) ≈ 0.0203 * Chl^0.8223 rtol = 1e-12
        end
        # The same numbers written out (four significant figures).
        @test absorption(Phytoplankton(Chl = 0.1), 440.0) ≈ 0.01205 rtol = 5e-4
        @test absorption(Phytoplankton(Chl = 10.0), 676.0) ≈ 0.1348 rtol = 5e-4
        # The blue-to-red ratio falls with chlorophyll (0.0520/0.0203·Chl^(0.6350−0.8223)).
        ratio(Chl) = (p = Phytoplankton(; Chl); absorption(p, 440.0) / absorption(p, 676.0))
        @test ratio(0.03) ≈ 4.94 rtol = 2e-3
        @test ratio(10.0) ≈ 1.66 rtol = 3e-3

        # No chlorophyll, no particulate absorption: zero at every tabulated
        # wavelength, also where E > 1 (698–700 nm).
        p0 = Phytoplankton(Chl = 0.0)
        @test all(iszero, absorption(p0, collect(300.0:2.0:1000.0)))
        @test absorption(Phytoplankton{Float32}(Chl = 0), 700f0) === 0f0
    end

    @testset "Gordon & Morel (1983) scattering: published values and parameters" begin
        # Fell (1997) Eq. 5.5: b_P(550) = 0.30 Chl^0.62. Values worked by hand:
        # 10^0.62 = 4.16869, 10^-0.62 = 0.239883.
        @test scattering(Phytoplankton(Chl = 10.0), 550.0) ≈ 1.25061 rtol = 1e-5
        @test scattering(Phytoplankton(Chl = 0.1), 550.0) ≈ 0.0719650 rtol = 1e-5
        @test scattering(Phytoplankton(Chl = 0.1), 440.0) ≈ 0.0719650 * 1.25 rtol = 1e-5
        gm = Phytoplankton(Chl = 1.0).scattering_model
        @test gm === GordonMorel1983() && (gm.b_ref, gm.chl_exponent, gm.λ_ref, gm.spectral_exponent) == (0.30, 0.62, 550.0, 1.0)
        # Morel & Maritorena (2001), Eq. 14 magnitude at 550 nm, as a non-default model.
        mm = Phytoplankton(Chl = 10.0, scattering_model = GordonMorel1983(b_ref = 0.416, chl_exponent = 0.766))
        @test scattering(mm, 550.0) ≈ 0.416 * 10^0.766 rtol = 1e-12
        @test_throws ArgumentError GordonMorel1983(b_ref = -0.1)
        @test_throws ArgumentError GordonMorel1983(λ_ref = 0.0)
        @test !isdefined(OceanOptics, :Gordon1992)
    end

    @testset "Construction sanity" begin
        @test_throws ArgumentError Phytoplankton{Float64, Bricaud1998}(Chl = -0.1)

        # Custom backscatter fraction routes into FF renormalization
        p = Phytoplankton{Float64, Bricaud1998}(Chl = 0.5, B = 0.015)
        @test backscatter_fraction(scattering_phase(p)) ≈ 0.015 rtol=0.05
    end
end
