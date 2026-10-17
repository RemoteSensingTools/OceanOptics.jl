@testset "PureWater" begin

    @testset "Pope-Fry (1997) absorption" begin
        # Reference: Pope-Fry 1997 table values (bundled from omlc.org), at the
        # table's 22 °C reference temperature (no Pegau correction).
        w = PureWater{Float64}(absorption_model = PopeFry1997(), temperature = 22.0)
        @test absorption(w, 440.0) ≈ 0.00635 rtol=1e-3
        @test absorption(w, 450.0) > absorption(w, 440.0)
        @test absorption(w, 380.0) ≈ 0.01137 rtol=1e-3
        # Outside 380–727.5 nm the table has no data: an error, not a clamp.
        @test_throws DomainError absorption(w, 375.0)
        @test_throws DomainError absorption(w, 750.0)
    end

    @testset "Smith-Baker (1981) absorption" begin
        w = PureWater{Float64}(absorption_model = SmithBaker1981())
        @test absorption(w, 500.0) > 0.0
        @test absorption(w, 500.0) < absorption(w, 200.0)   # UV is very absorbing
        # No temperature correction: Smith & Baker give no reference temperature.
        @test absorption(w, 700.0) == absorption(PureWater{Float64}(absorption_model = SmithBaker1981(),
                                                                     temperature = 5.0), 700.0)
        @test_throws DomainError absorption(w, 850.0)
    end

    @testset "IOCCG (2018) composite, the default" begin
        w = PureWater{Float64}(temperature = 22.0)
        @test w.absorption_model isa IOCCG2018
        @test absorption(w, 350.0) ≈ 0.0071 rtol=1e-4            # Morel et al. (2007) segment
        @test absorption(w, 440.0) ≈ 0.00635 rtol=1e-3           # Pope & Fry segment
        # Continued with Kou et al. (1993) above 750 nm (OCO/SIF wavelengths).
        kou = PureWater{Float64}(absorption_model = Kou1993(), temperature = 22.0)
        @test absorption(w, 760.0) == absorption(kou, 760.0)
        @test 2.0 < absorption(w, 760.0) < 3.5
        @test_throws DomainError absorption(w, 340.0)
    end

    @testset "Temperature correction relative to the table temperature" begin
        cold = PureWater{Float64}(absorption_model = PopeFry1997(), temperature = 12.0)
        ref  = PureWater{Float64}(absorption_model = PopeFry1997(), temperature = 22.0)
        @test absorption(cold, 500.0) == absorption(ref, 500.0)   # outside 600–800 nm
        @test absorption(cold, 720.0) < absorption(ref, 720.0)    # water absorbs more when warm
    end

    @testset "Morel (1974) scattering" begin
        w = PureWater{Float64}(scattering_model = Morel1974())
        # Reference: b(λ) = 0.00193 · (550/λ)^4.32, salinity-independent
        @test scattering(w, 550.0) ≈ 0.00193  rtol=1e-6
        @test scattering(w, 500.0) ≈ 0.00193 * (550/500)^4.32  rtol=1e-6
        @test scattering(w, 700.0) < scattering(w, 500.0)
    end

    @testset "Zhang-Hu (2009) scattering — faithful port" begin
        # Values checked against the Zhang et al. (2009) reference MATLAB
        # `betasw_ZHH2009.m` (ion-functions mirror). The port replaces the
        # Phase-1 calibrated approximation, so the numeric targets changed.
        w = PureWater{Float64}(scattering_model = ZhangHu2009(),
                               temperature = 20.0, salinity = 35.0)
        @test scattering(w, 500.0) ≈ 0.002547 rtol=1e-3

        # Pure-water limit (S = 0): β(90°) reduces to the density-only term.
        w_pure = PureWater{Float64}(scattering_model = ZhangHu2009(),
                                    temperature = 20.0, salinity = 0.0)
        @test scattering(w_pure, 500.0) ≈ 0.001958 rtol=1e-3

        # Morel (1966) measurement at 546 nm, S = 38.4 ‰, T = 20 °C — the
        # paper's Fig. 3 validation point. Agreement is 1 % per the paper.
        w_morel = PureWater{Float64}(scattering_model = ZhangHu2009(),
                                     temperature = 20.0, salinity = 38.4)
        @test scattering(w_morel, 546.0) ≈ 0.001791 rtol=1e-3

        # Temperature sensitivity: scattering decreases with warming near 20 °C
        w_warm = PureWater{Float64}(scattering_model = ZhangHu2009(),
                                    temperature = 30.0, salinity = 35.0)
        @test scattering(w_warm, 500.0) < scattering(w, 500.0)

        # Spectral slope is steeper than Rayleigh λ⁻⁴ at short λ (n(λ)
        # dispersion); roughly λ⁻⁴·³ per the paper conclusion.
        @test scattering(w, 400.0) / scattering(w, 700.0) > (700/400)^4.0
    end

    @testset "Backscatter: molecular is exactly b/2" begin
        w = PureWater{Float64}()
        for λ in (400.0, 500.0, 600.0, 700.0)
            @test backscattering(w, λ) ≈ 0.5 * scattering(w, λ)
        end
    end

    @testset "Röttgers et al. (2014) temperature correction" begin
        w₂₂ = PureWater{Float64}(temperature = 22.0, absorption_model = PopeFry1997())
        w₃₂ = PureWater{Float64}(temperature = 32.0, absorption_model = PopeFry1997())
        # Below 550 nm the coefficient is below its uncertainty: no correction.
        @test absorption(w₂₂, 450.0) == absorption(w₃₂, 450.0)
        # The O–H overtone near 740 nm: Ψ_T ≈ 0.016 1/(m K) (Röttgers Fig. 5),
        # far above the old five-point table's 0.0013.
        k = PureWater{Float64}(temperature = 32.0, absorption_model = Kou1993())
        k₀ = PureWater{Float64}(temperature = 22.0, absorption_model = Kou1993())
        @test 0.012 < (absorption(k, 740.0) - absorption(k₀, 740.0)) / 10 < 0.02
        # Ψ_T changes sign: slightly negative at 700 nm.
        @test absorption(w₃₂, 700.0) < absorption(w₂₂, 700.0)
    end

    @testset "Float32 parametric support" begin
        # Reference tables ship as Float64; arithmetic through them promotes
        # the result to Float64 regardless of the input type. We just check
        # that a `PureWater{Float32}` constructs and evaluates without
        # throwing, not that the output scalar type is preserved.
        w = PureWater{Float32}()
        @test isfinite(absorption(w, 440.0f0))
        @test isfinite(scattering(w, 500.0f0))
    end

    @testset "Vectorized wavelength grid" begin
        w = PureWater{Float64}(absorption_model = PopeFry1997(),
                               scattering_model = Morel1974())
        λ = [400.0, 450.0, 500.0, 550.0, 600.0]
        a = absorption(w, λ)
        b = scattering(w, λ)
        @test length(a) == length(λ)
        @test length(b) == length(λ)
        @test a[2] ≈ absorption(w, 450.0)
        @test b[3] ≈ scattering(w, 500.0)
    end
end
