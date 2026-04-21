@testset "PureWater" begin

    @testset "Pope-Fry (1997) absorption" begin
        # Reference: Pope-Fry 1997 table values (bundled from omlc.org).
        # DESIGN §11 quotes a_w(440 nm) ≈ 0.00635 m⁻¹.
        w = PureWater{Float64}(absorption_model = PopeFry1997())
        @test absorption(w, 440.0) ≈ 0.00635 rtol=1e-3
        # Sanity: absorption is monotone in the mid-visible window
        @test absorption(w, 450.0) > absorption(w, 440.0)
        # Endpoints (Flat() extrapolation): Pope-Fry spans 380–727.5 nm.
        @test absorption(w, 380.0) ≈ 0.01137 rtol=1e-3
    end

    @testset "Smith-Baker (1981) absorption" begin
        w = PureWater{Float64}(absorption_model = SmithBaker1981())
        # Source value at 500 nm is 0.0204 1/cm * 100 = not in the SB81 table;
        # SB81 is 10 nm grid. At 500 nm exactly: 0.00257 1/m. Check a table pt.
        @test absorption(w, 500.0) > 0.0
        @test absorption(w, 500.0) < absorption(w, 200.0)   # UV is very absorbing
    end

    @testset "Combined piecewise (fallback to Smith-Baker <380 nm)" begin
        w = PureWater{Float64}(absorption_model = CombinedPopeFryMason())
        # 440 nm is in the Pope-Fry segment
        @test absorption(w, 440.0) ≈ 0.00635 rtol=1e-3
        # 350 nm falls back to Smith-Baker (since Mason not bundled)
        @test absorption(w, 350.0) > 0.0
        # 760 nm is in the Smith-Baker segment
        @test absorption(w, 760.0) > 0.5
    end

    @testset "MasonConeFry2016 errors helpfully without data" begin
        # The Mason-Cone-Fry CSV is not bundled (DESIGN §8.1).
        w = PureWater{Float64}(absorption_model = MasonConeFry2016())
        @test_throws ErrorException absorption(w, 300.0)
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

    @testset "Pegau temperature correction is zero outside 600–800 nm" begin
        w₂₀ = PureWater{Float64}(temperature = 20.0, absorption_model = PopeFry1997())
        w₃₀ = PureWater{Float64}(temperature = 30.0, absorption_model = PopeFry1997())
        # Outside calibrated window — temperature has no effect
        @test absorption(w₂₀, 450.0) == absorption(w₃₀, 450.0)
        # Inside calibrated window — 10 K gap ⇒ measurable shift
        @test absorption(w₃₀, 700.0) > absorption(w₂₀, 700.0)
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
