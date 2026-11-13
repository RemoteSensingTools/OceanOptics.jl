using QuadGK

# Every exposed physical parameter is set to a non-default value and the result is
# compared with the closed form in that parameter: a parameter that is stored but not
# used fails here.
@testset "Exposed parameters reach the results" begin

    @testset "PureWater" begin
        w = PureWater{Float64}()
        @test (w.temperature, w.salinity, w.depolarization) == (20.0, 35.0, 0.039)
        # temperature: a(T) = a(22 °C) + Ψ_T (T − 22); Ψ_T(740 nm) ≈ 0.016 1/(m K) (Röttgers 2014)
        a22 = absorption(PureWater{Float64}(temperature = 22.0), 740.0)
        a32 = absorption(PureWater{Float64}(temperature = 32.0), 740.0)
        @test (a32 - a22) / 10 ≈ 0.016 atol = 1e-3
        @test absorption(PureWater{Float64}(temperature = 12.0), 740.0) - a22 ≈ -(a32 - a22) rtol = 1e-10
        # salinity: seawater scatters 30 % more than fresh water (Zhang, Hu & He 2009)
        @test scattering(w, 500.0) / scattering(PureWater{Float64}(salinity = 0.0), 500.0) ≈ 1.30 atol = 0.005
        @test scattering(PureWater{Float64}(temperature = 0.0), 500.0) != scattering(w, 500.0)
        # depolarization: magnitude and phase matrix
        w0 = PureWater{Float64}(depolarization = 0.0)
        @test scattering(w, 500.0) / scattering(w0, 500.0) ≈
              (6 + 6 * 0.039) / (6 - 7 * 0.039) * (2.039 / 1.039) / 2 rtol = 1e-12
        @test phase_function_moments(scattering_phase(w0), 2)[3] == 0.5
        # model choice by type
        @test scattering(PureWater{Float64}(scattering_model = Morel1974()), 550.0) ≈ 0.00193 rtol = 1e-12
        @test absorption(PureWater{Float64}(absorption_model = SmithBaker1981()), 500.0) ≈ 0.0257 rtol = 1e-6
    end

    @testset "RayleighWaterPhase" begin
        for ρ in (0.0, 0.039, 0.09, 0.2906)
            @test phase_function_moments(RayleighWaterPhase{Float64}(depolarization = ρ), 2)[3] ≈
                  (1 - ρ) / (2 + ρ) rtol = 1e-14
        end
    end

    @testset "Phytoplankton and GordonMorel1983" begin
        ref = Phytoplankton(Chl = 2.0)
        sc(; kw...) = Phytoplankton(Chl = 2.0, scattering_model = GordonMorel1983(; kw...))
        @test scattering(sc(b_ref = 0.45), 500.0) ≈ 1.5 * scattering(ref, 500.0) rtol = 1e-12
        @test scattering(sc(chl_exponent = 1.0), 500.0) ≈ scattering(ref, 500.0) * 2^0.38 rtol = 1e-12
        @test scattering(sc(λ_ref = 500.0), 500.0) ≈ 0.30 * 2^0.62 rtol = 1e-12
        @test scattering(sc(spectral_exponent = 2.0), 440.0) ≈ 0.30 * 2^0.62 * (550 / 440)^2 rtol = 1e-12
        @test scattering(sc(spectral_exponent = 0.0), 440.0) ≈ 0.30 * 2^0.62 rtol = 1e-12
        @test backscattering(Phytoplankton(Chl = 2.0, B = 0.02), 500.0) ≈ 0.02 * scattering(ref, 500.0) rtol = 1e-6
        @test backscattering(ref, 500.0) ≈ 0.01 * scattering(ref, 500.0) rtol = 1e-6
        # Float32 model, double-precision parameters: same law
        @test scattering(Phytoplankton{Float32}(Chl = 2, scattering_model = GordonMorel1983(b_ref = 0.45)), 500f0) ≈
              Float32(0.45 * 2^0.62 * 1.1) rtol = 1f-6
    end

    @testset "CDOM" begin
        c = CDOM{Float64}(a_ref = 0.2, λ_ref = 400.0, slope = 0.018)
        @test absorption(c, 400.0) == 0.2
        @test absorption(c, 450.0) ≈ 0.2 * exp(-0.018 * 50) rtol = 1e-13
        d = CDOM{Float64}()
        @test (d.a_ref, d.λ_ref, d.slope) == (0.03, 440.0, 0.014)       # slope: Fell (1997) Eq. 5.4
        @test absorption(d, 500.0) ≈ 0.03 * exp(-0.014 * 60) rtol = 1e-13
    end

    @testset "NonAlgalParticles" begin
        n = NonAlgalParticles{Float64}(a_ref = 0.1, b_ref = 0.6, λ_ref = 500.0, slope = 0.009, B = 0.03,
                                       spectral_exponent = 0.5)
        @test absorption(n, 550.0) ≈ 0.1 * exp(-0.009 * 50) rtol = 1e-13
        @test scattering(n, 600.0) ≈ 0.6 * sqrt(500 / 600) rtol = 1e-13
        @test backscattering(n, 600.0) ≈ 0.03 * scattering(n, 600.0) rtol = 1e-6
        d = NonAlgalParticles{Float64}()
        @test scattering(d, 550.0) ≈ 0.45 * 443 / 550 rtol = 1e-13
        @test absorption(d, 500.0) ≈ 0.05 * exp(-0.011 * 57) rtol = 1e-13
        f = NonAlgalParticles{Float64, FellMineral}(c_ref = 1.0, ω_nap = 0.8, λ_ref = 500.0, B = 0.05,
                                                    spectral_exponent = 2.0)
        @test absorption(f, 600.0) ≈ 0.2 * (500 / 600)^2 rtol = 1e-13
        @test scattering(f, 600.0) ≈ 0.8 * (500 / 600)^2 rtol = 1e-13
        @test backscattering(f, 600.0) ≈ 0.05 * scattering(f, 600.0) rtol = 1e-6
        # Defaults of Fell (1997) Eqs. 5.10–5.11: c ∝ 440/λ, ω₀ = 0.95, b_b/b = 0.036
        g = NonAlgalParticles{Float64, FellMineral}()
        @test absorption(g, 550.0) ≈ 0.05 * 0.40 * 440 / 550 rtol = 1e-12
        @test backscattering(g, 440.0) ≈ 0.95 * 0.036 * 0.40 rtol = 1e-6
    end

    @testset "Phase functions" begin
        @test phase_function_moments(HenyeyGreensteinPhase{Float64}(g = 0.6), 2) ≈ [1, 3 * 0.6, 5 * 0.36]
        @test HenyeyGreensteinPhase{Float64}().g == 0.85
        a = FournierForandPhase{Float64}(1.10, 3.5835)
        b = FournierForandPhase{Float64}(1.15, 3.5835)
        c = FournierForandPhase{Float64}(1.10, 4.0)
        @test backscatter_fraction(a) ≈ 0.0183 atol = 1e-4             # Mobley et al. (2002)
        @test backscatter_fraction(b) > backscatter_fraction(a)
        @test backscatter_fraction(c) > backscatter_fraction(a)
        # Mobley, Sundman & Boss (2002) family: n = 1.01 + 0.1542 (μ − 3)
        ff = FournierForandPhase{Float64}(backscatter_fraction = 0.05)
        @test ff.n ≈ 1.01 + 0.1542 * (ff.μ_j - 3) rtol = 1e-12
    end

    @testset "Refractive index" begin
        @test refractive_index(FixedRefractiveIndex(1.36), 500.0) == 1.36
        q(T, S) = refractive_index(QuanFry1995(temperature = T, salinity = S), 589.3)
        @test q(20.0, 35.0) - q(20.0, 0.0) ≈ 35 * 1.83e-4 rtol = 0.01      # dn/dS (Quan & Fry 1995)
        @test q(30.0, 35.0) < q(10.0, 35.0)
    end

    @testset "WaterRaman" begin
        r = WaterRaman{Float64}(shift_cm_inv = 3000.0, σ_cm_inv = 40.0, cross_section_ref = 5e-4,
                                cross_section_λ_ref = 500.0, cross_section_n = 4.0,
                                linear_depolarization = 0.3, excitation_range = (300.0, 600.0),
                                conversion = PhotonConversion())
        @test raman_peak_wavelength(r, 440.0) ≈ 1e7 / (1e7 / 440 - 3000) rtol = 1e-14
        λp = raman_peak_wavelength(r, 440.0)
        @test emission(r, 440.0, λp) ≈ 1e7 / λp^2 / (40 * sqrt(2π)) rtol = 1e-12
        @test excitation_absorption(r, 500.0) == 5e-4
        @test excitation_absorption(r, 400.0) ≈ 5e-4 * 1.25^4 rtol = 1e-13
        @test phase_function_moments(raman_phase(r), 2)[3] ≈ 0.7 / (2 * 1.6) rtol = 1e-14
        @test excitation_range(r) == (300.0, 600.0)
        @test radiance_source_kernel(r, 440.0, λp) ≈ inelastic_coefficient(r, 440.0, λp) * 440 / λp
        e = WaterRaman{Float64}(shift_cm_inv = 3000.0, σ_cm_inv = 40.0, cross_section_ref = 5e-4,
                                cross_section_λ_ref = 500.0, cross_section_n = 4.0)
        @test radiance_source_kernel(e, 440.0, λp) ≈ inelastic_coefficient(r, 440.0, λp)
        # Shared band defaults of the redistribution type and the process
        d = WaterRaman{Float64}(); s = WavenumberShiftEmission{Float64}()
        @test (d.shift_cm_inv, d.σ_cm_inv) == (s.shift_cm_inv, s.σ_cm_inv) == (3400.0, 85.0)
        @test redistribution(WavenumberShiftEmission{Float64}(shift_cm_inv = 3000.0, σ_cm_inv = 40.0), 440.0, λp) ==
              emission(r, 440.0, λp)
    end

    @testset "Emission bands" begin
        g = GaussianEmission{Float64}(center_nm = 700.0, σ_nm = 5.0)
        @test redistribution(g, 440.0, 700.0) ≈ 1 / (5 * sqrt(2π)) rtol = 1e-14
        @test redistribution(g, 440.0, 705.0) ≈ exp(-0.5) / (5 * sqrt(2π)) rtol = 1e-14
        d = GaussianEmission{Float64}()
        @test (d.center_nm, d.σ_nm) == (685.0, 10.6)                     # Fell (1997) §2.2.9
    end

    @testset "Fluorescence and absorbed light" begin
        phyto = Phytoplankton(Chl = 1.0, pigments = ProspectPigments(carotenoid_ratio = 0.3))
        mk(; φ = 0.02, car = 0.5, range = (400.0, 700.0)) =
            chlorophyll_fluorescence(phyto; quantum_yield = φ, carotenoid_transfer = car, excitation_range = range)
        f = mk()
        @test excitation_range(chlorophyll_fluorescence(phyto; quantum_yield = 0.02, carotenoid_transfer = 0.5)) ==
              (400.0, 700.0)
        @test excitation_absorption(mk(φ = 0.05), 450.0) ≈ 2.5 * excitation_absorption(f, 450.0) rtol = 1e-12
        a_chl, a_car, _ = absorption_by_part(phyto, 450.0)
        @test excitation_absorption(mk(car = 1.0), 450.0) - excitation_absorption(mk(car = 0.0), 450.0) ≈
              0.02 * a_car rtol = 1e-10
        @test absorption_part(Phytoplankton(Chl = 1.0, pigments = ProspectPigments(carotenoid_ratio = 0.6)),
                              :carotenoids, 450.0) > a_car
        λ = collect(350.0:5.0:750.0)
        layer(p) = layer_optics(OceanLayer(0.0, 2.0, AbstractOceanConstituent{Float64}[PureWater{Float64}(), phyto];
                                           inelastic = AbstractOceanInelasticProcess{Float64}[p]), λ; ℓ_max = 2)
        E = ones(length(λ)); λ_em = collect(650.0:2.0:800.0)
        P_par = inelastic_emission(layer(f), 1, E, λ_em)
        P_uv  = inelastic_emission(layer(mk(range = (350.0, 700.0))), 1, E, λ_em)
        P_set = inelastic_emission(layer(f), 1, E, λ_em; range = (350.0, 700.0))
        @test all(P_uv .> P_par) && P_set ≈ P_uv
        o = layer(f)
        # absorbed_par: range limits (water row; unit irradiance)
        par(r) = absorbed_par(o, E; range = r)[1]
        @test par((400, 700)) == absorbed_par(o, E)[1]
        @test par((400, 600)) + par((600, 700)) ≈ par((400, 700)) rtol = 1e-12
        @test par((350, 700)) > par((400, 700))
        # round-off allowance for negative irradiance
        noisy = copy(E); noisy[end] = -1e-6
        @test_throws ArgumentError absorbed_power(o, noisy)
        @test absorbed_power(o, noisy; roundoff_rtol = 1e-5)[:, 1:end-1] == absorbed_power(o, E)[:, 1:end-1]
        @test all(absorbed_power(o, noisy; roundoff_rtol = 1e-5)[:, end] .== 0)
        @test absorbed_par(o, noisy; roundoff_rtol = 1e-5) ≈ absorbed_par(o, E) rtol = 1e-6
        @test absorbed_photons(o, noisy; roundoff_rtol = 1e-5)[1, 1] == absorbed_photons(o, E)[1, 1]
        @test inelastic_emission(o, 1, noisy, λ_em; roundoff_rtol = 1e-5) == P_par
        tiny = copy(E); tiny[end] = -1e-12
        @test_throws ArgumentError absorbed_power(o, tiny; roundoff_rtol = 0)
        @test_throws ArgumentError absorbed_power(o, E; roundoff_rtol = -1)
        # CDOM fluorescence: yield and window
        cdom = CDOM{Float64}()
        cf(φ, r) = CDOMFluorescence(cdom; emission = GaussianEmission{Float64}(center_nm = 430.0, σ_nm = 20.0),
                                    quantum_yield = φ, excitation_range = r)
        @test excitation_absorption(cf(0.02, (300.0, 400.0)), 350.0) ≈ 0.02 * absorption(cdom, 350.0)
        @test excitation_range(cf(0.02, (310.0, 390.0))) == (310.0, 390.0)
    end
end
