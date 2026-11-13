using QuadGK

@testset "Water Raman" begin

    @testset "Default Haltrin-Kattawar parameterization" begin
        ram = WaterRaman{Float64}()

        # Defaults: Δν = 3400 cm⁻¹, σ_ν = 85, cross_section 2.7e-4 at 488 nm, n=5.5
        @test ram.shift_cm_inv        == 3400.0
        @test ram.σ_cm_inv            == 85.0
        @test ram.cross_section_ref   == 2.7e-4
        @test ram.cross_section_λ_ref == 488.0
        @test ram.cross_section_n     == 5.5

        # Peak wavelength: λ' = 440 nm → ≈ 517.4 nm (Stokes-shifted O-H stretch)
        @test raman_peak_wavelength(ram, 440.0) ≈ 517.4 rtol=1e-3

        # λ' = 532 nm (common Nd:YAG-doubled lidar band) → ≈ 649.5 nm
        @test raman_peak_wavelength(ram, 532.0) ≈ 649.5 rtol=1e-3
    end

    @testset "excitation_absorption: λ^{-5.5} cross-section" begin
        ram = WaterRaman{Float64}()

        # At the reference wavelength the cross-section equals `cross_section_ref`.
        @test excitation_absorption(ram, 488.0) ≈ ram.cross_section_ref rtol=1e-12

        # Spectral slope: ratio at 400 vs 500 nm equals (500/400)^5.5
        ratio = excitation_absorption(ram, 400.0) / excitation_absorption(ram, 500.0)
        @test ratio ≈ (ram.cross_section_λ_ref/400.0)^5.5 /
                      (ram.cross_section_λ_ref/500.0)^5.5  rtol=1e-12

        # UV is stronger, NIR weaker
        @test excitation_absorption(ram, 400.0) > excitation_absorption(ram, 700.0)
    end

    @testset "Emission spectrum properties" begin
        ram = WaterRaman{Float64}()

        # Normalization in wavelength space (Jacobian compensation)
        λ_prime = 440.0
        I, _ = quadgk(λ -> emission(ram, λ_prime, λ), 450.0, 650.0; rtol = 1e-9)
        @test I ≈ 1.0 rtol=1e-4

        # Emission peak is at the Raman-shifted wavelength
        λ_peak = raman_peak_wavelength(ram, λ_prime)
        λ_grid = range(λ_peak - 10.0, λ_peak + 10.0; length = 2001)
        f_grid = [emission(ram, λ_prime, λ) for λ in λ_grid]
        (_, i) = findmax(f_grid)
        @test λ_grid[i] ≈ λ_peak atol = 0.05
    end

    # b_R is an energy conversion coefficient (Mobley et al. 1993, Problem 7 and
    # Eqs. A3–A4): the energy kernel integrates over emission wavelength to b_R(λ′).
    # The reference b_R is written out here, not taken from the package.
    @testset "Energy convention: ∫ K_E dλ = b_R(λ′)" begin
        ram = WaterRaman{Float64}()
        @test conversion_units(ram) === EnergyConversion()
        @test ram isa WaterRaman{Float64, EnergyConversion}
        for λ′ in (350.0, 417.0, 488.0, 600.0)
            b_R = 2.7e-4 * (488.0 / λ′)^5.5
            λp  = 1e7 / (1e7 / λ′ - 3400.0)
            IE, _ = quadgk(λ -> radiance_source_kernel(ram, λ′, λ), 0.8λp, λp, 1.25λp; rtol = 1e-10)
            @test IE ≈ b_R rtol = 1e-8
            # Photon kernel = energy kernel × λ/λ′: more photons leave than b_R counts,
            # by the mean of λ/λ′ over the band (≈ λ_peak/λ′).
            IN, _ = quadgk(λ -> inelastic_coefficient(ram, λ′, λ), 0.8λp, λp, 1.25λp; rtol = 1e-10)
            @test IN > b_R
            @test IN / b_R ≈ λp / λ′ rtol = 2e-3
            @test inelastic_coefficient(ram, λ′, λp) ≈ radiance_source_kernel(ram, λ′, λp) * λp / λ′
        end
    end

    @testset "Photon convention: ∫ K_N dλ = b_R(λ′)" begin
        ram = WaterRaman{Float64}(conversion = PhotonConversion())
        @test conversion_units(ram) === PhotonConversion()
        # The energy form of the photon convention integrates to λ′/λ_peak of b_R, the
        # ratio for Δν = 3400 cm⁻¹ (0.858 = 417/486 at 417 nm excitation).
        for (λ′, low) in ((350.0, 0.881), (417.0, 0.858), (488.0, 0.834), (600.0, 0.796))
            b_R = 2.7e-4 * (488.0 / λ′)^5.5
            λp  = 1e7 / (1e7 / λ′ - 3400.0)
            IN, _ = quadgk(λ -> inelastic_coefficient(ram, λ′, λ), 0.8λp, λp, 1.25λp; rtol = 1e-10)
            IE, _ = quadgk(λ -> radiance_source_kernel(ram, λ′, λ), 0.8λp, λp, 1.25λp; rtol = 1e-10)
            @test IN ≈ b_R rtol = 1e-8
            @test IE / b_R ≈ low atol = 1e-3
        end
    end

    # Mobley et al. (1993), Problem 7: b_R = b_w(417) = 0.0063 1/m, all light Raman
    # scattered at 417 nm appears at 486 nm; the band-integrated source of radiance is
    # b_R E₀ / (4π), with no photon-energy ratio.
    @testset "Mobley et al. (1993) Problem 7 source" begin
        ram = WaterRaman{Float64}(cross_section_ref = 0.0063, cross_section_λ_ref = 417.0,
                                  linear_depolarization = 0.17, excitation_range = (410.0, 424.0))
        @test raman_peak_wavelength(ram, 417.0) ≈ 486.0 atol = 0.2
        Δz    = 2.0
        layer = OceanLayer(0.0, Δz, AbstractOceanConstituent{Float64}[PureWater{Float64}()];
                           inelastic = AbstractOceanInelasticProcess{Float64}[ram])
        λ  = collect(410.0:1.0:424.0)
        o  = layer_optics(layer, λ; ℓ_max = 2)
        # A line of unit irradiance (1 W m⁻²) at 417 nm: triangle of base 2 nm, height 1.
        E₀ = [x == 417.0 ? 1.0 : 0.0 for x in λ]
        λ_em = collect(460.0:0.02:515.0)
        P    = inelastic_emission(o, 1, E₀, λ_em)                   # W m⁻² nm⁻¹, all directions
        source = sum(P) * 0.02 / (4π * Δz)                          # W m⁻³ sr⁻¹
        @test source ≈ 0.0063 * 1.0 / (4π) rtol = 1e-6
        # The photon convention is 417/486 = 14.2 % lower for the same coefficient.
        ramN = WaterRaman{Float64}(cross_section_ref = 0.0063, cross_section_λ_ref = 417.0,
                                   excitation_range = (410.0, 424.0), conversion = PhotonConversion())
        oN = layer_optics(OceanLayer(0.0, Δz, AbstractOceanConstituent{Float64}[PureWater{Float64}()];
                                     inelastic = AbstractOceanInelasticProcess{Float64}[ramN]), λ; ℓ_max = 2)
        PN = inelastic_emission(oN, 1, E₀, λ_em)
        @test sum(PN) / sum(P) ≈ 417 / 486 rtol = 1e-3
    end

    # Mobley et al. (1993) Eq. 9 with ρ = 0.17, typed here from the paper:
    # p_R = (3/16π) (1+3ρ)/(1+2ρ) [1 + (1−ρ)/(1+3ρ) cos²Θ].
    @testset "Raman phase function and matrix" begin
        ram = WaterRaman{Float64}()
        @test !is_isotropic(ram)
        @test ram.linear_depolarization == 0.17
        ρ   = 0.17
        ρ_n = natural_depolarization(ram)
        @test ρ_n ≈ 2ρ / (1 + ρ) rtol = 1e-14
        @test ρ_n ≈ 0.2906 atol = 5e-5
        pf = raman_phase(ram)
        @test pf isa RayleighWaterPhase{Float64} && pf.depolarization == ρ_n
        mobley(μ) = 3 / (16π) * (1 + 3ρ) / (1 + 2ρ) * (1 + (1 - ρ) / (1 + 3ρ) * μ^2)
        for μ in (-1.0, -0.4, 0.0, 0.3, 0.8, 1.0)
            @test phase_function_value(pf, μ) / (4π) ≈ mobley(μ) rtol = 1e-13
        end
        I, _ = quadgk(mobley, -1, 1)
        @test 2π * I ≈ 1 rtol = 1e-12
        # β₂ of RayleighWaterPhase for ρ_n equals (1−ρ)/(2(1+2ρ)) = 0.3097.
        β = phase_function_moments(pf, 4)
        @test β[3] ≈ (1 - ρ) / (2 * (1 + 2ρ)) rtol = 1e-14
        @test β[3] ≈ 0.3097 atol = 5e-5
        @test β[1] == 1 && β[2] == 0 && β[4] == 0 && β[5] == 0
        # Phase matrix: Rayleigh matrix with depolarization ρ_n (Hansen & Travis 1974).
        g = greek_expansion(pf, 2; polarized = true)
        Δ = (1 - ρ_n) / (1 + ρ_n / 2)
        @test g.β[3] ≈ Δ / 2 && g.α[3] ≈ 3Δ && g.γ[3] ≈ sqrt(1.5) * Δ
        @test g.δ[2] ≈ 1.5 * Δ * (1 - 2ρ_n) / (1 - ρ_n)
        # A different ρ changes the phase function (the stored value is the one used).
        @test phase_function_moments(raman_phase(WaterRaman{Float64}(linear_depolarization = 0.0)), 2)[3] == 0.5
        @test phase_function_moments(raman_phase(WaterRaman{Float32}()), 2)[3] isa Float32
    end

    @testset "excitation_range + constructor validation" begin
        ram = WaterRaman{Float64}(excitation_range = (400.0, 650.0))
        @test excitation_range(ram) == (400.0, 650.0)

        @test_throws ArgumentError WaterRaman{Float64}(σ_cm_inv = 0.0)
        @test_throws ArgumentError WaterRaman{Float64}(cross_section_ref = -1e-4)
        @test_throws ArgumentError WaterRaman{Float64}(linear_depolarization = 0.8)
        @test_throws ArgumentError WaterRaman{Float64}(excitation_range = (700.0, 400.0))
    end

    @testset "Forward-mode AD through WaterRaman" begin
        ram = WaterRaman{Float64}()
        # ∂a^R/∂λ' = -n · a^R / λ' (analytic from λ^{-n} power law)
        for λ′ in (400.0, 500.0, 600.0)
            d_ad  = ForwardDiff.derivative(x -> excitation_absorption(ram, x), λ′)
            d_ref = -ram.cross_section_n * excitation_absorption(ram, λ′) / λ′
            @test d_ad ≈ d_ref rtol=1e-10
        end
    end
end
