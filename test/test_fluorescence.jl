@testset "Fluorescence processes" begin

    @testset "IsotropicFluorescence with a CDOM absorber" begin
        # Use CDOM as the driving absorber so we get a concrete absorption
        # without relying on Bricaud data. Physics is still meaningful
        # (conceptual CDOM → visible fluorescence, Coble-style).
        cdom = CDOM{Float64}(a_ref = 0.1, λ_ref = 440.0, slope = 0.014)
        em   = GaussianEmission{Float64}(center_nm = 520.0, σ_nm = 20.0)
        sif  = IsotropicFluorescence{Float64}(cdom;
                  emission         = em,
                  quantum_yield    = 0.005,
                  excitation_range = (350.0, 500.0))

        # excitation_absorption = φ · a_absorber(λ')
        @test excitation_absorption(sif, 440.0) ≈ 0.005 * 0.1     rtol=1e-12
        @test excitation_absorption(sif, 500.0) ≈ 0.005 * absorption(cdom, 500.0)

        # emission delegates to the redistribution
        @test emission(sif, 440.0, 520.0) ≈ redistribution(em, 440.0, 520.0)

        # Integration window
        @test excitation_range(sif) == (350.0, 500.0)
        @test is_isotropic(sif)
        @test quantum_yield(sif)    == 0.005

        # inelastic_coefficient = a^I · f
        @test inelastic_coefficient(sif, 440.0, 520.0) ≈
              excitation_absorption(sif, 440.0) * emission(sif, 440.0, 520.0)

        # Constructor validation
        @test_throws ArgumentError IsotropicFluorescence{Float64}(cdom;
            emission      = em, quantum_yield = -0.1)
        @test_throws ArgumentError IsotropicFluorescence{Float64}(cdom;
            emission = em, quantum_yield = 1.5)
        @test_throws ArgumentError IsotropicFluorescence{Float64}(cdom;
            emission = em, excitation_range = (500.0, 400.0))
    end

    @testset "chlorophyll_fluorescence: Fell (1997) parameterization" begin
        phyto = Phytoplankton{Float64, Bricaud1998}(Chl = 0.5)
        sif   = chlorophyll_fluorescence(phyto)

        # Type-parameter plumbing
        @test sif isa IsotropicFluorescence{Float64}
        @test sif.quantum_yield == 0.003
        @test excitation_range(sif) == (400.0, 700.0)
        @test is_isotropic(sif)

        # Emission peaks at 685 nm with σ = 10.6 nm, independent of λ'
        @test sif.emission isa GaussianEmission{Float64}
        @test sif.emission.center_nm == 685.0
        @test sif.emission.σ_nm      == 10.6

        # Emission density at the peak and at the half-max
        f_peak = emission(sif, 500.0, 685.0)
        f_half = emission(sif, 500.0, 685.0 + sqrt(2 * log(2)) * 10.6)
        @test f_peak > 0
        @test f_half ≈ f_peak / 2 rtol=1e-10

        # excitation_absorption now works with the Bricaud 1998 data bundled.
        # a^I(λ') = φ · a_absorber(λ'), so positive and much smaller than a_φ.
        a_phyto = absorption(phyto, 440.0)
        @test excitation_absorption(sif, 440.0) ≈ 0.003 * a_phyto rtol=1e-12

        # Custom parameterization
        sif_ens = chlorophyll_fluorescence(phyto;
                      peak_nm = 683.0, σ_nm = 12.0, quantum_yield = 0.005)
        @test sif_ens.emission.center_nm == 683.0
        @test sif_ens.emission.σ_nm      == 12.0
        @test sif_ens.quantum_yield      == 0.005
    end

    @testset "CDOMFluorescence" begin
        cdom = CDOM{Float64}(a_ref = 0.5, λ_ref = 440.0, slope = 0.017)
        em   = GaussianEmission{Float64}(center_nm = 440.0, σ_nm = 40.0)
        fluo = CDOMFluorescence{Float64}(cdom;
                  emission         = em,
                  quantum_yield    = 0.01,
                  excitation_range = (300.0, 400.0))

        @test excitation_absorption(fluo, 350.0) ≈ 0.01 * absorption(cdom, 350.0)
        @test emission(fluo, 350.0, 440.0)       == redistribution(em, 350.0, 440.0)
        @test excitation_range(fluo)             == (300.0, 400.0)
        @test is_isotropic(fluo)
        @test quantum_yield(fluo)                == 0.01

        # Validation
        @test_throws ArgumentError CDOMFluorescence{Float64}(cdom;
            emission = em, quantum_yield = -0.01)
    end
end
