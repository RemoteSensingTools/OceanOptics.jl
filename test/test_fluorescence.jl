using QuadGK

@testset "Fluorescence processes" begin

    @testset "IsotropicFluorescence with a CDOM absorber" begin
        # CDOM as a concrete driving absorber; the physics (UV-excited
        # visible fluorescence) is still meaningful.
        cdom = CDOM{Float64}(a_ref = 0.1, λ_ref = 440.0, slope = 0.014)
        em   = GaussianEmission{Float64}(center_nm = 520.0, σ_nm = 20.0)
        sif  = IsotropicFluorescence(cdom; emission = em, quantum_yield = 0.005,
                                     excitation_range = (350.0, 500.0))
        # excitation_absorption = φ · a_absorber(λ')
        @test excitation_absorption(sif, 440.0) ≈ 0.005 * 0.1 rtol = 1e-12
        @test excitation_absorption(sif, 500.0) ≈ 0.005 * absorption(cdom, 500.0)
        @test emission(sif, 440.0, 520.0) ≈ redistribution(em, 440.0, 520.0)
        @test excitation_range(sif) == (350.0, 500.0)
        @test is_isotropic(sif) && quantum_yield(sif) == 0.005
        @test only(excitation_pathways(sif)) == ExcitationPathway(cdom)
        @test inelastic_coefficient(sif, 440.0, 520.0) ≈
              excitation_absorption(sif, 440.0) * emission(sif, 440.0, 520.0)
        # Constructor validation; the yield and range are always explicit.
        kw = (emission = em, excitation_range = (350.0, 500.0))
        @test_throws ArgumentError IsotropicFluorescence(cdom; kw..., quantum_yield = -0.1)
        @test_throws ArgumentError IsotropicFluorescence(cdom; kw..., quantum_yield = 1.5)
        @test_throws ArgumentError IsotropicFluorescence(cdom; emission = em, quantum_yield = 0.01,
                                                         excitation_range = (500.0, 400.0))
        @test_throws UndefKeywordError IsotropicFluorescence(cdom; kw...)
        @test_throws ArgumentError IsotropicFluorescence(ExcitationPathway[]; kw..., quantum_yield = 0.01)
    end

    @testset "Excitation pathways" begin
        phyto = Phytoplankton(Chl = 0.5, pigments = ProspectPigments(carotenoid_ratio = 0.3))
        @test ExcitationPathway(phyto, :carotenoids, 0.4).efficiency == 0.4
        @test_throws ArgumentError ExcitationPathway(phyto, :total)            # not a part of it
        @test_throws ArgumentError ExcitationPathway(phyto, :chlorophyll, 1.2)
        f = IsotropicFluorescence([ExcitationPathway(phyto, :chlorophyll, 1.0),
                                   ExcitationPathway(phyto, :carotenoids, 0.4)];
                                  emission = chloroplast_emission(), quantum_yield = 0.02,
                                  excitation_range = (400.0, 700.0))
        for λ in (440.0, 490.0, 650.0)
            achl, acar, _ = absorption_by_part(phyto, λ)
            @test excitation_absorption(f, λ) ≈ 0.02 * (achl + 0.4acar) rtol = 1e-12
        end
    end

    @testset "chloroplast_emission (Fluspect)" begin
        em = chloroplast_emission()
        @test em isa TabulatedEmission{Float64} && em.λ[[1, end]] == [640.0, 851.0]
        x = 630.0:0.25:860.0
        area(r) = sum(redistribution(r, 0.0, xi) for xi in x) * step(x)
        @test area(em) ≈ 1 rtol = 1e-10
        # Combined spectrum: maxima at 685 and 739 nm, as in Fluspect's phi.
        f = [redistribution(em, 0.0, λ) for λ in 640.0:851.0]
        @test 640 + argmax(f) - 1 == 685
        @test 720 + argmax(f[81:end]) - 1 == 739        # second maximum, beyond the 685 nm band
        # Photosystem mixture: PSI peaks at 728 nm, PSII at 682 nm.
        psi, psii = chloroplast_emission(psi_share = 1), chloroplast_emission(psi_share = 0)
        @test 640 + argmax([redistribution(psi, 0.0, λ) for λ in 640.0:851.0]) - 1 == 728
        @test 640 + argmax([redistribution(psii, 0.0, λ) for λ in 640.0:851.0]) - 1 == 682
        mix = chloroplast_emission(psi_share = 0.3)
        @test redistribution(mix, 0.0, 700.0) ≈
              0.3redistribution(psi, 0.0, 700.0) + 0.7redistribution(psii, 0.0, 700.0) rtol = 1e-6
        @test area(mix) ≈ 1 rtol = 1e-10
        @test chloroplast_emission(Float32) isa TabulatedEmission{Float32}
        @test_throws ArgumentError chloroplast_emission(psi_share = 1.5)
    end

    @testset "chlorophyll_fluorescence" begin
        # Without a pigment partition: all phytoplankton absorption excites.
        phyto = Phytoplankton(Chl = 0.5)
        sif   = chlorophyll_fluorescence(phyto; quantum_yield = 0.02)
        @test sif isa IsotropicFluorescence{Float64}
        @test sif.emission isa TabulatedEmission && excitation_range(sif) == (400.0, 700.0)
        @test excitation_absorption(sif, 440.0) ≈ 0.02 * absorption(phyto, 440.0) rtol = 1e-12
        @test_throws UndefKeywordError chlorophyll_fluorescence(phyto)
        @test_throws ArgumentError chlorophyll_fluorescence(phyto; quantum_yield = 0.02,
                                                            carotenoid_transfer = 0.5)
        # The single-band Gaussian of Fell (1997) remains available.
        fell = chlorophyll_fluorescence(phyto; quantum_yield = 0.003,
                                        emission = GaussianEmission(center_nm = 685.0, σ_nm = 10.6))
        @test emission(fell, 500.0, 685.0) ≈ 1 / (10.6 * sqrt(2π))

        # Land analogy: chlorophyll only vs. carotenoids as (partial) antenna.
        p = Phytoplankton(Chl = 0.5, pigments = ProspectPigments(carotenoid_ratio = 0.3))
        @test_throws ArgumentError chlorophyll_fluorescence(p; quantum_yield = 0.02)   # choose η
        only_chl = chlorophyll_fluorescence(p; quantum_yield = 0.02, carotenoid_transfer = 0.0)
        antenna  = chlorophyll_fluorescence(p; quantum_yield = 0.02, carotenoid_transfer = 1.0)
        half     = chlorophyll_fluorescence(p; quantum_yield = 0.02, carotenoid_transfer = 0.5)
        @test [x.efficiency for x in excitation_pathways(half)] == [1.0, 0.5, 0.0]
        for λ in (450.0, 490.0, 520.0)
            achl, acar, aun = absorption_by_part(p, λ)
            @test excitation_absorption(only_chl, λ) ≈ 0.02achl rtol = 1e-12
            @test excitation_absorption(antenna, λ) ≈ 0.02(achl + acar) rtol = 1e-12
            @test excitation_absorption(half, λ) ≈ 0.02(achl + 0.5acar) rtol = 1e-12
        end
        @test excitation_absorption(only_chl, 490.0) < excitation_absorption(antenna, 490.0)
        # Beyond the carotenoid band the choice does not matter.
        @test excitation_absorption(only_chl, 650.0) ≈ excitation_absorption(antenna, 650.0)

        # Ocean pigment groups: separate efficiencies for photosynthetic and
        # photoprotective carotenoids; phycobilins default to 1.
        pig = Bidigare1990Pigments(chl_a = 0.5, chl_c = 0.1, psc = 0.3, ppc = 0.15, pe_peb = 0.05)
        pb  = Phytoplankton(Chl = 0.5, pigments = pig)
        f   = chlorophyll_fluorescence(pb; quantum_yield = 0.02,
                                       transfer = (photosynthetic_carotenoids = 0.8,
                                                   photoprotective_carotenoids = 0.0))
        η   = Dict(x.part => x.efficiency for x in excitation_pathways(f))
        @test η[:chlorophyll_a] == η[:chlorophyll_c] == 1
        @test η[:photosynthetic_carotenoids] == 0.8 && η[:photoprotective_carotenoids] == 0
        @test η[:phycoerythrin_peb] == 1 && η[:unassigned] == 0
        a = absorption_by_part(pb, 470.0)
        @test excitation_absorption(f, 470.0) ≈ 0.02(a[1] + a[2] + a[3] + 0.8a[4] + a[6] + a[7]) rtol = 1e-12

        # Defaults: carotenoid_transfer reaches the photosynthetic carotenoids only;
        # photoprotective carotenoids dissipate (η = 0) unless photoprotective_transfer says otherwise.
        fd = chlorophyll_fluorescence(pb; quantum_yield = 0.02, carotenoid_transfer = 0.7)
        ηd = Dict(x.part => x.efficiency for x in excitation_pathways(fd))
        @test ηd[:photosynthetic_carotenoids] == 0.7 && ηd[:photoprotective_carotenoids] == 0
        for λ in (360.0, 440.0, 470.0)
            a = absorption_by_part(pb, λ)
            @test a[5] > 0
            @test excitation_absorption(fd, λ) ≈ 0.02(a[1] + a[2] + a[3] + 0.7a[4] + a[6] + a[7]) rtol = 1e-12
        end
        fp = chlorophyll_fluorescence(pb; quantum_yield = 0.02, carotenoid_transfer = 0.7,
                                      photoprotective_transfer = 0.25)
        a = absorption_by_part(pb, 440.0)
        @test excitation_absorption(fp, 440.0) - excitation_absorption(fd, 440.0) ≈ 0.02 * 0.25a[5] rtol = 1e-12
        # carotenoid_transfer is required while a photosynthetic carotenoid has no η ...
        @test_throws ArgumentError chlorophyll_fluorescence(pb; quantum_yield = 0.02)
        # ... and not when `transfer` gives it (photoprotective carotenoids default to 0).
        ft = chlorophyll_fluorescence(pb; quantum_yield = 0.02,
                                      transfer = (photosynthetic_carotenoids = 0.8,))
        @test Dict(x.part => x.efficiency for x in excitation_pathways(ft))[:photoprotective_carotenoids] == 0
        @test_throws ArgumentError chlorophyll_fluorescence(pb; quantum_yield = 0.02,
                                                            carotenoid_transfer = 0.5,
                                                            transfer = (fucoxanthin = 1.0,))
        g = chlorophyll_fluorescence(pb; quantum_yield = 0.02, carotenoid_transfer = 0.5,
                                     phycobilin_transfer = 0.9)
        @test Dict(x.part => x.efficiency for x in excitation_pathways(g))[:phycoerythrin_pub] == 0.9
    end

    @testset "CDOMFluorescence" begin
        cdom = CDOM{Float64}(a_ref = 0.5, λ_ref = 440.0, slope = 0.017)
        em   = GaussianEmission{Float64}(center_nm = 440.0, σ_nm = 40.0)
        fluo = CDOMFluorescence(cdom; emission = em, quantum_yield = 0.01)
        @test excitation_absorption(fluo, 350.0) ≈ 0.01 * absorption(cdom, 350.0)
        @test emission(fluo, 350.0, 440.0) == redistribution(em, 350.0, 440.0)
        @test excitation_range(fluo) == (300.0, 400.0)
        @test is_isotropic(fluo) && quantum_yield(fluo) == 0.01
        @test only(excitation_pathways(fluo)).absorber === cdom
        @test_throws ArgumentError CDOMFluorescence(cdom; emission = em, quantum_yield = -0.01)
        @test_throws UndefKeywordError CDOMFluorescence(cdom; emission = em)
    end

    # Excitation below the 400 nm start of the pigment tables: the shares at 400 nm
    # apply, so a partition whose parts all transfer with η = 1 excites exactly like
    # the unpartitioned phytoplankton (reference: the total absorption itself).
    @testset "UV excitation with a pigment partition" begin
        plain = Phytoplankton(Chl = 0.5)
        whole = chlorophyll_fluorescence(plain; quantum_yield = 0.02, excitation_range = (350.0, 700.0))
        for pig in (ProspectPigments(carotenoid_ratio = 0.3),
                    Bidigare1990Pigments(chl_a = 0.5, chl_c = 0.1, psc = 0.3, ppc = 0.15))
            phyto = Phytoplankton(Chl = 0.5, pigments = pig)
            f = chlorophyll_fluorescence(phyto; quantum_yield = 0.02, carotenoid_transfer = 1.0,
                                         photoprotective_transfer = 1.0,
                                         excitation_range = (350.0, 700.0))
            # Continuous across 400 nm (a_p changes by < 0.1 % per 0.001 nm).
            @test excitation_absorption(f, 399.999) ≈ excitation_absorption(f, 400.0) rtol = 1e-4
            @test excitation_absorption(f, 399.0) ≈ 0.02 * absorption(plain, 399.0) rtol = 1e-12
            @test excitation_absorption(f, 399.0) > 0.9 * excitation_absorption(f, 400.0)
            I, _  = quadgk(x -> excitation_absorption(f, x), 350.0, 400.0, 700.0; rtol = 1e-8)
            I0, _ = quadgk(x -> excitation_absorption(whole, x), 350.0, 400.0, 700.0; rtol = 1e-8)
            @test I ≈ I0 rtol = 1e-7
            Iref, _ = quadgk(x -> 0.02 * absorption(plain, x), 350.0, 400.0, 700.0; rtol = 1e-8)
            @test I ≈ Iref rtol = 1e-7
        end
        # Strict rule: no excitation outside the table, so the 350–700 nm integral is
        # that of the total absorption over 400–700 nm only.
        strict = Phytoplankton(Chl = 0.5, pigments = ProspectPigments(carotenoid_ratio = 0.3,
                                                                    outside = UnassignedOutsideTable()))
        fs = chlorophyll_fluorescence(strict; quantum_yield = 0.02, carotenoid_transfer = 1.0,
                                      excitation_range = (350.0, 700.0))
        @test excitation_absorption(fs, 399.999) == 0
        Is, _ = quadgk(x -> excitation_absorption(fs, x), 350.0, 400.0, 700.0; rtol = 1e-8)
        Ivis, _ = quadgk(x -> 0.02 * absorption(plain, x), 400.0, 700.0; rtol = 1e-8)
        @test Is ≈ Ivis rtol = 1e-7
    end

    @testset "Photon and energy kernels" begin
        phyto = Phytoplankton(Chl = 1.0)
        sif   = chlorophyll_fluorescence(phyto; quantum_yield = 0.02)
        λ′, λ = 440.0, 685.0
        # Quantum yield is part of the excitation coefficient, applied once.
        @test excitation_absorption(sif, λ′) ≈ quantum_yield(sif) * absorption(phyto, λ′)
        # Fluorescence converts photons: the energy kernel is φ a f λ'/λ, written out
        # from its parts (φ = 0.02, Fluspect emission density at 685 nm).
        @test conversion_units(sif) === PhotonConversion()
        @test radiance_source_kernel(sif, λ′, λ) ≈
              0.02 * absorption(phyto, λ′) * redistribution(chloroplast_emission(), λ′, λ) * 440 / 685
        @test inelastic_coefficient(sif, λ′, λ) ≈
              0.02 * absorption(phyto, λ′) * redistribution(chloroplast_emission(), λ′, λ)
        # Photon budget: emission integrates to the excitation coefficient.
        em = 600.0:0.25:900.0
        @test sum(inelastic_coefficient(sif, λ′, x) for x in em) * step(em) ≈
              excitation_absorption(sif, λ′) rtol = 1e-9
        # Zero yield, zero source.
        dark = chlorophyll_fluorescence(phyto; quantum_yield = 0.0)
        @test radiance_source_kernel(dark, λ′, λ) == 0
        @test isempty(excitation_pathways(WaterRaman{Float64}()))
    end
end
