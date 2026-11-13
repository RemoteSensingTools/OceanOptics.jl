@testset "Pigment partitions" begin

    @testset "Bundled spectra" begin
        # Bidigare et al. (1990) chlorophyll a: Soret maximum 0.0265 m²/mg at 440 nm.
        tbl = OceanOptics._BIDIGARE_SPECTRA[1]
        @test tbl.λ[[1, end]] == [400.0, 750.0] && length(tbl.λ) == 176
        @test maximum(tbl.y) == 0.0265 && tbl.λ[argmax(tbl.y)] == 440
        # PROSPECT-D: Kab extends to 780 nm, Kca to 557 nm.
        @test OceanOptics.linterp(OceanOptics._PROSPECT_KAB, 780.0) > 0
        @test OceanOptics.linterp(OceanOptics._PROSPECT_KAB, 781.0) == 0
        @test OceanOptics.linterp(OceanOptics._PROSPECT_KCA, 558.0) == 0
    end

    @testset "ProspectPigments" begin
        pig   = ProspectPigments(carotenoid_ratio = 0.3)
        phyto = Phytoplankton(Chl = 0.5, pigments = pig)
        @test absorption_parts(phyto) == (:chlorophyll, :carotenoids, :unassigned)
        @test pigment_class(pig, :carotenoids) == :carotenoid
        @test_throws ArgumentError pigment_class(pig, :fucoxanthin)
        for λ in (350.0, 420.0, 450.0, 490.0, 560.0, 676.0, 790.0)
            a = absorption_by_part(phyto, λ)
            @test sum(a) ≈ absorption(phyto, λ) rtol = 1e-12
            @test all(≥(0), a)
        end
        # Shares follow Kab : r·Kca.
        kab, kca = OceanOptics.linterp(OceanOptics._PROSPECT_KAB, 480.0),
                   OceanOptics.linterp(OceanOptics._PROSPECT_KCA, 480.0)
        f = pigment_fractions(pig, 480.0)
        @test f[2] / f[1] ≈ 0.3kca / kab rtol = 1e-12
        @test f[3] == 0
        # No carotenoid absorption beyond 557 nm; UV outside the table is unassigned.
        @test absorption_part(phyto, :carotenoids, 600.0) == 0
        @test absorption_part(phyto, :unassigned, 350.0) == absorption(phyto, 350.0)
        @test pigment_fractions(ProspectPigments(carotenoid_ratio = 0.0), 450.0) == (1.0, 0.0, 0.0)
        @test_throws ArgumentError ProspectPigments(carotenoid_ratio = -0.1)
        @test_throws ArgumentError absorption_part(phyto, :total, 450.0)
    end

    @testset "Bidigare1990Pigments" begin
        pig = Bidigare1990Pigments(chl_a = 0.5, chl_b = 0.05, chl_c = 0.1, psc = 0.3, ppc = 0.15)
        phyto = Phytoplankton(Chl = 0.5, pigments = pig)
        @test length(absorption_parts(phyto)) == 8 && last(absorption_parts(phyto)) == :unassigned
        @test pigment_class(pig, :photoprotective_carotenoids) == :carotenoid
        @test pigment_class(pig, :phycoerythrin_peb) == :phycobilin
        @test pigment_class(pig, :chlorophyll_c) == :chlorophyll
        for λ in (380.0, 440.0, 470.0, 550.0, 676.0, 740.0)
            @test sum(absorption_by_part(phyto, λ)) ≈ absorption(phyto, λ) rtol = 1e-12
        end
        # Shares are c_k a*_k / Σ c_j a*_j at a tabulated wavelength (Bidigare 1990 at 440 nm).
        tbl = OceanOptics._BIDIGARE_SPECTRA
        s = [c * t.y[findfirst(==(440.0), t.λ)]
             for (c, t) in zip(OceanOptics._concentrations(pig), tbl)]
        @test collect(pigment_fractions(pig, 440.0))[1:7] ≈ s ./ sum(s) rtol = 1e-12
        # Only ratios matter.
        @test all(pigment_fractions(Bidigare1990Pigments(chl_a = 5.0, chl_b = 0.5, chl_c = 1.0,
                                                         psc = 3.0, ppc = 1.5), 470.0) .≈
                  pigment_fractions(pig, 470.0))
        # Red band: chlorophyll a only; beyond it and outside the table: unassigned.
        @test pigment_fractions(Bidigare1990Pigments(chl_a = 1.0, psc = 1.0), 676.0)[1] == 1
        @test last(pigment_fractions(pig, 740.0)) == 1 && last(pigment_fractions(pig, 390.0)) == 1
        @test_throws ArgumentError Bidigare1990Pigments(chl_a = 0.0)
        @test_throws ArgumentError Bidigare1990Pigments(chl_a = 1.0, ppc = -0.1)
        @test_throws UndefKeywordError Bidigare1990Pigments(psc = 1.0)
    end

    @testset "Phytoplankton without a partition" begin
        phyto = Phytoplankton(Chl = 0.5)
        @test absorption_parts(phyto) == (:total,)
        @test absorption_part(phyto, :total, 450.0) == absorption(phyto, 450.0)
        @test absorption_parts(CDOM{Float64}()) == (:total,)
    end

    @testset "Partitions take the phytoplankton's precision" begin
        p32 = Phytoplankton{Float32, Bricaud1998}(Chl = 0.5, pigments = ProspectPigments(carotenoid_ratio = 0.3))
        @test p32.pigments isa ProspectPigments{Float32}
        @test absorption_by_part(p32, 450.0f0) isa NTuple{3, typeof(absorption(p32, 450.0f0))}
        b32 = Phytoplankton{Float32, Bricaud1998}(Chl = 0.5, pigments = Bidigare1990Pigments(chl_a = 0.5, ppc = 0.1))
        @test b32.pigments isa Bidigare1990Pigments{Float32}
    end

    @testset "AD through a partition" begin
        d = ForwardDiff.derivative(r -> absorption_part(
                Phytoplankton(Chl = 0.5, pigments = ProspectPigments(carotenoid_ratio = r)),
                :carotenoids, 480.0), 0.3)
        @test d > 0
        d = ForwardDiff.derivative(c -> absorption_part(
                Phytoplankton(Chl = 0.5, pigments = Bidigare1990Pigments(chl_a = 0.5, ppc = c)),
                :photoprotective_carotenoids, 470.0), 0.1)
        @test d > 0
    end
end
