@testset "Tabulated inputs" begin

    @testset "TabulatedGreekPhase" begin
        # Rayleigh matrix (Hansen & Travis 1974), vSmartMOM convention.
        kw = (α = [0.0, 0.0, 3.0], γ = [0.0, 0.0, sqrt(1.5)], δ = [0.0, 1.5, 0.0],
              ϵ = zeros(3), ζ = zeros(3))
        pf = TabulatedGreekPhase([1.0, 0.0, 0.5]; kw..., source = "test")
        @test is_polarizable(pf) && !is_wavelength_dependent(pf)
        @test greek_expansion(pf, 4; polarized = true) ==
              greek_expansion(RayleighWaterPhase{Float64}(depolarization = 0.0), 4; polarized = true)
        @test phase_function_moments(pf, 1) == [1.0, 0.0]
        @test phase_function_moments(pf, 5) == [1.0, 0.0, 0.5, 0.0, 0.0, 0.0]
        @test !is_polarizable(TabulatedGreekPhase([1.0, 0.0, 0.5]))

        # Normalization is checked, and repaired only on explicit request.
        @test_throws ArgumentError TabulatedGreekPhase(1.0001 .* [1.0, 0.0, 0.5])
        c  = 1.0001
        pn = TabulatedGreekPhase(c .* [1.0, 0.0, 0.5]; (; (k => c .* v for (k, v) in pairs(kw))...)...,
                                 normalize = true)
        @test greek_expansion(pn, 2; polarized = true).γ ≈ kw.γ
        # A strongly peaked expansion truncated far before convergence
        # (HG g = 0.95 at ℓ = 33) has negative backscatter: rejected.
        @test_throws ArgumentError TabulatedGreekPhase(
            phase_function_moments(HenyeyGreensteinPhase{Float64}(g = 0.95), 33))
        # Moment-based diagnostics agree with the closed forms of the source.
        for g in (-0.4, 0.5, 0.85)
            hg = HenyeyGreensteinPhase{Float64}(g = g)
            tp = TabulatedGreekPhase(phase_function_moments(hg, 400))
            @test backscatter_fraction(tp) ≈ backscatter_fraction(hg) rtol = 1e-10
            @test asymmetry_parameter(tp) ≈ g atol = 1e-14
            @test phase_function_value(tp, -0.3) ≈ phase_function_value(hg, -0.3) rtol = 1e-10
        end
    end

    @testset "Exact-wavelength constituents" begin
        a = TabulatedAbsorber([350.0, 450.0], [0.02, 0.01]; source = "test")
        @test absorption(a, 450.0) == 0.01 && absorption(a, 450) == 0.01
        @test absorption(a, Float32(450)) == 0.01
        @test absorption(a, [350.0, 450.0]) == [0.02, 0.01]
        @test scattering(a, 450.0) == 0.0 && !is_total_absorption(a)
        @test_throws ArgumentError absorption(a, 400.0)
        @test_throws ArgumentError TabulatedAbsorber([450.0, 350.0], [0.01, 0.02])
        @test_throws ArgumentError TabulatedAbsorber([350.0, 450.0], [0.01, -0.02])
        @test_throws DimensionMismatch TabulatedAbsorber([350.0, 450.0], [0.01])
        grid = collect(400.0:0.1:401.0)
        @test absorption(TabulatedAbsorber(grid, fill(0.1, length(grid))), 400.3) == 0.1
        # 1e7/ν round trips match the tabulated wavelengths.
        @test absorption(a, 1e7 / (1e7 / 450.0)) == 0.01

        s = TabulatedScatterer([350.0, 450.0], [0.01, 0.005], HenyeyGreensteinPhase{Float64}(g = 0.9))
        @test absorption(s, 350.0) == 0.0
        @test backscattering(s, 350.0) ≈ 0.01 * backscatter_fraction(HenyeyGreensteinPhase{Float64}(g = 0.9))
        # A wavelength-dependent phase must cover every tabulated wavelength.
        pf = WavelengthDependentPhase([450.0, 650.0], [HenyeyGreensteinPhase{Float64}(g = g) for g in (0.3, 0.8)])
        @test_throws ArgumentError TabulatedScatterer([450.0, 550.0], [0.2, 0.1], pf)
    end

    @testset "Chowdhary et al. (2020) Table 2 ocean bodies" begin
        # Inputs and derived values transcribed from JQSRT 242, 106717,
        # Table 2 (p. 4), checked against the rendered table. a_blk is bulk
        # absorption (Eqs. 11a/11b: a_blk = a_p + a_y + a_w), water scattering
        # has δ_w = 0, and the body is Δz = 100 m thick. The printed ω_blk and
        # τ_blk follow from Eqs. (10) and (12); compared within half a unit of
        # their last printed digit.
        λ   = [350.0, 450.0, 550.0, 650.0]
        b_w = [0.0134, 0.0045, 0.0019, 0.0010]
        src = "Chowdhary et al. 2020, JQSRT 242, 106717, Table 2"
        water_scat = TabulatedScatterer(λ, b_w, RayleighWaterPhase(depolarization = 0.0); source = src)
        printed(x, decimals) = (x, 0.5 * 10.0^-decimals)
        cases = (
            (a_blk = [0.0204, 0.0092, 0.0565, 0.3400], b_p = zeros(4),
             ω = printed.([0.39644970, 0.32846715, 0.03253425, 0.002932551], [8, 8, 8, 9]),
             τ = printed.([3.38, 1.37, 5.84, 34.1], [2, 2, 2, 1])),
            (a_blk = [0.0215, 0.0144, 0.1065, 0.3787], b_p = [0.0422, 0.0335, 0.8050, 0.8050],
             ω = printed.([0.72114137, 0.72519084, 0.88340267, 0.68034101], 8),
             τ = printed.([7.71, 5.24, 91.34, 118.47], 2)),
        )
        for case in cases
            constituents = AbstractOceanConstituent{Float64}[
                TabulatedAbsorber(λ, case.a_blk; total = true, source = src), water_scat]
            if any(>(0), case.b_p)
                # τ and ϖ do not depend on the phase function; this stand-in is
                # NOT the AOS-IV hydrosol matrix (supplement Table S1).
                push!(constituents, TabulatedScatterer(λ, case.b_p, HenyeyGreensteinPhase(g = 0.95)))
            end
            o = layer_optics(OceanLayer(0.0, 100.0, constituents), λ; ℓ_max = 8)
            for i in eachindex(λ)
                @test abs(single_scattering_albedo(o)[i] - case.ω[i][1]) ≤ case.ω[i][2]
                @test abs(optical_depth(o)[i] - case.τ[i][1]) ≤ case.τ[i][2]
            end
        end
        # Pure-water body, polarized: the δ_w = 0 Rayleigh matrix needs no approximation.
        o = layer_optics(OceanLayer(0.0, 100.0, AbstractOceanConstituent{Float64}[
                TabulatedAbsorber(λ, cases[1].a_blk; total = true), water_scat]), λ;
                ℓ_max = 8, polarized = true)
        c = only(o.components)
        @test c.label == src && c.polarization == :full
        @test greek_expansion(c, 1).γ[3] ≈ sqrt(1.5) && greek_expansion(c, 1).δ[2] ≈ 1.5
        @test backscatter_fraction(o) ≈ fill(0.5, 4)
    end
end
