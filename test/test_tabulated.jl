@testset "Tabulated inputs" begin

    # Rayleigh phase matrix without depolarization, entered by hand from the
    # closed form (Hansen & Travis 1974): β₂ = 1/2, α₂ = 3, γ₂ = √(3/2), δ₁ = 3/2.
    rayleigh_greek() = (β = [1.0, 0.0, 0.5], α = [0.0, 0.0, 3.0],
                        γ = [0.0, 0.0, sqrt(1.5)], δ = [0.0, 1.5, 0.0],
                        ϵ = zeros(3), ζ = zeros(3))

    @testset "TabulatedGreekPhase construction" begin
        g = rayleigh_greek()
        pf = TabulatedGreekPhase(g.β; g.α, g.γ, g.δ, g.ϵ, g.ζ, source = "test")
        @test is_polarizable(pf) && !is_wavelength_dependent(pf)
        @test phase_matrix_moments(pf, 4) ==
              phase_matrix_moments(RayleighWaterPhase{Float64}(depolarization = 0.0), 4)
        # Truncation and exact zero padding.
        @test phase_function_moments(pf, 1) == [1.0, 0.0]
        @test phase_function_moments(pf, 5) == [1.0, 0.0, 0.5, 0.0, 0.0, 0.0]

        scalar = TabulatedGreekPhase(g.β)
        @test !is_polarizable(scalar)
        @test_throws ArgumentError phase_matrix_moments(scalar, 4)

        # All or none of the polarized coefficients.
        @test_throws ArgumentError TabulatedGreekPhase(g.β; α = g.α)
        @test_throws DimensionMismatch TabulatedGreekPhase(g.β; g.α, g.γ, g.δ, g.ϵ,
                                                           ζ = zeros(2))
        # Normalization is checked, and only repaired on explicit request.
        @test_throws ArgumentError TabulatedGreekPhase(1.0001 .* g.β)
        c = 1.0001
        pf_n = TabulatedGreekPhase(c .* g.β; α = c .* g.α, γ = c .* g.γ, δ = c .* g.δ,
                                   ϵ = c .* g.ϵ, ζ = c .* g.ζ, normalize = true)
        @test phase_function_moments(pf_n, 2) ≈ g.β
        @test phase_matrix_moments(pf_n, 2).γ ≈ g.γ
        # (2ℓ+1) applied twice is caught.
        @test_throws ArgumentError TabulatedGreekPhase([(2ℓ + 1)^2 * 0.8^ℓ for ℓ in 0:4])
        # A strongly forward-peaked expansion truncated far before convergence
        # (HG g = 0.95 at L = 33) has negative backscatter: rejected at once.
        hg95 = phase_function_moments(HenyeyGreensteinPhase{Float64}(g = 0.95), 33)
        @test_throws ArgumentError TabulatedGreekPhase(hg95)
        # Row-shifted polarized table (α₀ = 3, γ₁ = √(3/2)) is rejected.
        @test_throws ArgumentError TabulatedGreekPhase(g.β; α = [3.0, 0.0, 0.0],
                                                       γ = [0.0, sqrt(1.5), 0.0], g.δ, g.ϵ, g.ζ)
    end

    @testset "Backscatter fraction from moments" begin
        # Exact for a finite expansion; compare with independent closed forms.
        g = rayleigh_greek()
        @test backscatter_fraction(TabulatedGreekPhase(g.β)) ≈ 0.5 atol = 1e-15
        for gHG in (-0.4, 0.0, 0.5, 0.85)
            hg = HenyeyGreensteinPhase{Float64}(g = gHG)
            pf = TabulatedGreekPhase(phase_function_moments(hg, 400))
            @test backscatter_fraction(pf) ≈ backscatter_fraction(hg) rtol = 1e-10
            @test asymmetry_parameter(pf) ≈ gHG atol = 1e-14
            @test phase_function_value(pf, -0.3) ≈ phase_function_value(hg, -0.3) rtol = 1e-10
        end
    end

    @testset "Wavelength-dependent TabulatedGreekPhase" begin
        λ  = [450.0, 650.0]
        gs = [0.3, 0.8]
        β  = [(2ℓ + 1) * g^ℓ for ℓ in 0:200, g in gs]   # 0.8^200 ≈ 4e-20: converged
        pf = TabulatedGreekPhase(λ, β; source = "two HG columns")
        @test is_wavelength_dependent(pf) && !is_polarizable(pf)
        @test phase_function_moments(pf, 2, 650.0) ≈ [1.0, 2.4, 3.2]
        @test backscatter_fraction(pf, 450.0) ≈
              backscatter_fraction(HenyeyGreensteinPhase{Float64}(g = 0.3)) rtol = 1e-12
        # No interpolation, and no λ-free query of a λ-dependent table.
        @test_throws ArgumentError phase_function_moments(pf, 2, 550.0)
        @test_throws ArgumentError phase_function_moments(pf, 2)
        @test_throws ArgumentError TabulatedGreekPhase([650.0, 450.0], β)

        # Through layer_optics: each wavelength uses its own column.
        s     = TabulatedScatterer(λ, [0.2, 0.1], pf)
        layer = OceanLayer(0.0, 1.0, AbstractOceanConstituent{Float64}[s])
        opt   = layer_optics(layer, λ; ℓ_max = 8)
        @test opt.β[2, :] ≈ 3 .* gs
        @test opt.B ≈ [backscatter_fraction(HenyeyGreensteinPhase{Float64}(g = g))
                       for g in gs] rtol = 1e-8
        @test opt.ϖ == [1.0, 1.0]
        # A scatterer cannot cover wavelengths its phase table lacks.
        @test_throws ArgumentError TabulatedScatterer([450.0, 550.0], [0.2, 0.1], pf)

        # Decreasing output grid and Float32 output with a per-λ phase table.
        opt_r = layer_optics(layer, reverse(λ); ℓ_max = 8)
        @test opt_r.β[2, :] ≈ 3 .* reverse(gs)
        opt_32 = layer_optics(layer, Float32.(λ); ℓ_max = 8, FT = Float32)
        @test opt_32 isa OceanLayerOptics{Float32}
        @test opt_32.β[2, :] ≈ Float32.(3 .* gs)
    end

    @testset "Exact-wavelength tabulated constituents" begin
        a = TabulatedAbsorber([350.0, 450.0], [0.02, 0.01]; source = "test")
        @test absorption(a, 450.0) == 0.01
        @test absorption(a, 450) == 0.01
        @test absorption(a, Float32(450)) == 0.01
        @test absorption(a, [350.0, 450.0]) == [0.02, 0.01]
        @test scattering(a, 450.0) == 0.0
        @test_throws ArgumentError absorption(a, 400.0)
        @test_throws ArgumentError TabulatedAbsorber([450.0, 350.0], [0.01, 0.02])
        @test_throws ArgumentError TabulatedAbsorber([350.0, 450.0], [0.01, -0.02])
        @test_throws DimensionMismatch TabulatedAbsorber([350.0, 450.0], [0.01])

        # Range-generated grids still match their decimal literals.
        grid = collect(400.0:0.1:401.0)
        t    = TabulatedAbsorber(grid, fill(0.1, length(grid)))
        @test absorption(t, 400.3) == 0.1

        s = TabulatedScatterer([350.0, 450.0], [0.01, 0.005],
                               HenyeyGreensteinPhase{Float64}(g = 0.9))
        @test absorption(s, 350.0) == 0.0
        @test backscattering(s, 350.0) ≈
              0.01 * backscatter_fraction(HenyeyGreensteinPhase{Float64}(g = 0.9))
    end

    @testset "Total absorption cannot be double counted" begin
        λ     = [350.0, 450.0]
        a_blk = TabulatedAbsorber(λ, [0.02, 0.01]; total = true)
        water = PureWater{Float64}()
        scat  = TabulatedScatterer(λ, [0.01, 0.005], RayleighWaterPhase{Float64}())
        @test is_total_absorption(a_blk)
        @test !is_total_absorption(water)

        ok = OceanLayer(0.0, 1.0, AbstractOceanConstituent{Float64}[a_blk, scat])
        @test layer_optics(ok, λ; ℓ_max = 4).τ ≈ [0.03, 0.015]

        twice = OceanLayer(0.0, 1.0, AbstractOceanConstituent{Float64}[a_blk, water])
        @test_throws ArgumentError layer_optics(twice, λ; ℓ_max = 4)
        two_totals = OceanLayer(0.0, 1.0, AbstractOceanConstituent{Float64}[a_blk, a_blk])
        @test_throws ArgumentError layer_optics(two_totals, λ; ℓ_max = 4)
        # A constituent-only absorber may be combined with anything.
        part = TabulatedAbsorber(λ, [0.02, 0.01])
        both = OceanLayer(0.0, 1.0, AbstractOceanConstituent{Float64}[part, water])
        @test layer_optics(both, λ; ℓ_max = 4) isa OceanLayerOptics
    end

    @testset "Chowdhary et al. (2020) Table 2 ocean bodies" begin
        # Inputs and derived values transcribed from the JQSRT 242, 106717
        # article, Table 2 (p. 4), and checked against the rendered table.
        # a_blk is bulk absorption (Eqs. 11a/11b: a_blk = a_p + a_y + a_w),
        # water scattering has δ_w = 0, and the body is Δz = 100 m thick.
        # The printed ω_blk and τ_blk follow from Eqs. (10) and (12); they are
        # compared to within half a unit of their last printed digit.
        λ   = [350.0, 450.0, 550.0, 650.0]
        b_w = [0.0134, 0.0045, 0.0019, 0.0010]
        src = "Chowdhary et al. 2020, JQSRT 242, 106717, Table 2"
        water_scat = TabulatedScatterer(λ, b_w, RayleighWaterPhase(depolarization = 0.0);
                                        source = src)
        printed(x, decimals) = (x, 0.5 * 10.0^-decimals)

        cases = (
            (name   = "AOS-II/III",
             a_blk  = [0.0204, 0.0092, 0.0565, 0.3400], b_p = zeros(4),
             ω      = printed.([0.39644970, 0.32846715, 0.03253425, 0.002932551], [8, 8, 8, 9]),
             τ      = printed.([3.38, 1.37, 5.84, 34.1], [2, 2, 2, 1])),
            (name   = "AOS-IV",
             a_blk  = [0.0215, 0.0144, 0.1065, 0.3787], b_p = [0.0422, 0.0335, 0.8050, 0.8050],
             ω      = printed.([0.72114137, 0.72519084, 0.88340267, 0.68034101], 8),
             τ      = printed.([7.71, 5.24, 91.34, 118.47], 2)),
        )
        for case in cases
            constituents = AbstractOceanConstituent{Float64}[
                TabulatedAbsorber(λ, case.a_blk; total = true, source = src), water_scat]
            if any(>(0), case.b_p)
                # τ and ϖ do not depend on the phase function. This stand-in
                # is NOT the AOS-IV hydrosol matrix (supplement Table S1).
                push!(constituents, TabulatedScatterer(λ, case.b_p,
                                                       HenyeyGreensteinPhase(g = 0.95)))
            end
            opt = layer_optics(OceanLayer(0.0, 100.0, constituents), λ; ℓ_max = 8)
            for i in eachindex(λ)
                @test abs(opt.ϖ[i] - case.ω[i][1]) ≤ case.ω[i][2]
                @test abs(opt.τ[i] - case.τ[i][1]) ≤ case.τ[i][2]
            end
        end

        # Pure-water body, polarized: δ_w = 0 Rayleigh matrix, no approximation needed.
        body = OceanLayer(0.0, 100.0, AbstractOceanConstituent{Float64}[
            TabulatedAbsorber(λ, cases[1].a_blk; total = true), water_scat])
        opt  = layer_optics(body, λ; ℓ_max = 8, polarized = true)
        @test opt.α[3, :] ≈ fill(3.0, 4)
        @test opt.γ[3, :] ≈ fill(sqrt(1.5), 4)
        @test opt.δ[2, :] ≈ fill(1.5, 4)
        @test opt.B ≈ fill(0.5, 4)

        # Adding PureWater() to the bulk absorption would double count water.
        @test_throws ArgumentError layer_optics(
            OceanLayer(0.0, 100.0, AbstractOceanConstituent{Float64}[
                TabulatedAbsorber(λ, cases[1].a_blk; total = true), water_scat,
                PureWater{Float64}()]), λ; ℓ_max = 8)
    end
end
