@testset "OceanLayer + layer_optics" begin

    @testset "Thickness and geometry" begin
        layer = OceanLayer(0.0, 5.0, AbstractOceanConstituent{Float64}[])
        @test thickness(layer)      == 5.0
        @test midpoint_depth(layer) == 2.5
        @test_throws ArgumentError OceanLayer(5.0, 3.0, AbstractOceanConstituent{Float64}[])
        @test_throws ArgumentError OceanLayer(-1.0, 2.0, AbstractOceanConstituent{Float64}[])
    end

    @testset "Single-constituent layer: pure seawater" begin
        water = PureWater{Float64}(scattering_model = Morel1974())
        layer = OceanLayer(0.0, 1.0, AbstractOceanConstituent{Float64}[water])
        λ     = [440.0, 500.0, 600.0]
        opt   = layer_optics(layer, λ; ℓ_max = 16)

        @test n_wavelengths(opt) == 3
        @test n_moments(opt)     == 17

        # Optical depth = (a + b) · Δz at each wavelength
        for (i, λi) in enumerate(λ)
            τ_ref = (absorption(water, λi) + scattering(water, λi)) * thickness(layer)
            @test opt.τ[i] ≈ τ_ref rtol=1e-12
        end

        # Single-scattering albedo consistency
        @test all(0.0 .≤ opt.ϖ .≤ 1.0)
        # Rayleigh water scatters forward/back symmetric → B = 0.5 exactly
        @test all(opt.B .≈ 0.5)

        # Legendre moments: β_0 is 1 everywhere by normalization
        @test all(opt.β[1, :] .≈ 1.0)
    end

    @testset "Multi-constituent layer mixes correctly" begin
        water = PureWater{Float64}(scattering_model = Morel1974())
        cdom  = CDOM{Float64}(a_ref = 0.1)
        layer = OceanLayer(0.0, 2.0, AbstractOceanConstituent{Float64}[water, cdom])
        λ     = [440.0, 550.0]
        opt   = layer_optics(layer, λ; ℓ_max = 8)

        # Absorption is water + cdom; CDOM contributes no scattering.
        for (i, λi) in enumerate(λ)
            a = absorption(water, λi) + absorption(cdom, λi)
            b = scattering(water, λi)
            τ_ref = (a + b) * thickness(layer)
            ϖ_ref = b / (a + b)
            @test opt.τ[i] ≈ τ_ref rtol=1e-12
            @test opt.ϖ[i] ≈ ϖ_ref rtol=1e-12
        end
    end

    @testset "Empty-scatterer layer yields unit isotropic phase moments" begin
        cdom  = CDOM{Float64}(a_ref = 0.1)
        layer = OceanLayer(0.0, 1.0, AbstractOceanConstituent{Float64}[cdom])
        opt   = layer_optics(layer, [500.0]; ℓ_max = 4)
        @test opt.ϖ[1] == 0.0
        @test opt.β[1, 1] ≈ 1.0
        @test all(opt.β[2:end, 1] .≈ 0.0)
    end

    @testset "Polarized layer_optics" begin
        # Rayleigh-dominated pure-water layer: α, γ, δ must be non-zero
        # at ℓ=1 (δ) and ℓ=2 (α, γ).
        water = PureWater{Float64}(scattering_model = Morel1974())
        layer = OceanLayer(0.0, 1.0, AbstractOceanConstituent{Float64}[water])
        opt   = layer_optics(layer, [500.0]; ℓ_max = 4, polarized = true)

        # β layout (always present). Default ρ = 0.039 depolarization gives
        # β₂ = 0.5·dpl_p ≈ 0.4713 — close to but not exactly 1/2.
        ρ     = 0.039
        dpl_p = (1 - ρ) / (1 + ρ/2)
        @test opt.β[1, 1] ≈ 1.0
        @test opt.β[3, 1] ≈ 0.5 * dpl_p rtol=1e-10

        # Polarized slots now populated
        @test !isempty(opt.α)
        @test size(opt.α) == size(opt.β)
        @test opt.α[3, 1] > 0                  # α₂ = 3·dpl_p > 0
        @test opt.γ[3, 1] > 0                  # γ₂ = √(3/2)·dpl_p > 0
        @test opt.δ[2, 1] > 0                  # δ₁ = 1.5·dpl_p·dpl_r > 0
        @test opt.ϵ[3, 1] == 0.0
        @test opt.ζ[3, 1] == 0.0

        # Scalar default still returns empty Greek slots
        opt_scalar = layer_optics(layer, [500.0]; ℓ_max = 4)
        @test isempty(opt_scalar.α)
    end

    @testset "OceanLayer with fluorophores (Phase 2)" begin
        cdom  = CDOM{Float64}(a_ref = 0.1)
        phyto = Phytoplankton{Float64, Bricaud1998}(Chl = 0.5)
        sif   = chlorophyll_fluorescence(phyto)
        ram   = WaterRaman{Float64}()

        # Default (no fluorophores) — backward compatible
        layer0 = OceanLayer(0.0, 5.0, AbstractOceanConstituent{Float64}[cdom])
        @test isempty(layer0.fluorophores)

        # With fluorophores
        layer = OceanLayer(0.0, 5.0,
                           AbstractOceanConstituent{Float64}[cdom];
                           fluorophores = AbstractOceanInelasticProcess{Float64}[sif, ram])
        @test length(layer.fluorophores) == 2
        @test layer.fluorophores[1] === sif
        @test layer.fluorophores[2] === ram

        # Layer-only fluorophores (no elastic constituents) — FT is inferred
        # from the fluorophore.
        layer_f = OceanLayer(0.0, 1.0,
                             AbstractOceanConstituent{Float64}[];
                             fluorophores = AbstractOceanInelasticProcess{Float64}[ram])
        @test layer_f isa OceanLayer{Float64}

        # Elastic layer_optics ignores fluorophores (two-pass solve is
        # downstream, per DESIGN §2.2.9)
        opt = layer_optics(layer, [500.0]; ℓ_max = 4)
        @test opt isa OceanLayerOptics
    end
end
