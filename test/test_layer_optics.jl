@testset "Layer and column optics" begin

    water = PureWater{Float64}(depolarization = 0.0)
    phyto = Phytoplankton{Float64, Bricaud1998}(Chl = 1.0)
    cdom  = CDOM{Float64}(a_ref = 0.05)
    λ     = [412.5, 443.0, 490.0, 555.0, 670.0]

    @testset "Components and mixed quantities" begin
        layer = OceanLayer(2.0, 7.0, AbstractOceanConstituent{Float64}[water, cdom, phyto])
        o = layer_optics(layer, λ; ℓ_max = 32)
        @test o isa OceanLayerOptics{Float64} && !is_polarized(o)
        @test thickness(o) == 5 && o.depth_top == 2 && n_moments(o) == 33
        @test [c.label for c in o.components] == ["PureWater", "Phytoplankton"]
        for (i, λi) in enumerate(λ)
            a = absorption(water, λi) + absorption(cdom, λi) + absorption(phyto, λi)
            bw, bp = scattering(water, λi), scattering(phyto, λi)
            @test absorption_optical_depth(o)[i] ≈ 5a rtol = 1e-12
            @test o.components[1].τ[i] ≈ 5bw rtol = 1e-12
            @test o.components[2].τ[i] ≈ 5bp rtol = 1e-12
            @test optical_depth(o)[i] ≈ 5(a + bw + bp) rtol = 1e-12
            @test single_scattering_albedo(o)[i] ≈ (bw + bp) / (a + bw + bp) rtol = 1e-12
            # Mixed expansion: scattering-weighted, so β₁ comes from the particles only.
            g = greek_expansion(o, i)
            @test g.β[1] ≈ 1
            @test g.β[2] ≈ bp / (bw + bp) * phase_function_moments(scattering_phase(phyto), 1)[2] rtol = 1e-12
            # Backscatter fraction from the phase functions, not from the truncated expansion.
            @test backscatter_fraction(o)[i] ≈ (0.5bw + bp * backscatter_fraction(scattering_phase(phyto))) /
                                               (bw + bp) rtol = 1e-12
        end
        # Wavelength-independent phase functions are expanded once.
        @test all(!is_wavelength_dependent(c) for c in o.components)
        @test all(c.polarization == :scalar for c in o.components)
    end

    @testset "Polarization policy and labels" begin
        layer = OceanLayer(0.0, 1.0, AbstractOceanConstituent{Float64}[water, phyto])
        @test_throws ArgumentError layer_optics(layer, λ; ℓ_max = 8, polarized = true)
        @test_throws ArgumentError layer_optics(layer, λ; ℓ_max = 8, polarized = true,
                                                scalar_phase_approximation = :rayleigh)
        o = layer_optics(layer, λ; ℓ_max = 8, polarized = true,
                         scalar_phase_approximation = :depolarizing)
        @test is_polarized(o)
        @test [c.polarization for c in o.components] == [:full, :depolarizing]
        # Particles contribute only to β: the mixed α₂ is the water share of 3.
        for i in eachindex(λ)
            w = o.components[1].τ[i] / scattering_optical_depth(o)[i]
            g = greek_expansion(o, i)
            @test g.α[3] ≈ 3w rtol = 1e-12
            @test g.γ[3] ≈ sqrt(1.5) * w rtol = 1e-12
        end
        # Polarization does not change the scalar part.
        os = layer_optics(layer, λ; ℓ_max = 8)
        @test all(greek_expansion(o, i).β ≈ greek_expansion(os, i).β for i in eachindex(λ))
        # A pure-water layer needs no approximation.
        ow = layer_optics(OceanLayer(0.0, 1.0, AbstractOceanConstituent{Float64}[water]), λ;
                          ℓ_max = 4, polarized = true)
        @test only(ow.components).polarization == :full
    end

    @testset "Zero scattering and empty layers" begin
        a = TabulatedAbsorber([450.0, 550.0], [0.2, 0.3])
        s = TabulatedScatterer([450.0, 550.0], [0.0, 0.1], RayleighWaterPhase{Float64}(depolarization = 0.0))
        o = layer_optics(OceanLayer(0.0, 1.0, AbstractOceanConstituent{Float64}[a, s]),
                         [450.0, 550.0]; ℓ_max = 4, polarized = true)
        @test single_scattering_albedo(o)[1] == 0 && backscatter_fraction(o)[1] == 0
        g = greek_expansion(o, 1)
        @test g.β == [1.0, 0, 0, 0, 0] && all(iszero, (g.α, g.γ, g.δ, g.ϵ, g.ζ))
        @test greek_expansion(o, 2).α[3] ≈ 3
        e = layer_optics(OceanLayer(0.0, 1.0, AbstractOceanConstituent{Float64}[]), [500.0]; ℓ_max = 2)
        @test isempty(e.components) && optical_depth(e) == [0.0]
        @test greek_expansion(e, 1).β == [1.0, 0.0, 0.0]
    end

    @testset "Wavelength grids and precision" begin
        layer = OceanLayer(0.0, 2.0, AbstractOceanConstituent{Float64}[water, phyto])
        o  = layer_optics(layer, λ; ℓ_max = 8)
        or = layer_optics(layer, reverse(λ); ℓ_max = 8)
        @test optical_depth(or) ≈ reverse(optical_depth(o))
        @test_throws ArgumentError layer_optics(layer, [550.0, 450.0, 500.0]; ℓ_max = 8)
        @test_throws ArgumentError layer_optics(layer, [450.0, 450.0]; ℓ_max = 8)
        # vSmartMOM-style wavenumber grid, evaluated at λ = 1e7/ν.
        ν = collect(range(1e7 / 670, 1e7 / 412.5; length = 7))
        @test n_wavelengths(layer_optics(layer, 1e7 ./ ν; ℓ_max = 8)) == 7
        o32 = layer_optics(layer, λ; ℓ_max = 8, FT = Float32)
        @test o32 isa OceanLayerOptics{Float32}
        @test optical_depth(o32) ≈ Float32.(optical_depth(o)) rtol = 1e-5
        # Outside a constituent's tabulated range: an error, never an extrapolation.
        @test_throws DomainError layer_optics(layer, [250.0]; ℓ_max = 2)
    end

    @testset "Wavelength-dependent phase matrices" begin
        λ2 = [450.0, 650.0]
        pf = WavelengthDependentPhase(λ2, [RayleighWaterPhase{Float64}(depolarization = ρ)
                                           for ρ in (0.0, 0.09)])
        s  = TabulatedScatterer(λ2, [0.1, 0.1], pf)
        o  = layer_optics(OceanLayer(0.0, 1.0, AbstractOceanConstituent{Float64}[s]), λ2;
                          ℓ_max = 4, polarized = true)
        c  = only(o.components)
        @test is_wavelength_dependent(c) && length(c.greek) == 2 && c.phase_index == [1, 2]
        @test greek_expansion(c, 1).γ[3] ≈ sqrt(1.5)
        @test greek_expansion(c, 2).γ[3] < sqrt(1.5)
        @test c.label == "TabulatedScatterer"

        # A fine solver grid with band-representative phase matrices: each
        # distinct expansion is computed and stored once.
        pn  = WavelengthDependentPhase(λ2, pf.phases; lookup = :nearest)
        λf  = collect(range(440.0, 660.0; length = 2000))
        sf  = TabulatedScatterer(λf, fill(0.1, 2000), pn)
        t   = @elapsed of = layer_optics(OceanLayer(0.0, 1.0, AbstractOceanConstituent{Float64}[sf]), λf;
                                         ℓ_max = 64, polarized = true)
        cf  = only(of.components)
        @test length(cf.greek) == 2 && length(cf.phase_index) == 2000
        @test greek_expansion(cf, 1) === cf.greek[1] && greek_expansion(cf, 2000) === cf.greek[2]
        @test cf.phase_index[findfirst(>(550.0), λf)] == 2
        @test t < 10
        @test_throws ArgumentError WavelengthDependentPhase(λ2, pf.phases; lookup = :linear)
    end

    @testset "Total absorption cannot be double counted" begin
        λ2 = [350.0, 450.0]
        bulk = TabulatedAbsorber(λ2, [0.02, 0.01]; total = true)
        scat = TabulatedScatterer(λ2, [0.01, 0.005], RayleighWaterPhase{Float64}())
        @test optical_depth(layer_optics(OceanLayer(0.0, 1.0, AbstractOceanConstituent{Float64}[bulk, scat]),
                                         λ2; ℓ_max = 2)) ≈ [0.03, 0.015]
        @test_throws ArgumentError layer_optics(OceanLayer(0.0, 1.0,
            AbstractOceanConstituent{Float64}[bulk, PureWater{Float64}()]), λ2; ℓ_max = 2)
        @test_throws ArgumentError layer_optics(OceanLayer(0.0, 1.0,
            AbstractOceanConstituent{Float64}[bulk, bulk]), λ2; ℓ_max = 2)
    end

    @testset "Constructor validation" begin
        g  = GreekCoefficients([1.0, 0.0, 0.5])
        sc(; τ = [0.1, 0.2], greek = [g], B = [0.5], pol = :scalar, label = "x", index = Int[]) =
            ScatteringComponent(label, τ, greek, B, pol; phase_index = index)
        @test sc() isa ScatteringComponent
        @test_throws ArgumentError sc(τ = [-0.1, 0.2])
        @test_throws ArgumentError sc(pol = :full)                # scalar expansion labelled :full
        @test_throws ArgumentError sc(pol = :maybe)
        @test_throws DimensionMismatch sc(B = [0.5, 0.5])
        @test_throws ArgumentError sc(B = [1.5])
        # Several expansions need a per-wavelength index into them.
        @test_throws ArgumentError sc(greek = [g, g], B = [0.5, 0.5])
        @test_throws DimensionMismatch sc(greek = [g, g], B = [0.5, 0.5], index = [1, 2, 1])
        @test_throws ArgumentError sc(greek = [g, g], B = [0.5, 0.5], index = [1, 3])
        @test is_wavelength_dependent(sc(greek = [g, g], B = [0.5, 0.5], index = [2, 1]))
        lay(; λ = [450.0, 550.0], top = 0.0, bot = 1.0, τa = [0.1, 0.1], comps = [sc()], n = 3, pol = false) =
            OceanLayerOptics{Float64}(λ, top, bot, τa, comps, n, pol)
        @test lay() isa OceanLayerOptics
        @test_throws ArgumentError lay(λ = [450.0, 450.0])
        @test_throws ArgumentError lay(top = 1.0)
        @test_throws DimensionMismatch lay(τa = [0.1])
        @test_throws ArgumentError lay(n = 4)
        @test_throws ArgumentError lay(pol = true)
        @test_throws DimensionMismatch lay(comps = [sc(τ = [0.1])])
    end

    @testset "Column optics" begin
        column = uniform_column(AbstractOceanConstituent{Float64}[water]; n_layers = 3, depth = 30.0)
        λ2 = [450.0, 550.0]
        col = column_optics(column, λ2; refractive_index = FixedRefractiveIndex(1.34),
                            bottom = BlackBottom(), ℓ_max = 8, polarized = true)
        @test col isa OceanColumnOptics && n_layers(col) == 3 && bottom_depth(col) == 30
        @test is_polarized(col) && n_moments(col) == 9 && wavelengths(col) == λ2
        @test col.refractive_index == [1.34, 1.34] && bottom_albedo(col) == [0.0, 0.0]
        @test [l.depth_top for l in col] == [0.0, 10.0, 20.0]
        @test interface_depths(col) == [0.0, 10.0, 20.0, 30.0]
        # Equal phase functions share one phase_id and the same expansion objects.
        @test allequal(only(l.components).phase_id for l in col)
        @test all(only(l.components).greek[1] === only(col[1].components).greek[1] for l in col)
        colQ = column_optics(column, λ2; refractive_index = QuanFry1995(temperature = 20, salinity = 35),
                             bottom = LambertianBottom([0.1, 0.2]), ℓ_max = 8)
        @test colQ.refractive_index ≈ refractive_index(QuanFry1995(temperature = 20, salinity = 35), λ2)
        @test bottom_albedo(colQ) == [0.1, 0.2]
        @test bottom_albedo(column_optics(column, λ2; refractive_index = FixedRefractiveIndex(1.34),
                                          bottom = LambertianBottom(0.3))) == [0.3, 0.3]
        # Required choices and validity domains.
        @test_throws UndefKeywordError column_optics(column, λ2; bottom = BlackBottom())
        @test_throws UndefKeywordError column_optics(column, λ2; refractive_index = FixedRefractiveIndex(1.34))
        @test_throws DomainError column_optics(column, [750.0]; bottom = BlackBottom(),
                                               refractive_index = QuanFry1995(temperature = 20, salinity = 35))
        @test_throws DimensionMismatch column_optics(column, λ2; refractive_index = FixedRefractiveIndex(1.34),
                                                     bottom = LambertianBottom([0.1, 0.2, 0.3]))
        @test_throws ArgumentError LambertianBottom(1.2)
        l1, l2 = col.layers[1], col.layers[2]
        @test_throws ArgumentError OceanColumnOptics([l2, l1], [1.34, 1.34], BlackBottom())
        @test_throws ArgumentError OceanColumnOptics([l1, col.layers[3]], [1.34, 1.34], BlackBottom())
    end

    @testset "Fluorophores are not elastic optics" begin
        layer = OceanLayer(0.0, 5.0, AbstractOceanConstituent{Float64}[cdom];
                           fluorophores = AbstractOceanInelasticProcess{Float64}[WaterRaman{Float64}()])
        @test isempty(layer_optics(layer, [500.0]; ℓ_max = 4).components)
    end
end
