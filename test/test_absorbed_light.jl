@testset "Absorbed light and emission" begin
    water = PureWater{Float64}()
    phyto = Phytoplankton(Chl = 1.0, pigments = ProspectPigments(carotenoid_ratio = 0.3))
    f     = chlorophyll_fluorescence(phyto; quantum_yield = 0.02, carotenoid_transfer = 0.5)
    layer = OceanLayer(5.0, 7.0, AbstractOceanConstituent{Float64}[water, phyto];
                       inelastic = AbstractOceanInelasticProcess{Float64}[f, WaterRaman{Float64}()])
    λ  = collect(400.0:2.0:750.0)
    o  = layer_optics(layer, λ; ℓ_max = 8)
    E₀ = [1.2 + 0.3sin(x / 40) for x in λ]          # W m⁻² nm⁻¹, layer mean

    @testset "Components and links" begin
        @test [(c.label, c.part) for c in o.absorbers] ==
              [("PureWater", :total), ("Phytoplankton", :chlorophyll),
               ("Phytoplankton", :carotenoids), ("Phytoplankton", :unassigned)]
        @test absorption_optical_depth(o) ≈ 2 .* (absorption(water, λ) .+ absorption(phyto, λ)) rtol = 1e-12
        c = o.inelastic[1]
        @test c.process === f && c.excitation == [2, 3, 4] && c.efficiency == [1.0, 0.5, 0.0]
        # The linked components reproduce the process's excitation coefficient.
        for k in (1, 40, 140)
            τI = quantum_yield(f) * sum(c.efficiency[i] * o.absorbers[c.excitation[i]].τ[k]
                                        for i in eachindex(c.excitation))
            @test τI / thickness(o) ≈ excitation_absorption(f, λ[k]) rtol = 1e-12
        end
        @test isempty(o.inelastic[2].excitation)
    end

    @testset "Absorbed power, photons and PAR" begin
        A = absorbed_power(o, E₀)
        @test size(A) == (4, length(λ))
        @test A[2, 10] ≈ o.absorbers[2].τ[10] * E₀[10]
        @test vec(sum(A; dims = 1)) ≈ absorption_optical_depth(o) .* E₀ rtol = 1e-12
        # 1 W at 550 nm is 4.598 µmol photons s⁻¹ (hc N_A = 0.119627 J m mol⁻¹).
        N = absorbed_photons(o, ones(length(λ)))
        k = findfirst(==(550.0), λ)
        @test N[1, k] / o.absorbers[1].τ[k] ≈ 550e-9 / 0.119626566 * 1e6 rtol = 1e-7
        par = absorbed_par(o, E₀)
        inpar = 400 .≤ λ .≤ 700
        Np = absorbed_photons(o, E₀)
        @test par[2] ≈ sum((Np[2, i] + Np[2, i+1]) / 2 * 2.0 for i in findall(inpar)[1:end-1]) rtol = 1e-12
        @test par[4] == 0                              # all absorption in 400–700 nm is attributed
        @test all(≥(0), par)
        # Decreasing grids give the same integrals.
        orev = layer_optics(layer, reverse(λ); ℓ_max = 8)
        @test absorbed_par(orev, reverse(E₀)) ≈ par rtol = 1e-12
        @test absorber_index(o, "Phytoplankton", :carotenoids) == 3
        @test absorber_index(o, "PureWater") == 1
        @test_throws ArgumentError absorber_index(o, "Phytoplankton")
        @test_throws DimensionMismatch absorbed_power(o, E₀[1:5])
        @test_throws ArgumentError absorbed_power(o, -E₀)
        @test_throws ArgumentError absorbed_par(o, E₀; range = (380.0, 700.0))
    end

    @testset "Integration windows" begin
        w(x, r) = OceanOptics._window_weights(x, r, "")
        # Linear integrands are integrated exactly over [lo, hi], also when the
        # grid has no points at the limits, and in either order.
        x = collect(398.0:5.0:703.0)
        @test sum(w(x, (400, 700))) ≈ 300 rtol = 1e-12
        @test sum(w(x, (400, 700)) .* x) ≈ (700^2 - 400^2) / 2 rtol = 1e-12
        @test w(reverse(x), (400, 700)) ≈ reverse(w(x, (400, 700)))
        @test sum(w(x, (401.5, 402.5))) ≈ 1 rtol = 1e-12        # inside a single interval
        @test w(collect(400.0:10.0:700.0), (400, 700))[[1, 2, end]] ≈ [5, 10, 5]
        @test_throws ArgumentError w(x, (390, 700))
        @test_throws ArgumentError w(x, (700, 400))
    end

    @testset "Excitation optical depth" begin
        c = o.inelastic[1]
        τI = excitation_optical_depth(o, 1)
        @test τI ≈ [thickness(o) * excitation_absorption(f, x) for x in λ] rtol = 1e-12
        τR = excitation_optical_depth(o, 2)                        # Raman: from the process
        @test τR ≈ [thickness(o) * excitation_absorption(o.inelastic[2].process, x) for x in λ]
    end

    @testset "Fluorescence photon budget" begin
        λ_em = collect(630.0:0.5:860.0)
        P    = inelastic_emission(o, 1, E₀, λ_em)
        @test all(≥(0), P)
        hc   = 6.62607015e-34 * 299_792_458.0
        emitted = sum(P[i] * λ_em[i] * 1e-9 / hc for i in eachindex(λ_em)) * 0.5
        # Photons delivered to chlorophyll in 400–700 nm (trapezoidal on the layer grid).
        Nd = absorbed_photons(o, E₀)
        c  = o.inelastic[1]
        w  = OceanOptics._window_weights(λ, (400.0, 700.0), "")
        delivered = sum(w[k] * sum(c.efficiency[i] * Nd[c.excitation[i], k]
                                   for i in eachindex(c.excitation)) for k in eachindex(λ))
        # µmol → photons: × N_A × 1e-6
        @test emitted ≈ quantum_yield(f) * delivered * 6.02214076e23 * 1e-6 rtol = 1e-6
        # Linear in the yield; zero yield emits nothing.
        f2 = chlorophyll_fluorescence(phyto; quantum_yield = 0.04, carotenoid_transfer = 0.5)
        o2 = layer_optics(OceanLayer(5.0, 7.0, AbstractOceanConstituent{Float64}[water, phyto];
                                     inelastic = AbstractOceanInelasticProcess{Float64}[f2]), λ; ℓ_max = 8)
        @test inelastic_emission(o2, 1, E₀, λ_em) ≈ 2P rtol = 1e-12
        f0 = chlorophyll_fluorescence(phyto; quantum_yield = 0.0, carotenoid_transfer = 0.5)
        o0 = layer_optics(OceanLayer(5.0, 7.0, AbstractOceanConstituent{Float64}[water, phyto];
                                     inelastic = AbstractOceanInelasticProcess{Float64}[f0]), λ; ℓ_max = 8)
        @test all(iszero, inelastic_emission(o0, 1, E₀, λ_em))
    end

    @testset "Excitation must attenuate" begin
        other = Phytoplankton(Chl = 2.0, pigments = ProspectPigments(carotenoid_ratio = 0.3))
        g = chlorophyll_fluorescence(other; quantum_yield = 0.02, carotenoid_transfer = 0.0)
        bad = OceanLayer(0.0, 1.0, AbstractOceanConstituent{Float64}[water, phyto];
                         inelastic = AbstractOceanInelasticProcess{Float64}[g])
        @test_throws ArgumentError layer_optics(bad, λ)
        # An equal (not identical) constituent is accepted.
        same = Phytoplankton(Chl = 1.0, pigments = ProspectPigments(carotenoid_ratio = 0.3))
        ok = OceanLayer(0.0, 1.0, AbstractOceanConstituent{Float64}[water, same];
                        inelastic = AbstractOceanInelasticProcess{Float64}[f])
        @test layer_optics(ok, λ; ℓ_max = 4).inelastic[1].excitation == [2, 3, 4]
        # Two equal constituents would make the excitation ambiguous.
        twice = OceanLayer(0.0, 1.0, AbstractOceanConstituent{Float64}[water, phyto, same];
                           inelastic = AbstractOceanInelasticProcess{Float64}[f])
        @test_throws ArgumentError layer_optics(twice, λ; ℓ_max = 4)
    end

    @testset "Column with a tabulated Lambertian floor" begin
        column = uniform_column(AbstractOceanConstituent{Float64}[water, phyto]; n_layers = 3,
                                depth = 15.0, inelastic = AbstractOceanInelasticProcess{Float64}[f])
        floor  = TabulatedLambertianBottom([350.0, 550.0, 800.0], [0.1, 0.3, 0.2]; source = "test")
        col = column_optics(column, λ; refractive_index = FixedRefractiveIndex(1.34),
                            bottom = floor, ℓ_max = 8)
        @test col.bottom isa LambertianBottom
        @test bottom_albedo(col)[findfirst(==(450.0), λ)] ≈ 0.2
        @test bottom_albedo(col)[findfirst(==(550.0), λ)] ≈ 0.3
        @test all(l -> length(l.inelastic) == 1, col)
        @test_throws DomainError column_optics(column, [340.0, 400.0];
                                               refractive_index = FixedRefractiveIndex(1.34),
                                               bottom = floor)
        @test_throws ArgumentError OceanColumnOptics(col.layers, col.refractive_index, floor)
        # Other bottom types pass through to the solver.
        struct TestBottom <: AbstractOceanBottom end
        @test column_optics(column, λ; refractive_index = FixedRefractiveIndex(1.34),
                            bottom = TestBottom(), ℓ_max = 4).bottom isa TestBottom
        @test_throws ArgumentError TabulatedLambertianBottom([400.0, 500.0], [0.1, 1.1])
        @test_throws ArgumentError TabulatedLambertianBottom([500.0, 400.0], [0.1, 0.1])
    end
end
