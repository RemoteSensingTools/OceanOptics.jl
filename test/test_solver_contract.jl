# What a radiative-transfer solver relies on when it consumes the layer and
# column optics: expansion keys, wavelength matching, irradiance round-off,
# bottoms, the interface index and the layout of scalar expansions.

struct UnknownTestBottom <: AbstractOceanBottom end

@testset "Solver contract" begin
    @testset "phase_id: content key of an expansion" begin
        # Published FNV-1a 64-bit test vectors (Fowler, Noll & Vo reference suite).
        @test OceanOptics._fnv1a(OceanOptics._FNV_OFFSET, "")       == 0xcbf29ce484222325
        @test OceanOptics._fnv1a(OceanOptics._FNV_OFFSET, "a")      == 0xaf63dc4c8601ec8c
        @test OceanOptics._fnv1a(OceanOptics._FNV_OFFSET, "foobar") == 0x85944171f73967e8

        # Keys computed outside Julia (Python `struct.pack`, byte-wise FNV-1a) over
        # the type name, then length and little-endian values of α, β, γ, δ, ϵ, ζ:
        # fixed numbers, so the key cannot depend on the session or Julia version.
        g = GreekCoefficients([1.0, 0.0, 0.5])
        @test phase_id(g) == 0x5be608ecb86bcfb4
        @test phase_id(GreekCoefficients(Float32[1, 0, 0.5])) == 0x10904bc3da6e3ce3
        @test phase_id(depolarized(g)) == 0xdc568201c219d257

        # Equal content in distinct arrays ⇒ equal key; −0.0 == 0.0.
        @test phase_id(GreekCoefficients([1.0, 0.0, 0.5])) == phase_id(g)
        @test phase_id(GreekCoefficients([1.0, -0.0, 0.5])) == phase_id(g)
        # Any difference in a coefficient, ℓ_max, element type or polarization ⇒ another key.
        others = (GreekCoefficients([1.0, 0.0, 0.5 + 1e-15]), GreekCoefficients([1.0, 0.0, 0.5, 0.0]),
                  GreekCoefficients(Float32[1, 0, 0.5]), depolarized(g))
        @test allunique([phase_id(g); [phase_id(x) for x in others]])

        λ     = [450.0, 550.0, 650.0]
        pf    = HenyeyGreensteinPhase{Float64}(g = 0.8)
        layer(pf, FT = Float64) = OceanLayer(FT(0), FT(1),
            AbstractOceanConstituent{FT}[TabulatedScatterer(FT.(λ), FT[0.3, 0.2, 0.1], pf)])
        sc(o) = only(o.scatterers)
        o16 = layer_optics(layer(pf), λ; ℓ_max = 16)
        @test sc(o16).phase_id == [phase_id(only(sc(o16).greek))]
        @test sc(o16).phase_id == sc(layer_optics(layer(HenyeyGreensteinPhase{Float64}(g = 0.8)), λ; ℓ_max = 16)).phase_id
        @test sc(o16).phase_id != sc(layer_optics(layer(pf), λ; ℓ_max = 32)).phase_id
        @test sc(o16).phase_id != sc(layer_optics(layer(pf), λ; ℓ_max = 16, polarized = true,
                                                  scalar_phase_approximation = :depolarizing)).phase_id
        pf32 = HenyeyGreensteinPhase{Float32}(g = 0.8f0)
        @test sc(o16).phase_id != sc(layer_optics(layer(pf32, Float32), Float32.(λ); ℓ_max = 16)).phase_id

        # Equal-valued tables held in different arrays share the key.
        β = [1.0, 1.5, 1.25, 0.875]                              # (2ℓ+1) gˡ, g = 1/2
        tab() = sc(layer_optics(layer(TabulatedGreekPhase(copy(β))), λ; ℓ_max = 3)).phase_id
        @test tab() == tab() == [phase_id(GreekCoefficients(β))]

        # A wavelength-dependent phase function: one key per distinct expansion.
        wd = WavelengthDependentPhase([450.0, 650.0],
                                      [HenyeyGreensteinPhase{Float64}(g = x) for x in (0.8, 0.9)];
                                      lookup = :nearest)
        c = sc(layer_optics(layer(wd), λ; ℓ_max = 16))
        @test length(c.phase_id) == 2 && allunique(c.phase_id)
        @test c.phase_id[1] == only(sc(o16).phase_id)            # the g = 0.8 entry
        @test [phase_id(c, k) for k in eachindex(λ)] == [phase_id(greek_expansion(c, k)) for k in eachindex(λ)]
        @test phase_id(c, 1) != phase_id(c, 3)
    end

    @testset "Tabulated constituents: wavelength match and interpolation" begin
        λ_tab = [350.0, 412.0, 443.0, 475.0, 550.0, 670.0, 865.0]
        a_tab = [0.07, 0.05, 0.04, 0.03, 0.06, 0.44, 4.6]
        a     = TabulatedAbsorber(λ_tab, a_tab)
        # Single-precision wavelengths obtained from single-precision wavenumbers.
        λ32 = 1f7 ./ (1f7 ./ Float32.(λ_tab))
        @test any(Float64.(λ32) .!= λ_tab)                       # the round trip is not exact
        @test absorption(a, λ32) == a_tab
        @test absorption(a, Float64.(λ32)) == a_tab
        @test absorption(a, 1e7 ./ (1e7 ./ λ_tab)) == a_tab
        o = layer_optics(OceanLayer(0.0, 2.0, AbstractOceanConstituent{Float64}[a]), λ32; ℓ_max = 2)
        @test only(o.absorbers).τ ≈ 2 .* a_tab rtol = 1e-15
        # The tolerance is 1e-6 relative: 2e-6 away is another wavelength.
        @test absorption(a, 550.0 * (1 + 0.9e-6)) == 0.06
        @test_throws "match tolerance 1.0e-6 relative" absorption(a, 550.0 * (1 + 2e-6))
        @test_throws "is not a tabulated wavelength" absorption(a, 500.0)
        @test_throws "neighbouring wavelengths must differ" TabulatedAbsorber(
            [550.0, 550.0 * (1 + 1e-6)], [0.1, 0.1])

        lin = TabulatedAbsorber(λ_tab, a_tab; lookup = :linear)
        @test absorption(lin, λ_tab) == a_tab
        @test absorption(lin, 512.5) ≈ (0.03 + 0.06) / 2 rtol = 1e-14       # midpoint of 475 and 550
        @test absorption(lin, 400.0) ≈ 0.07 + (0.05 - 0.07) * 50 / 62 rtol = 1e-14
        @test_throws DomainError absorption(lin, 340.0)
        @test_throws DomainError absorption(lin, 900.0)
        @test_throws "lookup must be :exact or :linear, got :nearest" TabulatedAbsorber(
            λ_tab, a_tab; lookup = :nearest)
        # The lookup is a type: names and lookup objects select the same one.
        @test a isa TabulatedAbsorber{Float64, ExactLookup}
        @test lin isa TabulatedAbsorber{Float64, LinearLookup}
        @test TabulatedAbsorber(λ_tab, a_tab; lookup = LinearLookup()) isa TabulatedAbsorber{Float64, LinearLookup}
        @test (@inferred absorption(lin, 512.5)) isa Float64
        # Slope of the segment 400–500 nm of a table (1, 2): 0.01 per nm.
        ramp = TabulatedAbsorber([400.0, 500.0], [1.0, 2.0]; lookup = :linear)
        @test ForwardDiff.derivative(x -> absorption(ramp, x), 450.0) ≈ 0.01 rtol = 1e-14
        # A 1 cm⁻¹ wavenumber grid inside the table runs through layer_optics.
        λν = 1e7 ./ (14_300.0:1.0:14_310.0)                       # 698.8–699.3 nm, decreasing
        τ  = only(layer_optics(OceanLayer(0.0, 1.0, AbstractOceanConstituent{Float64}[lin]), λν;
                               ℓ_max = 2).absorbers).τ
        @test τ ≈ [0.44 + (4.6 - 0.44) * (x - 670) / (865 - 670) for x in λν] rtol = 1e-13

        pf = RayleighWaterPhase{Float64}()
        s  = TabulatedScatterer([400.0, 500.0], [0.02, 0.01], pf; lookup = :linear)
        @test scattering(s, 425.0) ≈ 0.0175 rtol = 1e-14
        @test backscattering(s, 425.0) ≈ 0.0175 / 2 rtol = 1e-12  # Rayleigh: B = 1/2
        @test_throws "is not a tabulated wavelength" scattering(
            TabulatedScatterer([400.0, 500.0], [0.02, 0.01], pf), 425.0)
        @test_throws "lookup must be :exact or :linear, got :cubic" TabulatedScatterer(
            [400.0, 500.0], [0.02, 0.01], pf; lookup = :cubic)
        # Interpolated values have the type the query and the table promote to,
        # on the matching and on the interpolating branch.
        s32 = TabulatedScatterer(Float32[400, 500], Float32[0.02, 0.01],
                                 RayleighWaterPhase{Float32}(); lookup = :linear)
        @test scattering(s32, 425.0) isa Float64 && scattering(s32, 400.0) isa Float64
        @test scattering(s32, 425f0) isa Float32 && scattering(s32, 400f0) isa Float32
    end

    @testset "Irradiance with round-off negatives" begin
        λ = collect(400.0:50.0:700.0)
        o = layer_optics(OceanLayer(0.0, 3.0, AbstractOceanConstituent{Float64}[PureWater{Float64}(), CDOM{Float64}()];
                                    inelastic = AbstractOceanInelasticProcess{Float64}[WaterRaman{Float64}()]),
                         λ; ℓ_max = 2)
        E₀    = [1.0, 0.8, 0.5, 0.2, 1e-9, 0.0, 0.0]
        noisy = E₀ .+ [0, 0, 0, 0, 0, -1e-14, -1e-9]            # down to −1e-9 of the peak
        @test absorbed_power(o, noisy) == absorbed_power(o, E₀)
        @test absorbed_photons(o, noisy) == absorbed_photons(o, E₀)
        @test absorbed_par(o, noisy) == absorbed_par(o, E₀)
        @test inelastic_emission(o, 1, noisy, [520.0, 560.0]; range = (400, 700)) ==
              inelastic_emission(o, 1, E₀, [520.0, 560.0]; range = (400, 700))
        @test all(≥(0), absorbed_power(o, noisy))
        # √eps(Float64) ≈ 1.5e-8 of the peak is the limit.
        @test absorbed_power(o, E₀ .+ [0, 0, 0, 0, 0, 0, -1e-8]) == absorbed_power(o, E₀)
        for bad in (-2e-8, -1e-3, -1.0)
            @test_throws "below the round-off allowance" absorbed_power(o, E₀ .+ [0, 0, 0, 0, 0, 0, bad])
            @test_throws "below the round-off allowance" absorbed_par(o, E₀ .+ [0, 0, 0, 0, 0, 0, bad])
            @test_throws "below the round-off allowance" inelastic_emission(
                o, 1, E₀ .+ [0, 0, 0, 0, 0, 0, bad], [520.0]; range = (400, 700))
        end
        # The allowance scales with the field: the same −1e-9 is not round-off of a 1e-6 peak.
        @test_throws "below the round-off allowance" absorbed_power(o, 1e-6 .* E₀ .+ [0, 0, 0, 0, 0, 0, -1e-9])
        @test_throws "must be finite" absorbed_power(o, [NaN; E₀[2:end]])
        # Single precision: √eps(Float32) ≈ 3.5e-4.
        @test absorbed_power(o, Float32.(E₀) .+ Float32[0, 0, 0, 0, 0, 0, -1f-4]) == absorbed_power(o, Float32.(E₀))
        @test_throws "below the round-off allowance" absorbed_power(
            o, Float32.(E₀) .+ Float32[0, 0, 0, 0, 0, 0, -1f-3])
    end

    @testset "Bottoms and interface index" begin
        column = uniform_column(AbstractOceanConstituent{Float64}[PureWater{Float64}()]; n_layers = 2, depth = 10.0)
        col = column_optics(column, [450.0, 550.0]; refractive_index = FixedRefractiveIndex(1.34),
                            bottom = UnknownTestBottom(), ℓ_max = 2)
        @test col.bottom === UnknownTestBottom()
        @test_throws ArgumentError bottom_albedo(col)
        @test_throws "UnknownTestBottom has no Lambertian albedo" bottom_albedo(col)

        @test refractive_index(FixedRefractiveIndex(1.0), 550.0) == 1.0     # index-matched interface
        @test FixedRefractiveIndex(1) isa FixedRefractiveIndex{Float64}
        for n in (0.9, 1 - 1e-12, 0.0, -1.34, Inf, NaN)
            @test_throws ArgumentError FixedRefractiveIndex(n)
        end
        @test_throws "n_water / n_air" FixedRefractiveIndex(1 / 1.34)
        # A column built directly applies the same rule.
        @test OceanColumnOptics(col.layers, [1.0, 1.34], BlackBottom()) isa OceanColumnOptics
        @test_throws "n_water / n_air" OceanColumnOptics(col.layers, [1.34, 1 / 1.34], BlackBottom())
        @test_throws "n_water / n_air" OceanColumnOptics(col.layers, [1.34, NaN], BlackBottom())
    end

    @testset "Composite pure-water absorption names its range" begin
        w = PureWater{Float64}()                                  # IOCCG2018: 350–2500 nm
        @test absorption(w, 350.0) > 0 && absorption(w, 2500.0) > 0
        for λ in (349.9, 2500.1)
            err = try absorption(w, λ) catch e; e end
            @test err isa DomainError && err.val == λ
            @test occursin("350.0–2500.0 nm", sprint(showerror, err))
        end
    end

    @testset "zero_padded: six equal-length coefficient vectors" begin
        g = GreekCoefficients([1.0, 0.3, 0.5])
        p = zero_padded(g)
        @test p.β == [1.0, 0.3, 0.5] && is_polarized(p)
        @test p.α == p.γ == p.δ == p.ϵ == p.ζ == [0.0, 0.0, 0.0]
        @test all(length(getfield(p, f)) == 3 for f in fieldnames(GreekCoefficients))
        @test zero_padded(GreekCoefficients(Float32[1, 0.3])) isa GreekCoefficients{Float32}
        full = greek_expansion(RayleighWaterPhase{Float64}(), 4; polarized = true)
        @test zero_padded(full) === full                          # already polarized: unchanged
        @test any(!iszero, full.γ)
    end
end
