@testset "Inference ($FT)" for FT in (Float64, Float32)
    λ = FT(550)
    for c in iop_constituents(FT)
        @test (@inferred absorption(c, λ)) isa AbstractFloat
        @test (@inferred scattering(c, λ)) isa AbstractFloat
        @test (@inferred backscattering(c, λ)) isa AbstractFloat
        @test (@inferred absorption_by_part(c, λ)) isa Tuple
    end
    bidigare = Bidigare1990Pigments(chl_a = FT(0.5), chl_b = FT(0.05), psc = FT(0.3))
    @test (@inferred pigment_fractions(bidigare, λ)) isa NTuple{8, FT}
    for pf in scalar_phases(FT)
        @test (@inferred greek_expansion(pf, 16)) isa GreekCoefficients{FT}
        @test (@inferred backscatter_fraction(pf)) isa FT
    end
    @test (@inferred greek_expansion(RayleighWaterPhase{FT}(), 16; polarized = true)) isa GreekCoefficients{FT}
    @test (@inferred refractive_index(FixedRefractiveIndex(FT(1.34)), λ)) isa FT
    @test (@inferred phase_id(greek_expansion(RayleighWaterPhase{FT}(), 4))) isa UInt64

    @testset "layer and column optics" begin
        layer = inelastic_layer(FT)
        λg    = collect(FT, 350:10:750)
        o     = @inferred layer_optics(FT, layer, λg; ℓ_max = 8)
        @test typeof(o) === OceanLayerOptics{FT, Vector{FT}, Vector{Int}}
        @test (@inferred layer_optics(layer, λg; ℓ_max = 8)) isa OceanLayerOptics{FT}
        @test (@inferred layer_optics(FT, layer, λg; ℓ_max = 8, polarized = true,
                                      scalar_phase_approximation = :depolarizing)) isa OceanLayerOptics{FT}
        column = uniform_column(layer.constituents; n_layers = 2, depth = FT(10))
        index  = FixedRefractiveIndex(FT(1.34))
        @test (@inferred column_optics(FT, column, λg; refractive_index = index,
                                       bottom = BlackBottom(), ℓ_max = 8)) isa
              OceanColumnOptics{FT, BlackBottom, Vector{FT}, Vector{Int}}
        @test (@inferred column_optics(FT, column, λg; refractive_index = index,
                                       bottom = LambertianBottom(FT(0.1)), ℓ_max = 8)) isa
              OceanColumnOptics{FT, LambertianBottom{FT, Vector{FT}}, Vector{FT}, Vector{Int}}

        E₀ = ones(FT, length(λg))
        @test (@inferred absorbed_power(o, E₀)) isa Matrix{FT}
        @test (@inferred absorbed_photons(o, E₀)) isa Matrix
        @test (@inferred absorbed_par(o, E₀)) isa Vector
        λ_em = collect(FT, 660:10:740)
        for j in eachindex(o.inelastic)          # fluorescence, water Raman
            @test (@inferred excitation_optical_depth(o, j)) isa Vector{FT}
            @test (@inferred inelastic_emission(o, j, E₀, λ_em)) isa Vector{FT}
        end
        @test all(c -> isconcretetype(typeof(c)), o.inelastic)
    end
end
