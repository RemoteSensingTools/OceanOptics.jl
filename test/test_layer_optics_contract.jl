@testset "OceanLayerOptics consumer contract" begin

    # Minimal valid scalar bundle: two wavelengths, Rayleigh-like β.
    valid(; kw...) = OceanLayerOptics(; λ = [450.0, 550.0], Δz = 1.0,
                                      τ = [0.1, 0.2], ϖ = [0.5, 0.4], B = [0.5, 0.5],
                                      β = [1.0 1.0; 0.0 0.0; 0.5 0.5], kw...)
    pol(; kw...) = (α = zeros(3, 2), γ = zeros(3, 2), δ = zeros(3, 2),
                    ϵ = zeros(3, 2), ζ = zeros(3, 2), kw...)

    @testset "Valid scalar and polarized bundles" begin
        o = valid()
        @test !is_polarized(o)
        @test n_moments(o) == 3 && n_wavelengths(o) == 2
        @test is_polarized(valid(; pol()...))
        # Decreasing wavelength order (e.g. wavenumber-ordered consumers) is fine.
        @test n_wavelengths(valid(; λ = [550.0, 450.0])) == 2
        # Zero optical depth and a non-scattering wavelength are physical.
        @test valid(; τ = [0.0, 0.2], ϖ = [0.0, 0.4]) isa OceanLayerOptics
    end

    @testset "Partial polarized bundles are rejected" begin
        # Previously `is_polarized` only looked at α, so these looked polarized.
        @test_throws ArgumentError valid(; α = zeros(3, 2))
        @test_throws ArgumentError valid(; pol(ζ = zeros(0, 0))...)
        @test_throws DimensionMismatch valid(; pol(γ = zeros(2, 2))...)
        @test_throws ArgumentError valid(; pol(δ = fill(NaN, 3, 2))...)
    end

    @testset "Wavelength grid" begin
        @test_throws ArgumentError valid(; λ = [450.0, 450.0])
        @test_throws ArgumentError valid(λ = [450.0, 550.0, 500.0], τ = [0.1, 0.2, 0.1],
                                         ϖ = [0.5, 0.4, 0.4], B = [0.5, 0.5, 0.5],
                                         β = [1.0 1.0 1.0])
        @test_throws ArgumentError valid(; λ = [-450.0, 550.0])
        @test_throws ArgumentError valid(; λ = [450.0, Inf])
        @test_throws DimensionMismatch valid(; τ = [0.1, 0.2, 0.3])
        @test_throws DimensionMismatch valid(; β = [1.0 1.0 1.0])
    end

    @testset "Coefficient ranges" begin
        @test_throws ArgumentError valid(; Δz = 0.0)
        @test_throws ArgumentError valid(; τ = [-0.1, 0.2])
        @test_throws ArgumentError valid(; τ = [NaN, 0.2])
        @test_throws ArgumentError valid(; ϖ = [1.1, 0.4])
        @test_throws ArgumentError valid(; B = [0.5, -0.01])
        @test_throws ArgumentError valid(; J_source = [1.0, -1.0])
        @test_throws DimensionMismatch valid(; J_source = [1.0])
        @test valid(; J_source = [1.0, 0.0]) isa OceanLayerOptics
    end

    @testset "Phase normalization and moment bounds" begin
        @test_throws ArgumentError valid(; β = [0.98 1.0; 0.0 0.0; 0.5 0.5])
        # Henyey-Greenstein g = 0.8 with the (2ℓ+1) factor applied twice:
        # β₁ = 3·(3g) = 7.2 > 3 is impossible for any non-negative phase function.
        g = 0.8
        β_twice = [(2ℓ + 1)^2 * g^ℓ for ℓ in 0:2]
        @test_throws ArgumentError valid(; β = [β_twice β_twice])
        β_once = [(2ℓ + 1) * g^ℓ for ℓ in 0:2]
        @test valid(; β = [β_once β_once]) isa OceanLayerOptics
    end

    @testset "Low-order polarized coefficients" begin
        # α, γ, ϵ, ζ multiply GSFs that vanish for ℓ < 2: a Rayleigh table
        # shifted up by one row must be rejected.
        shifted = pol(α = [0.0 0.0; 3.0 3.0; 0.0 0.0], γ = [0.0 0.0; sqrt(1.5) sqrt(1.5); 0.0 0.0])
        @test_throws ArgumentError valid(; shifted...)
        # δ uses P_ℓ like β, so δ₁ ≠ 0 is legitimate (Rayleigh δ₁ = 3/2).
        @test is_polarized(valid(; pol(δ = [0.0 0.0; 1.5 1.5; 0.0 0.0])...))
    end

    @testset "Array types and constructors" begin
        @test_throws ArgumentError valid(; τ = Float32[0.1, 0.2])
        @test_throws ArgumentError valid(; α = zeros(Float32, 3, 2))
        # The implicit positional constructor of earlier releases still works.
        o = valid()
        @test OceanLayerOptics(o.λ, o.Δz, o.τ, o.ϖ, o.B, o.β, o.α, o.γ, o.δ, o.ϵ, o.ζ,
                               nothing) isa OceanLayerOptics
        # The element type is not restricted: a differentiable bundle constructs.
        D  = ForwardDiff.Dual{Nothing}(0.1, 1.0)
        z  = zero(D)
        od = OceanLayerOptics(; λ = [450.0 + z, 550.0 + z], Δz = 1.0, τ = [D, 2D],
                              ϖ = [0.5 + z, 0.4 + z], B = [0.5 + z, 0.5 + z],
                              β = [1.0+z 1.0+z; z z; 0.5+z 0.5+z])
        @test od.τ[2] == 2D
    end

    @testset "Every shipped constituent passes the contract" begin
        layer = OceanLayer(0.0, 2.0, AbstractOceanConstituent{Float64}[
            PureWater{Float64}(), CDOM{Float64}(a_ref = 0.05),
            Phytoplankton{Float64, Bricaud1998}(Chl = 2.0), NonAlgalParticles{Float64}(),
            NonAlgalParticles{Float64, FellMineral}()])
        λ = collect(350.0:25.0:750.0)
        for FT in (Float64, Float32), ℓ_max in (8, 128)
            @test layer_optics(layer, λ; ℓ_max, FT) isa OceanLayerOptics{FT}
            @test is_polarized(layer_optics(layer, λ; ℓ_max, FT, polarized = true,
                                            scalar_phase_approximation = :depolarizing))
        end
        @test_throws ArgumentError layer_optics(layer, [550.0, 450.0, 500.0]; ℓ_max = 8)
        @test layer_optics(layer, reverse(λ); ℓ_max = 8).τ ≈ reverse(layer_optics(layer, λ; ℓ_max = 8).τ)
    end

    @testset "Polarized zero-scattering wavelength" begin
        λ     = [450.0, 550.0]
        s     = TabulatedScatterer(λ, [0.0, 0.1], RayleighWaterPhase{Float64}())
        a     = TabulatedAbsorber(λ, [0.2, 0.2])
        opt   = layer_optics(OceanLayer(0.0, 1.0, AbstractOceanConstituent{Float64}[a, s]),
                             λ; ℓ_max = 4, polarized = true)
        @test opt.ϖ[1] == 0 && opt.B[1] == 0
        @test opt.β[:, 1] == [1.0, 0, 0, 0, 0]
        @test all(iszero, (opt.α[:, 1], opt.γ[:, 1], opt.δ[:, 1], opt.ϵ[:, 1], opt.ζ[:, 1]))
        @test opt.α[3, 2] > 0
    end

    @testset "Polarized layer_optics requires an explicit particle choice" begin
        water = PureWater{Float64}(depolarization = 0.0)
        phyto = Phytoplankton{Float64, Bricaud1998}(Chl = 1.0)
        layer = OceanLayer(0.0, 1.0, AbstractOceanConstituent{Float64}[water, phyto])
        λ     = [450.0, 550.0]

        # Fournier-Forand particles have no phase matrix: no silent fallback.
        @test_throws ArgumentError layer_optics(layer, λ; ℓ_max = 8, polarized = true)
        @test_throws ArgumentError layer_optics(layer, λ; ℓ_max = 8, polarized = true,
                                                scalar_phase_approximation = :rayleigh)

        opt = layer_optics(layer, λ; ℓ_max = 8, polarized = true,
                           scalar_phase_approximation = :depolarizing)
        @test is_polarized(opt)
        # The depolarizing particles contribute only to β, so the polarized
        # coefficients are the ρ = 0 Rayleigh values (Hansen & Travis 1974:
        # α₂ = 3, γ₂ = √(3/2), δ₁ = 3/2) weighted by the water scattering share.
        for (i, λi) in enumerate(λ)
            w = scattering(water, λi) / (scattering(water, λi) + scattering(phyto, λi))
            @test opt.α[3, i] ≈ 3w             rtol = 1e-12
            @test opt.γ[3, i] ≈ sqrt(1.5) * w  rtol = 1e-12
            @test opt.δ[2, i] ≈ 1.5w           rtol = 1e-12
            @test all(iszero, opt.ϵ[:, i]) && all(iszero, opt.ζ[:, i])
        end
        # The scalar result is unchanged by the polarization choice.
        @test opt.β ≈ layer_optics(layer, λ; ℓ_max = 8).β
        # The keyword has no effect on scalar output.
        @test layer_optics(layer, λ; ℓ_max = 8, scalar_phase_approximation = :depolarizing).β ==
              layer_optics(layer, λ; ℓ_max = 8).β
    end
end
