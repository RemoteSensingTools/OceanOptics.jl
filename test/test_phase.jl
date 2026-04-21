using QuadGK

@testset "Phase functions" begin

    @testset "Rayleigh (depolarized) — de Haan convention" begin
        # Pure-Rayleigh (ρ=0): β₀ = 1, β₂ = 1/2, all other moments zero.
        pf0 = RayleighWaterPhase{Float64}(depolarization = 0.0)
        β   = phase_function_moments(pf0, 8)
        @test β[1] ≈ 1.0                    # β₀ (normalization)
        @test β[2] ≈ 0.0                    # β₁ (asymmetry = 0)
        @test β[3] ≈ 0.5                    # β₂ (pure Rayleigh, de Haan)
        @test all(β[4:end] .≈ 0.0)

        # Normalization of the continuous form: (1/2) ∫₋₁¹ p(μ) dμ = 1
        integrand(μ) = phase_function_value(pf0, μ) / 2
        ∫p = 0.0
        for i in 1:200
            μ  = -1 + 2 * (i - 0.5) / 200
            ∫p += integrand(μ) * (2 / 200)
        end
        @test ∫p ≈ 1.0 rtol=1e-3

        # Forward/back symmetry: B = 1/2 exactly
        @test backscatter_fraction(pf0) == 0.5

        # With depolarization ρ, β₂ = ½·dpl_p where dpl_p = (1-ρ)/(1+ρ/2).
        ρ     = 0.09
        pf_ρ  = RayleighWaterPhase{Float64}(depolarization = ρ)
        dpl_p = (1 - ρ) / (1 + ρ/2)
        β_ρ   = phase_function_moments(pf_ρ, 2)
        @test β_ρ[3] ≈ 0.5 * dpl_p          rtol=1e-12
    end

    @testset "Rayleigh phase_matrix_moments (polarized Greek)" begin
        # Reference values from Sanghavi (2014) Eq. 16 / vSmartMOM
        # `get_greek_rayleigh`, reproduced here for the pure-Rayleigh
        # (ρ = 0) and typical-water (ρ = 0.039) cases.
        for ρ in (0.0, 0.039, 0.09)
            pf     = RayleighWaterPhase{Float64}(depolarization = ρ)
            dpl_p  = (1 - ρ) / (1 + ρ/2)
            dpl_r  = (1 - 2ρ) / (1 - ρ)
            m      = phase_matrix_moments(pf, 4)

            # ℓ = 0 row
            @test m.β[1] ≈ 1.0
            @test iszero(m.α[1]) && iszero(m.γ[1]) &&
                  iszero(m.δ[1]) && iszero(m.ϵ[1]) && iszero(m.ζ[1])

            # ℓ = 1: only δ non-zero
            @test m.δ[2] ≈ 1.5 * dpl_p * dpl_r
            @test iszero(m.α[2]) && iszero(m.β[2]) && iszero(m.γ[2])
            @test iszero(m.ϵ[2]) && iszero(m.ζ[2])

            # ℓ = 2: α, β, γ non-zero
            @test m.α[3] ≈ 3 * dpl_p
            @test m.β[3] ≈ 0.5 * dpl_p
            @test m.γ[3] ≈ dpl_p * sqrt(1.5)
            @test iszero(m.δ[3]) && iszero(m.ϵ[3]) && iszero(m.ζ[3])

            # ℓ ≥ 3: all zero
            @test all(iszero, m.α[4:end])
            @test all(iszero, m.β[4:end])
            @test all(iszero, m.γ[4:end])
            @test all(iszero, m.δ[4:end])
            @test all(iszero, m.ϵ[4:end])
            @test all(iszero, m.ζ[4:end])
        end
    end

    @testset "Scalar phase_matrix_moments fallback" begin
        # Non-polarizable phase functions should populate β from the
        # scalar moments and zero the other five Greek slots.
        hg = HenyeyGreensteinPhase{Float64}(g = 0.7)
        m  = phase_matrix_moments(hg, 6)
        @test m.β == phase_function_moments(hg, 6)
        @test all(iszero, m.α) && all(iszero, m.γ) && all(iszero, m.δ)
        @test all(iszero, m.ϵ) && all(iszero, m.ζ)
        @test !is_polarizable(hg)
    end

    @testset "Henyey-Greenstein closed form (de Haan convention)" begin
        # In the de Haan convention β_ℓ = (2ℓ+1) · g^ℓ.
        for g in (-0.5, 0.0, 0.3, 0.85)
            pf = HenyeyGreensteinPhase{Float64}(g = g)
            β  = phase_function_moments(pf, 12)
            for ℓ in 0:12
                @test β[ℓ + 1] ≈ (2ℓ + 1) * g^ℓ
            end
            @test asymmetry_parameter(pf) == g
        end

        # Normalization: (1/2) ∫ p dμ = 1 for any |g| < 1. Use QuadGK to
        # cope with strong forward peaks (g=0.85 has ~(1+g)^(-3/2) cusp
        # near μ=1). The phase function is analytic for |g|<1, so adaptive
        # integration converges cheaply.
        for g in (0.0, 0.5, -0.3, 0.85)
            pf = HenyeyGreensteinPhase{Float64}(g = g)
            ∫p, _ = QuadGK.quadgk(μ -> phase_function_value(pf, μ), -1.0, 1.0;
                                  rtol = 1e-8)
            @test ∫p / 2 ≈ 1.0 rtol=1e-6
        end

        # Backscatter fraction: closed form vs numerical integration
        for g in (0.0, 0.3, 0.85, -0.2)
            pf      = HenyeyGreensteinPhase{Float64}(g = g)
            B_cf    = backscatter_fraction(pf)
            # Numerical reference on back hemisphere
            B_ref = 0.0
            N = 500
            for i in 1:N
                μ  = -1 + (i - 0.5) / N                      # on [-1, 0)
                B_ref += phase_function_value(pf, μ) * (1 / N) / 2
            end
            @test B_cf ≈ B_ref rtol=5e-3
        end

        # Range validation
        @test_throws ArgumentError HenyeyGreensteinPhase{Float64}(g = 1.0)
        @test_throws ArgumentError HenyeyGreensteinPhase{Float64}(g = -1.0)
    end

    @testset "Fournier-Forand normalization and Mobley inversion" begin
        pf = FournierForandPhase{Float64}(backscatter_fraction = 0.018)

        # Normalization: (1/2) ∫ p dμ = 1 to the renorm factor's tolerance.
        # FF has p ~ 1/sin²(θ/2) near μ = 1 — use QuadGK with the same
        # forward-peak subdivisions the renorm factor uses.
        integrand(μ) = phase_function_value(pf, μ)
        I, _  = QuadGK.quadgk(integrand,
                              -1.0, 0.9, 0.999, 0.99999, 1.0; rtol = 1e-8)
        @test I / 2 ≈ 1.0 rtol=1e-6

        # Strongly forward-peaked: value at μ=0.99 much greater than μ=0
        @test phase_function_value(pf, 0.99) > 10 * phase_function_value(pf, 0.0)

        # Legendre moments: after post-normalization in numerical_moments,
        # β₀ = 1 exactly; β₁ > 0 (forward-peaked).
        β = phase_function_moments(pf, 64)
        @test β[1] ≈ 1.0 rtol=1e-12
        @test β[2]   > 0.5
    end

    @testset "Phase-moment mixing" begin
        pf1 = HenyeyGreensteinPhase{Float64}(g = 0.3)
        pf2 = HenyeyGreensteinPhase{Float64}(g = 0.85)
        β1  = phase_function_moments(pf1, 6)
        β2  = phase_function_moments(pf2, 6)

        # Equal weights → simple average
        β_mix = mix_phase_moments([1.0, 1.0], [β1, β2])
        @test β_mix ≈ (β1 .+ β2) ./ 2

        # Weighted mix: b₂ ≫ b₁ → β_mix close to β2
        β_mix = mix_phase_moments([0.01, 10.0], [β1, β2])
        @test β_mix[2] ≈ β2[2] rtol=5e-3

        # Degenerate: all zero scatterings → unit isotropic
        β_iso = mix_phase_moments([0.0, 0.0], [β1, β2])
        @test β_iso[1] ≈ 1.0
        @test all(β_iso[2:end] .≈ 0.0)
    end

    @testset "PetzoldPhase (bundled, CC-BY via Ocean Optics Web Book)" begin
        pf = PetzoldPhase{Float64}()
        # Normalized: (1/2) ∫₋₁¹ p dμ should equal 1 after the loader's
        # trapezoid normalization
        β = phase_function_moments(pf, 8)
        @test β[1] ≈ 1.0 rtol=5e-2

        # Strongly forward-peaked ocean particulate scattering
        @test phase_function_value(pf, 0.99) > 50 * phase_function_value(pf, 0.0)
    end
end
