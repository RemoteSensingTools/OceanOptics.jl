using QuadGK

# Generic quadrature path, bypassing closed-form overrides.
numeric_moments(pf, ℓ_max) =
    invoke(phase_function_moments, Tuple{AbstractOceanPhaseFunction, Integer}, pf, ℓ_max)
numeric_backscatter(pf) =
    invoke(backscatter_fraction, Tuple{AbstractOceanPhaseFunction}, pf)

@testset "Phase functions" begin

    @testset "Rayleigh water: closed forms (Hansen & Travis 1974)" begin
        pf0 = RayleighWaterPhase{Float64}(depolarization = 0.0)
        β   = phase_function_moments(pf0, 8)
        @test β[1:3] ≈ [1.0, 0.0, 0.5] && all(iszero, β[4:end])
        @test backscatter_fraction(pf0) == 0.5
        # The generic θ quadrature reproduces the closed form.
        @test numeric_moments(pf0, 8) ≈ β atol = 1e-12

        for ρ in (0.0, 0.039, 0.09)
            pf    = RayleighWaterPhase{Float64}(depolarization = ρ)
            dpl_p = (1 - ρ) / (1 + ρ / 2)
            dpl_r = (1 - 2ρ) / (1 - ρ)
            g     = greek_expansion(pf, 4; polarized = true)
            @test g isa GreekCoefficients && is_polarized(g)
            @test g.β ≈ [1, 0, 0.5dpl_p, 0, 0]
            @test g.α ≈ [0, 0, 3dpl_p, 0, 0]
            @test g.γ ≈ [0, 0, sqrt(1.5) * dpl_p, 0, 0]
            @test g.δ ≈ [0, 1.5dpl_p * dpl_r, 0, 0, 0]
            @test all(iszero, g.ϵ) && all(iszero, g.ζ)
            @test !is_polarized(greek_expansion(pf, 4))
        end
    end

    @testset "Henyey-Greenstein: closed forms and generic quadrature" begin
        for g in (-0.5, 0.0, 0.3, 0.85)
            pf = HenyeyGreensteinPhase{Float64}(g = g)
            β  = phase_function_moments(pf, 12)
            @test β ≈ [(2ℓ + 1) * g^ℓ for ℓ in 0:12]
            @test asymmetry_parameter(pf) == g
        end
        # The quadrature used for every phase function without a closed form
        # reproduces HG moments and backscatter, also for a strong forward peak.
        for g in (0.5, 0.95, 0.99)
            pf = HenyeyGreensteinPhase{Float64}(g = g)
            @test maximum(abs, (numeric_moments(pf, 400) .- phase_function_moments(pf, 400)) ./
                               (2 .* (0:400) .+ 1)) < 1e-10
            @test numeric_backscatter(pf) ≈ backscatter_fraction(pf) rtol = 1e-10
        end
        @test_throws ArgumentError HenyeyGreensteinPhase{Float64}(g = 1.0)
    end

    @testset "Fournier-Forand" begin
        # Mobley, Sundman & Boss (2002): n = 1.10, μ = 3.5835 gives B = 0.0183
        # and lies on their family n = 1.01 + 0.1542 (μ − 3).
        pf = FournierForandPhase{Float64}(backscatter_fraction = 0.0183)
        @test pf.n ≈ 1.10 atol = 1e-3
        @test pf.μ_j ≈ 3.5835 atol = 1e-3
        @test backscatter_fraction(FournierForandPhase{Float64}(n = 1.10, μ_j = 3.5835)) ≈ 0.0183 atol = 5e-5

        # Exact normalization, checked by an independent adaptive integration
        # in θ (the forward peak is p ~ θ^(2|ν|-2)), no renormalization involved.
        for B in (0.0183, 0.01, 0.001)
            pf  = FournierForandPhase{Float64}(backscatter_fraction = B)
            brk = [0.0; 10.0 .^ (-14:-1); π]
            # Evaluate through the angle (cos θ cannot resolve θ ≲ 1e-8 rad).
            fθ(θ) = OceanOptics._phase_at_angle(pf, θ) * sin(θ) / 2
            I = sum(quadgk(fθ, brk[i], brk[i+1]; rtol = 1e-12)[1] for i in 1:length(brk)-1)
            @test I ≈ 1 atol = 1e-9
            # β₁ = 3⟨cos Θ⟩ from the same independent integration.
            g1 = sum(quadgk(θ -> fθ(θ) * cos(θ), brk[i], brk[i+1]; rtol = 1e-12)[1]
                     for i in 1:length(brk)-1)
            @test phase_function_moments(pf, 256)[2] ≈ 3g1 rtol = 1e-9
            @test numeric_backscatter(pf) ≈ backscatter_fraction(pf) rtol = 1e-8
        end
        @test phase_function_value(pf, 0.99) > 10 * phase_function_value(pf, 0.0)
        @test isfinite(phase_function_value(pf, 1.0))
        # δ = 1 (θ ≈ 9.94° for n = 1.10) is a removable 0/0 of the closed form:
        # compare with the closed form evaluated in extended precision.
        p1 = FournierForandPhase{Float64}(n = 1.10, μ_j = 3.5835)
        function ff_big(n, μj, μ)
            n, μj, s = big(n), big(μj), (1 - big(μ)) / 2
            ν = (3 - μj) / 2; δ = 4s / (3(n - 1)^2); δ180 = 4 / (3(n - 1)^2)
            main = (ν * (1 - δ) - (1 - δ^ν) + (δ * (1 - δ^ν) - ν * (1 - δ)) / s) / (4π * (1 - δ)^2 * δ^ν)
            back = (1 - δ180^ν) / (16π * (δ180 - 1) * δ180^ν) * (3 * (1 - 2s)^2 - 1)
            return Float64(4π * (main + back))
        end
        setprecision(BigFloat, 256) do
            for μ in (0.985, 0.98497, 0.98498, 0.9849795, 0.99, 0.98)
                @test phase_function_value(p1, μ) ≈ ff_big(1.10, 3.5835, μ) rtol = 1e-12
            end
        end
        μc = 1 - 2 * 3 * 0.1^2 / 4                        # δ = 1 exactly
        @test isfinite(phase_function_value(p1, μc))
        @test phase_function_value(p1, μc) ≈ phase_function_value(p1, μc + 1e-9) rtol = 1e-6
        # Very small B: most scattering lies below the quadrature floor and is
        # carried by the analytic forward tail; normalization stays exact.
        pt = FournierForandPhase{Float64}(backscatter_fraction = 1e-5)
        βt = phase_function_moments(pt, 64)
        @test βt[1] == 1 && 0.99 < βt[2] / 3 < 1
        @test_throws ArgumentError FournierForandPhase{Float64}(n = 1.9, μ_j = 3.5)
        @test_throws ArgumentError FournierForandPhase{Float64}(backscatter_fraction = 0.6)
        @test_throws ArgumentError FournierForandPhase{Float64}(n = 0.9, μ_j = 3.5)
        @test_throws ArgumentError FournierForandPhase{Float64}(n = 1.1, μ_j = 5.5)
        @test_throws ArgumentError FournierForandPhase{Float64}()
    end

    @testset "Petzold (San Diego Harbor)" begin
        pf = PetzoldPhase{Float64}()
        # Independent check of the tabulated VSF and its extrapolation:
        # Mobley (1994) Table 3.10 gives b = 1.824 m⁻¹ for San Diego Harbor.
        @test pf.b ≈ 1.824 rtol = 1e-2
        @test phase_function_moments(PetzoldPhase{Float32}(), 16)[1] ≈ 1
        @test 0.015 < backscatter_fraction(pf) < 0.025
        @test 0.90 < asymmetry_parameter(pf) < 0.94
        @test phase_function_value(pf, cosd(0.1)) ≈ pf.p[1]
        @test phase_function_value(pf, 0.99) > 50 * phase_function_value(pf, 0.0)
    end

    @testset "Scalar-only phase functions have no phase matrix" begin
        for pf in (HenyeyGreensteinPhase{Float64}(g = 0.7),
                   FournierForandPhase{Float64}(backscatter_fraction = 0.01),
                   PetzoldPhase{Float64}())
            @test !is_polarizable(pf)
            @test_throws ArgumentError greek_expansion(pf, 6; polarized = true)
            g = depolarizing_greek_expansion(pf, 6)
            @test is_polarized(g) && g.β == phase_function_moments(pf, 6)
            @test all(iszero, (g.α, g.γ, g.δ, g.ϵ, g.ζ))
        end
    end

    @testset "Wavelength argument is ignored by wavelength-independent phases" begin
        hg = HenyeyGreensteinPhase{Float64}(g = 0.7)
        rw = RayleighWaterPhase{Float64}()
        @test !is_wavelength_dependent(hg)
        @test phase_function_moments(hg, 6, 550.0) == phase_function_moments(hg, 6)
        @test backscatter_fraction(hg, 550.0) == backscatter_fraction(hg)
        @test greek_expansion(rw, 4, 550.0; polarized = true) ==
              greek_expansion(rw, 4; polarized = true)
    end
end
