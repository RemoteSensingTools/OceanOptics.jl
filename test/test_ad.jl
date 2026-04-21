using ForwardDiff

@testset "Forward-mode AD transparency (DESIGN §12.5)" begin

    @testset "CDOM ∂a/∂λ" begin
        cdom = CDOM{Float64}(a_ref = 0.03, λ_ref = 440.0, slope = 0.014)
        # Analytic derivative: d/dλ [a_ref · exp(-S(λ-λ_ref))] = -S · a(λ)
        for λ in (420.0, 440.0, 500.0)
            da_dλ_ad  = ForwardDiff.derivative(x -> absorption(cdom, x), λ)
            da_dλ_ref = -cdom.slope * absorption(cdom, λ)
            @test da_dλ_ad ≈ da_dλ_ref rtol=1e-10
        end
    end

    @testset "PureWater absorption through interpolation + Pegau correction" begin
        w = PureWater{Float64}(absorption_model = PopeFry1997(),
                               temperature = 20.0)
        # Smooth (piecewise-linear interp) — AD should return finite derivative
        for λ in (450.0, 500.0, 650.0)
            d = ForwardDiff.derivative(x -> absorption(w, x), λ)
            @test isfinite(d)
        end
    end

    @testset "Morel 1974 scattering ∂b/∂λ" begin
        w = PureWater{Float64}(scattering_model = Morel1974())
        # b ∝ λ^{-4.32} → ∂b/∂λ = -4.32 · b / λ
        for λ in (400.0, 500.0, 700.0)
            db_dλ_ad  = ForwardDiff.derivative(x -> scattering(w, x), λ)
            db_dλ_ref = -4.32 * scattering(w, λ) / λ
            @test db_dλ_ad ≈ db_dλ_ref rtol=1e-10
        end
    end

    @testset "Fournier-Forand phase ∂p/∂cosθ" begin
        pf = FournierForandPhase{Float64}(backscatter_fraction = 0.018)
        for μ in (-0.5, 0.0, 0.5, 0.9)
            d = ForwardDiff.derivative(x -> phase_function_value(pf, x), μ)
            @test isfinite(d)
        end
    end

    @testset "Henyey-Greenstein ∂p/∂cosθ matches closed form" begin
        g = 0.7
        pf = HenyeyGreensteinPhase{Float64}(g = g)
        # p = (1-g²)(1+g²-2gμ)^{-3/2} ⇒ ∂p/∂μ = 3g(1-g²)(1+g²-2gμ)^{-5/2}
        for μ in (-0.5, 0.0, 0.5)
            dp_ad  = ForwardDiff.derivative(x -> phase_function_value(pf, x), μ)
            dp_ref = 3 * g * (1 - g^2) / (1 + g^2 - 2g*μ)^2.5
            @test dp_ad ≈ dp_ref rtol=1e-10
        end
    end

    @testset "CDOM ∂a/∂slope (parametric gradient)" begin
        # Gradient of the summed CDOM absorption w.r.t. the spectral slope.
        # `CDOM` itself requires `FT<:AbstractFloat`, so we differentiate
        # through the analytic formula directly rather than instantiating
        # the struct with a `Dual` type parameter. Integration window kept
        # above `λ_ref` so increasing `S` unambiguously shrinks the sum
        # (a mixed integration window straddling `λ_ref` can flip the sign).
        a_ref, λ_ref = 0.03, 440.0
        f(S) = sum(a_ref * exp(-S * (λ - λ_ref))
                   for λ in (500.0, 600.0, 700.0))
        d = ForwardDiff.derivative(f, 0.014)
        @test isfinite(d) && d < 0
    end
end
