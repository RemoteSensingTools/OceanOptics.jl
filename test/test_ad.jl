using ForwardDiff

@testset "Forward-mode automatic differentiation" begin

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

    @testset "CDOM ∂a/∂slope and ∂a/∂a_ref match the closed forms" begin
        # a = a_ref exp(−S (λ − λ_ref)) ⇒ ∂a/∂S = −(λ − λ_ref) a and ∂a/∂a_ref = a / a_ref.
        a_ref, λ_ref, S = 0.03, 440.0, 0.014
        for λ in (412.5, 500.0, 600.0, 700.0)
            a = absorption(CDOM{Float64}(; a_ref, λ_ref, slope = S), λ)
            @test ForwardDiff.derivative(x -> absorption(CDOM(; a_ref, λ_ref, slope = x), λ), S) ≈
                  -(λ - λ_ref) * a rtol = 1e-12
            @test ForwardDiff.derivative(x -> absorption(CDOM(; a_ref = x, λ_ref, slope = S), λ), a_ref) ≈
                  a / a_ref rtol = 1e-12
        end
    end

    @testset "Absorbed light ∂/∂E₀" begin
        λ = collect(400.0:50.0:700.0)
        o = layer_optics(OceanLayer(0.0, 3.0, AbstractOceanConstituent{Float64}[PureWater{Float64}(),
                                                                               CDOM{Float64}()]),
                         λ; ℓ_max = 2)
        E₀ = [1.0, 0.9, 0.8, 0.6, 0.5, 0.3, 0.2]
        τ  = absorption_optical_depth(o)
        # A[c, k] = τ_c[k] E₀[k]: the Jacobian of Σ_c A[c, ·] is diag(τ).
        J = ForwardDiff.jacobian(E -> vec(sum(absorbed_power(o, E); dims = 1)), E₀)
        @test J ≈ [i == j ? τ[i] : 0.0 for i in eachindex(τ), j in eachindex(τ)] rtol = 1e-13
        # Photon units and absorbed PAR are linear in the irradiance, so the
        # derivative along a scaling x·E₀ is the value at E₀.
        @test ForwardDiff.derivative(x -> sum(absorbed_photons(o, x .* E₀)), 1.3) ≈
              sum(absorbed_photons(o, E₀)) rtol = 1e-13
        @test ForwardDiff.derivative(x -> absorbed_par(o, x .* E₀), 0.7) ≈ absorbed_par(o, E₀) rtol = 1e-13
    end

    # -------------------------------------------------------------------------
    # Derivatives with respect to model parameters
    # -------------------------------------------------------------------------
    # A model built with a dual parameter has the dual element type. The
    # reference is a central difference of the same call with plain `Float64`
    # parameters, step h = 1e-5·|x|: truncation ~h²f‴/6 and round-off ~ε/h give
    # 1e-8 or better relative to the derivative for every case below (the
    # largest observed is 6e-9). Stated tolerance: rtol = 1e-6, and atol = 1e-12
    # where the derivative is exactly zero.
    central(f, x; h = 1e-5 * abs(x)) = (f(x + h) .- f(x - h)) ./ (2h)
    agree(ad, fd; rtol = 1e-6, atol = 1e-12) = all(isapprox.(ad, fd; rtol, atol))
    check(f, x; atol = 1e-12) = agree(ForwardDiff.derivative(f, x), central(f, x); atol)
    λs  = (412.5, 443.7, 551.3, 672.1)
    one_dual(x) = ForwardDiff.Dual(x, one(x))
    D   = typeof(one_dual(0.5))

    @testset "Constructors take the element type of their parameters" begin
        @test Phytoplankton(Chl = 1) isa Phytoplankton{Float64}
        @test CDOM(a_ref = 0.5f0) isa CDOM{Float64}
        @test (@inferred Phytoplankton(Chl = one_dual(0.5))) isa Phytoplankton{D}
        @test (@inferred CDOM(slope = one_dual(0.014))) isa CDOM{D}
        @test (@inferred PureWater(temperature = one_dual(20.0))) isa PureWater{D}
        @test (@inferred NonAlgalParticles(a_ref = one_dual(0.05))) isa NonAlgalParticles{D}
        @test (@inferred FournierForandPhase(backscatter_fraction = one_dual(0.018))) isa
              FournierForandPhase{D}
        @test (@inferred CDOM{D}(a_ref = one_dual(0.03))) isa CDOM{D}
    end

    @testset "Phytoplankton ∂/∂Chl" begin
        for λ in λs
            @test check(x -> absorption(Phytoplankton(Chl = x), λ), 0.5)
            @test check(x -> scattering(Phytoplankton(Chl = x), λ), 0.5)
            @test check(x -> backscattering(Phytoplankton(Chl = x), λ), 0.5)
            # b_p = 0.30 Chl^0.62 (550/λ) ⇒ ∂b_p/∂Chl = 0.62 b_p / Chl
            @test ForwardDiff.derivative(x -> scattering(Phytoplankton(Chl = x), λ), 0.5) ≈
                  0.62 * 0.30 * 0.5^(0.62 - 1) * 550 / λ rtol = 1e-12
            # Bricaud a_p = A Chl^E is not linear in Chl: E = Chl ∂a/∂Chl / a lies in (0, 1).
            E = 0.5 * ForwardDiff.derivative(x -> absorption(Phytoplankton(Chl = x), λ), 0.5) /
                absorption(Phytoplankton(Chl = 0.5), λ)
            @test 0.5 < E < 1
        end
        pig = ProspectPigments(carotenoid_ratio = 0.3)
        @test check(x -> collect(absorption_by_part(Phytoplankton(Chl = x, pigments = pig), 443.7)), 0.5)
        @test check(x -> backscattering(Phytoplankton(Chl = 0.5, B = x), 551.3), 0.01)
        @test check(x -> scattering(Phytoplankton{typeof(x)}(
                             Chl = 0.5, scattering_model = GordonMorel1983(b_ref = x)), 551.3), 0.30)
        p = Phytoplankton(Chl = one_dual(0.5))
        @test (@inferred absorption(p, 443.7)) isa D
        @test (@inferred scattering(p, 443.7)) isa D
        @test (@inferred backscattering(p, 443.7)) isa D
    end

    @testset "Non-algal particles ∂/∂(a_ref, slope, b_ref, γ, B, c_ref, ω)" begin
        for λ in λs
            @test check(x -> absorption(NonAlgalParticles(a_ref = x), λ), 0.05)
            @test check(x -> absorption(NonAlgalParticles(slope = x), λ), 0.011)
            @test check(x -> scattering(NonAlgalParticles(b_ref = x), λ), 0.45)
            @test check(x -> scattering(NonAlgalParticles(spectral_exponent = x), λ), 1.0)
            @test check(x -> backscattering(NonAlgalParticles(B = x), λ), 0.02)
            @test check(x -> absorption(NonAlgalParticles{typeof(x), FellMineral}(c_ref = x), λ), 0.40)
            @test check(x -> scattering(NonAlgalParticles{typeof(x), FellMineral}(ω_nap = x), λ), 0.90)
            # Fell: a = (1 − ω) c_ref (440/λ) ⇒ ∂a/∂ω = −c_ref (440/λ)
            @test ForwardDiff.derivative(
                      x -> absorption(NonAlgalParticles{typeof(x), FellMineral}(ω_nap = x), λ), 0.90) ≈
                  -0.40 * 440 / λ rtol = 1e-12
        end
        @test (@inferred absorption(NonAlgalParticles(a_ref = one_dual(0.05)), 443.7)) isa D
    end

    @testset "Pure water ∂/∂(T, S, δ) and refractive index ∂n/∂(T, S)" begin
        # Absorption: a(T) = a(T_ref) + Ψ_T (T − T_ref) above 550 nm, no dependence below.
        for λ in (650.0, 720.0, 1250.0)
            @test check(x -> absorption(PureWater(temperature = x), λ), 20.0)
            @test ForwardDiff.derivative(x -> absorption(PureWater(temperature = x), λ), 20.0) ≈
                  absorption(PureWater{Float64}(temperature = 23.0), λ) -
                  absorption(PureWater{Float64}(temperature = 22.0), λ) rtol = 1e-9
        end
        @test ForwardDiff.derivative(x -> absorption(PureWater(temperature = x), 443.7), 20.0) == 0
        @test ForwardDiff.derivative(x -> absorption(PureWater(salinity = x), 650.0), 35.0) == 0
        for λ in λs
            @test check(x -> scattering(PureWater(temperature = x), λ), 20.0)
            @test check(x -> scattering(PureWater(salinity = x), λ), 35.0)
            @test check(x -> scattering(PureWater(depolarization = x), λ), 0.039)
            @test check(x -> refractive_index(QuanFry1995(temperature = x, salinity = 35.0), λ), 20.0)
            # Quan & Fry (1995): ∂n/∂S = n₁ + n₂T + n₃T² + n₆/λ
            @test ForwardDiff.derivative(
                      x -> refractive_index(QuanFry1995(temperature = 20.0, salinity = x), λ), 30.0) ≈
                  1.779e-4 - 1.05e-6 * 20 + 1.6e-8 * 20^2 + 0.01155 / λ rtol = 1e-12
        end
        # Salinity raises the scattering of seawater, warming lowers it at 20 °C.
        @test ForwardDiff.derivative(x -> scattering(PureWater(salinity = x), 500.0), 35.0) > 0
        w = PureWater(temperature = one_dual(20.0))
        @test (@inferred absorption(w, 650.0)) isa D
        @test (@inferred scattering(w, 500.0)) isa D
        @test (@inferred refractive_index(QuanFry1995(temperature = one_dual(20.0), salinity = 35.0),
                                          500.0)) isa D
    end

    @testset "Phase-function parameters" begin
        ff(x) = FournierForandPhase(backscatter_fraction = x)
        for B in (0.005, 0.018, 0.05)
            # The constructor inverts B ↦ (n, μ_j); the closed-form B of the result returns it.
            @test ForwardDiff.derivative(x -> backscatter_fraction(ff(x)), B) ≈ 1 rtol = 1e-10
            for μ in (-0.9, -0.5, 0.3, 0.9, 0.999)
                @test check(x -> phase_function_value(ff(x), μ), B)
            end
            # Carried through the Legendre quadrature; β₀ = 1 for every B.
            dβ = ForwardDiff.derivative(x -> phase_function_moments(ff(x), 16), B)
            @test agree(dβ, central(x -> phase_function_moments(ff(x), 16), B); atol = 1e-7)
            @test abs(dβ[1]) < 1e-7
            @test all(<(0), dβ[2:end])          # more backscatter, less forward peak
        end
        # The value of a dual-parameter model is that of the plain model.
        @test ForwardDiff.value(ff(one_dual(0.018)).μ_j) ≈ ff(0.018).μ_j rtol = 1e-13
        @test ForwardDiff.value.(phase_function_moments(ff(one_dual(0.018)), 16)) ≈
              phase_function_moments(ff(0.018), 16) rtol = 1e-10
        @test (@inferred phase_function_moments(ff(one_dual(0.018)), 8)) isa Vector{D}
        # Henyey–Greenstein: β_ℓ = (2ℓ+1) g^ℓ ⇒ ∂β_ℓ/∂g = (2ℓ+1) ℓ g^(ℓ−1)
        g = 0.8
        @test ForwardDiff.derivative(x -> phase_function_moments(HenyeyGreensteinPhase(g = x), 6), g) ≈
              [(2ℓ + 1) * ℓ * g^(ℓ - 1) for ℓ in 0:6] rtol = 1e-12
        @test check(x -> backscatter_fraction(HenyeyGreensteinPhase(g = x)), g)
        # Rayleigh water: β₂ = (1 − ρ)/(2 + ρ) ⇒ ∂β₂/∂ρ = −3/(2 + ρ)²
        ρ = 0.039
        @test ForwardDiff.derivative(
                  x -> phase_function_moments(RayleighWaterPhase(depolarization = x), 2)[3], ρ) ≈
              -3 / (2 + ρ)^2 rtol = 1e-12
    end

    # A mixed layer: water, CDOM, phytoplankton and non-algal particles, all at
    # the element type that the parameters promote to.
    function mixed_layer(; Chl = 0.5, a_cdom = 0.03, S_cdom = 0.014, a_nap = 0.05, T = 20.0,
                         S = 35.0, B = 0.01, Δz = 5.0)
        FT = promote_type(map(typeof, (Chl, a_cdom, S_cdom, a_nap, T, S, B, Δz))...)
        return OceanLayer(0.0, Δz, AbstractOceanConstituent{FT}[
            PureWater{FT}(temperature = T, salinity = S),
            CDOM{FT}(a_ref = a_cdom, slope = S_cdom),
            Phytoplankton{FT}(Chl = Chl, B = B),
            NonAlgalParticles{FT}(a_ref = a_nap)])
    end
    λ_layer = [412.5, 443.7, 551.3, 672.1, 720.0]
    optics(; kw...) = layer_optics(mixed_layer(; kw...), λ_layer; ℓ_max = 16)
    defaults = (Chl = 0.5, a_cdom = 0.03, S_cdom = 0.014, a_nap = 0.05, T = 20.0, S = 35.0,
                B = 0.01, Δz = 5.0)

    @testset "layer_optics: ∂(τ_a, τ, ω̄)/∂$name" for name in keys(defaults)
        x₀ = defaults[name]
        at(x) = optics(; name => x)
        @test check(x -> absorption_optical_depth(at(x)), x₀)
        @test check(x -> optical_depth(at(x)), x₀)
        # ω̄ and B are ratios of order 0.1–1 that some parameters leave unchanged (Δz); the
        # central difference of a constant is round-off, ε/h ≈ 1e-11, hence atol = 1e-10.
        @test check(x -> single_scattering_albedo(at(x)), x₀; atol = 1e-10)
        @test check(x -> backscatter_fraction(at(x)), x₀; atol = 1e-10)
        # The values are those of the plain layer.
        o_dual, o = at(one_dual(x₀)), at(x₀)
        @test ForwardDiff.value.(optical_depth(o_dual)) ≈ optical_depth(o) rtol = 1e-13
        @test ForwardDiff.value.(single_scattering_albedo(o_dual)) ≈ single_scattering_albedo(o) rtol = 1e-13
    end

    @testset "layer_optics: which optical depths respond to which parameter" begin
        dτa(name) = ForwardDiff.derivative(x -> absorption_optical_depth(optics(; name => x)), defaults[name])
        dτs(name) = ForwardDiff.derivative(x -> scattering_optical_depth(optics(; name => x)), defaults[name])
        @test all(>(0), dτa(:Chl)[1:4]) && all(>(0), dτs(:Chl))
        @test all(>(0), dτa(:a_cdom)) && all(==(0), dτs(:a_cdom))
        @test dτa(:a_cdom) ≈ [5.0 * exp(-0.014 * (λ - 440)) for λ in λ_layer] rtol = 1e-12
        @test all(==(0), dτa(:S)) && all(>(0), dτs(:S))
        @test all(==(0), dτa(:B)) && all(==(0), dτs(:B))
        @test dτa(:T)[1:2] == [0, 0] && dτa(:T)[5] > 0      # Röttgers correction from 550 nm
        # Component optical depths: only the phytoplankton component depends on Chl.
        o = optics(Chl = one_dual(0.5))
        for c in o.absorbers
            sensitive = any(!iszero, ForwardDiff.partials.(c.τ, 1))
            @test sensitive == (c.label == "Phytoplankton")
        end
    end

    @testset "layer_optics: Greek coefficients ∂β_ℓ/∂B (Fournier-Forand)" begin
        β(x, iλ) = greek_expansion(optics(B = x), iλ).β
        for iλ in (1, 3)
            dβ = ForwardDiff.derivative(x -> β(x, iλ), 0.01)
            @test agree(dβ, central(x -> β(x, iλ), 0.01); atol = 1e-7)
            @test maximum(abs, dβ) > 10          # a silent zero derivative fails here
        end
        # Phase parameters that are plain: the expansion has zero partials.
        o = optics(Chl = one_dual(0.5))
        @test all(c -> all(g -> all(iszero, ForwardDiff.partials.(g.β, 1)), c.greek), o.scatterers)
        # Polarized: the Rayleigh matrix of water follows its depolarization ratio,
        # α₂ = 3 (1 − ρ)/(1 + ρ/2) ⇒ ∂α₂/∂ρ = −4.5/(1 + ρ/2)².
        dα = ForwardDiff.derivative(0.039) do ρ
            FT = typeof(ρ)
            layer = OceanLayer(0.0, 1.0, AbstractOceanConstituent{FT}[PureWater{FT}(depolarization = ρ)])
            greek_expansion(layer_optics(layer, λ_layer; ℓ_max = 4, polarized = true), 1).α[3]
        end
        @test dα ≈ -4.5 / (1 + 0.039 / 2)^2 rtol = 1e-12
    end

    @testset "layer_optics: type stability with dual parameters" begin
        layer = mixed_layer(Chl = one_dual(0.5))
        @test layer isa OceanLayer{D}
        o = @inferred layer_optics(layer, λ_layer; ℓ_max = 8)
        @test o isa OceanLayerOptics{D}
        @test (@inferred absorption_optical_depth(o)) isa Vector{D}
        @test (@inferred optical_depth(o)) isa Vector{D}
        @test (@inferred single_scattering_albedo(o)) isa Vector{D}
        @test (@inferred greek_expansion(o, 2)) isa GreekCoefficients{D}
        # A plain layer evaluated at a dual element type: wavelength derivatives.
        plain = mixed_layer()
        @test (@inferred layer_optics(D, plain, λ_layer; ℓ_max = 8)) isa OceanLayerOptics{D}
        τ₁(x) = absorption_optical_depth(layer_optics(typeof(x), plain, [x, 551.3]; ℓ_max = 2))[1]
        @test check(τ₁, 437.3)
    end

    @testset "column_optics with dual parameters" begin
        function column(; Chl = 0.5, T = 20.0)
            top, deep = mixed_layer(; Chl, T, Δz = 5.0), mixed_layer(; Chl = 0.1 * one(Chl), T)
            lower = OceanLayer(5.0, 20.0, deep.constituents)
            return column_optics(OceanColumn([top, lower]), λ_layer[1:4];   # Quan–Fry: 400–700 nm
                                 refractive_index = QuanFry1995(temperature = T, salinity = 35.0),
                                 bottom = LambertianBottom(0.1), ℓ_max = 8)
        end
        @test check(x -> optical_depth(column(Chl = x).layers[1]), 0.5)
        @test all(iszero, ForwardDiff.derivative(x -> optical_depth(column(Chl = x).layers[2]), 0.5))
        @test check(x -> column(T = x).refractive_index, 20.0)
        @test check(x -> absorption_optical_depth(column(T = x).layers[2]), 20.0)
        col = column(Chl = one_dual(0.5))
        @test col isa OceanColumnOptics{D}
        @test all(iszero, ForwardDiff.partials.(bottom_albedo(col), 1))
    end

    @testset "Tabulated coefficients and IOP carry dual numbers" begin
        λ_tab, a_tab = [400.0, 500.0, 600.0], [0.10, 0.20, 0.40]
        @test ForwardDiff.derivative(x -> absorption(TabulatedAbsorber(λ_tab, x .* a_tab), 500.0), 1.0) ≈ 0.20
        b(x) = TabulatedScatterer(λ_tab .* one(x), x .* a_tab, HenyeyGreensteinPhase(g = 0.8 * one(x)))
        @test ForwardDiff.derivative(x -> backscattering(b(x), 600.0), 1.0) ≈
              0.40 * backscatter_fraction(HenyeyGreensteinPhase(g = 0.8)) rtol = 1e-12
        # ω = b/(a + b) ⇒ ∂ω/∂a = −b/(a + b)²
        @test ForwardDiff.derivative(x -> single_scattering_albedo(IOP(x, 0.3, 0.01) + IOP(0.1, 0.2, 0.01)),
                                     0.1) ≈ -0.5 / 0.7^2 rtol = 1e-12
        @test IOP(1, 2, 0.5) isa IOP{Float64}
    end

    @testset "What is rejected for dual numbers" begin
        x = one_dual(0.5)
        # A model of fixed floating-point type cannot store a dual parameter.
        @test_throws MethodError Phytoplankton{Float64}(Chl = x)
        @test_throws MethodError CDOM{Float64}(a_ref = x)
        @test_throws MethodError FournierForandPhase{Float64}(backscatter_fraction = one_dual(0.018))
        # Measured phase tables have no parameter and exist in floating point only.
        @test_throws "has no parameter to differentiate" PetzoldPhase{D}()
        @test_throws ArgumentError PetzoldPhase{D}()
        # Inelastic processes are floating point only.
        @test_throws TypeError WaterRaman{D}()
        @test_throws MethodError chlorophyll_fluorescence(Phytoplankton(Chl = x); quantum_yield = 0.02)
        # A backscatter fraction outside the Mobley et al. (2002) family, also as a dual.
        @test_throws "outside the range" FournierForandPhase(backscatter_fraction = one_dual(0.7))
    end
    @testset "Range checks accept dual numbers at their bounds" begin
        # ForwardDiff 1.x orders duals by value, then partials: Dual(35, 1) > 35.
        nS(S) = refractive_index(QuanFry1995(temperature = 20, salinity = S), 500)
        # n is linear in S (Quan & Fry 1995), so the backward difference is exact.
        @test ForwardDiff.derivative(nS, 35.0) ≈ nS(35.0) - nS(34.0) rtol = 1e-9
        nT(T) = refractive_index(QuanFry1995(temperature = T, salinity = 10), 500)
        # n is quadratic in T: the three-point backward difference is exact.
        @test ForwardDiff.derivative(nT, 30.0) ≈ (3nT(30.0) - 4nT(29.0) + nT(28.0)) / 2 rtol = 1e-8
        nλ(λ) = refractive_index(QuanFry1995(temperature = 20, salinity = 35), λ)
        @test ForwardDiff.derivative(nλ, 700.0) ≈ (3nλ(700.0) - 4nλ(699.99) + nλ(699.98)) / 0.02 rtol = 1e-6
        # Fell mineral b = ω c_ref (λ_ref/λ): ∂b/∂ω = c_ref λ_ref/λ, also at ω = 1.
        bω(ω) = scattering(NonAlgalParticles{typeof(ω), FellMineral}(ω_nap = ω), 500.0)
        @test ForwardDiff.derivative(bω, 1.0) ≈ 0.40 * 440 / 500 rtol = 1e-14
        # A slope at the edge of the advisory range does not warn.
        aS(S) = absorption(CDOM(slope = S), 500.0)
        d = @test_logs ForwardDiff.derivative(aS, 0.030)
        @test d ≈ -60 * 0.03 * exp(-0.030 * 60) rtol = 1e-14
        # A wavelength at the last node of an absorption table.
        w = PureWater(absorption_model = PopeFry1997())
        @test isfinite(ForwardDiff.derivative(x -> absorption(w, x), 727.5))
    end

    @testset "Zhang-Hu scattering: ∂b/∂S at S = 0 is the one-sided derivative" begin
        b(S) = scattering(PureWater(temperature = 20.0, salinity = S), 500.0)
        d0 = ForwardDiff.derivative(b, 0.0)
        # b(S) − b(0) = d S + a S^1.5 + O(S²): D(h) = (b(h) − b(0))/h = d + a√h + O(h), and
        # 2D(h/4) − D(h) cancels the √h term.
        D(h) = (b(h) - b(0.0)) / h
        @test isfinite(d0) && d0 > 0
        @test d0 ≈ 2D(2.5e-5) - D(1e-4) rtol = 1e-4
        @test d0 ≈ ForwardDiff.derivative(b, 1e-10) rtol = 1e-4
    end

    @testset "Fournier-Forand: nested derivatives with respect to B" begin
        μj(B) = FournierForandPhase(backscatter_fraction = B).μ_j
        d1(B) = ForwardDiff.derivative(μj, B)
        d2(B) = ForwardDiff.derivative(d1, B)
        h = 1e-5
        # Richardson-extrapolated central difference of a derivative g, error O(h⁴).
        rich(g, B) = (4 * (g(B + h / 2) - g(B - h / 2)) / h - (g(B + h) - g(B - h)) / 2h) / 3
        for B in (0.005, 0.02, 0.05)
            @test d2(B) ≈ rich(d1, B) rtol = 1e-8
            @test ForwardDiff.derivative(d2, B) ≈ rich(d2, B) rtol = 1e-8
            @test ForwardDiff.hessian(x -> μj(x[1]), [B])[1] ≈ d2(B) rtol = 1e-12
            # The family relation B(n(μ_j), μ_j) = B holds to second order: its second
            # derivative along B vanishes.
            Bback(x) = backscatter_fraction(FournierForandPhase(backscatter_fraction = x))
            @test abs(ForwardDiff.derivative(y -> ForwardDiff.derivative(Bback, y), B)) < 1e-8
        end
        β2(B) = phase_function_moments(FournierForandPhase(backscatter_fraction = B), 4)[3]
        dβ2(B) = ForwardDiff.derivative(β2, B)
        for B in (0.01, 0.02)
            d = ForwardDiff.derivative(dβ2, B)
            @test isfinite(d)
            @test d ≈ rich(dβ2, B) rtol = 1e-6
        end
        # Second derivatives in n and μ_j directly.
        for (f, x0) in ((x -> phase_function_moments(FournierForandPhase(n = x, μ_j = 3.5), 4)[3], 1.1),
                        (x -> phase_function_moments(FournierForandPhase(n = 1.1, μ_j = x), 4)[3], 3.5))
            df(x) = ForwardDiff.derivative(f, x)
            @test ForwardDiff.derivative(df, x0) ≈ rich(df, x0) rtol = 1e-6
        end
    end

    @testset "Phytoplankton without an element type follows a dual scattering model" begin
        # b_p = b_ref Chl^0.62 (550/λ): ∂b_p/∂b_ref = Chl^0.62 · 550/λ.
        bp(x) = scattering(Phytoplankton(Chl = 0.5, scattering_model = GordonMorel1983(b_ref = x)), 500.0)
        @test ForwardDiff.derivative(bp, 0.3) ≈ 0.5^0.62 * 550 / 500 rtol = 1e-14
        τ(x) = optical_depth(layer_optics(OceanLayer(0.0, 2.0,
                   [PureWater{typeof(x)}(), Phytoplankton(Chl = 0.5, scattering_model = GordonMorel1983(b_ref = x))]),
                   [500.0]))[1]
        @test ForwardDiff.derivative(τ, 0.3) ≈ 2 * 0.5^0.62 * 550 / 500 rtol = 1e-12
        x = ForwardDiff.Dual(0.3, 1.0)
        @test_throws "cannot hold" Phytoplankton{Float32}(Chl = 0.5, scattering_model = GordonMorel1983(b_ref = x))
        @test Phytoplankton{Float32}(Chl = 0.5).scattering_model isa GordonMorel1983{Float64}
    end
end
