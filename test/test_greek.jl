@testset "GreekCoefficients and angular scattering matrices" begin

    # Rayleigh scattering without depolarization (Hansen & Travis 1974):
    # Greek coefficients in vSmartMOM's convention (`get_greek_rayleigh(0)`)
    # and the standard angular scattering matrix (F12 < 0 at 90°).
    rayleigh = GreekCoefficients([1.0, 0.0, 0.5]; α = [0.0, 0.0, 3.0],
                                 γ = [0.0, 0.0, sqrt(1.5)], δ = [0.0, 1.5, 0.0],
                                 ϵ = zeros(3), ζ = zeros(3))
    rayleigh_F(μ) = (F11 = 0.75 .* (1 .+ μ .^ 2), F12 = -0.75 .* (1 .- μ .^ 2),
                     F22 = 0.75 .* (1 .+ μ .^ 2), F33 = 1.5 .* μ, F34 = zero(μ), F44 = 1.5 .* μ)

    @testset "Validation" begin
        @test !is_polarized(GreekCoefficients([1.0, 0.5]))
        @test is_polarized(rayleigh) && n_moments(rayleigh) == 3
        @test_throws ArgumentError GreekCoefficients([1.0, 0.5]; α = [0.0, 0.0])
        @test_throws DimensionMismatch GreekCoefficients([1.0, 0.0, 0.5]; α = zeros(2), γ = zeros(3),
                                                         δ = zeros(3), ϵ = zeros(3), ζ = zeros(3))
        @test_throws ArgumentError GreekCoefficients([0.98, 0.5])
        @test_throws ArgumentError GreekCoefficients([1.0, NaN])
        # Doubled (2ℓ+1) factor for HG g = 0.8: β₁ = 7.2 > 3.
        @test_throws ArgumentError GreekCoefficients([(2ℓ + 1)^2 * 0.8^ℓ for ℓ in 0:3])
        # Row-shifted Rayleigh table (α₀ = 3, γ₁ = √1.5).
        @test_throws ArgumentError GreekCoefficients([1.0, 0.0, 0.5]; α = [3.0, 0.0, 0.0],
                                                     γ = [0.0, sqrt(1.5), 0.0], δ = [0.0, 1.5, 0.0],
                                                     ϵ = zeros(3), ζ = zeros(3))
    end

    @testset "Operations" begin
        @test resize_moments(rayleigh, 1).β == [1.0, 0.0]
        @test resize_moments(rayleigh, 5).α == [0.0, 0.0, 3.0, 0.0, 0.0, 0.0]
        d = depolarized(GreekCoefficients([1.0, 2.4]))
        @test is_polarized(d) && all(iszero, (d.α, d.γ, d.δ, d.ϵ, d.ζ))
        hg(g, L) = GreekCoefficients([(2ℓ + 1) * g^ℓ for ℓ in 0:L])
        m = mix_greek([1.0, 3.0], [hg(0.2, 4), hg(0.6, 4)])
        @test m.β ≈ 0.25 .* hg(0.2, 4).β .+ 0.75 .* hg(0.6, 4).β
        @test_throws ArgumentError mix_greek([1.0, 1.0], [hg(0.2, 4), hg(0.2, 5)])
        @test_throws ArgumentError mix_greek([0.0, 0.0], [hg(0.2, 4), hg(0.6, 4)])
        # Moment-series backscatter is exact for the expansion: HG closed form.
        g = 0.5
        @test backscatter_fraction(hg(g, 400)) ≈
              (1 - g^2) / (2g) * (1 / sqrt(1 + g^2) - 1 / (1 + g)) rtol = 1e-12
        @test asymmetry_parameter(hg(0.7, 3)) ≈ 0.7
    end

    @testset "evaluate_scattering_matrix: standard convention" begin
        μ = collect(-1.0:0.1:1.0)
        F = evaluate_scattering_matrix(rayleigh, μ)
        ref = rayleigh_F(μ)
        for k in keys(ref)
            @test getfield(F, k) ≈ getfield(ref, k) atol = 1e-12
        end
        @test F.F12[11] ≈ -0.75                     # μ = 0
        @test keys(evaluate_scattering_matrix(GreekCoefficients([1.0, 0.5]), μ)) == (:F11,)
    end

    @testset "TabulatedMuellerPhase: Rayleigh table" begin
        θ  = collect(0.0:1.0:180.0)
        pf = TabulatedMuellerPhase(θ; rayleigh_F(cosd.(θ))...)
        @test is_polarizable(pf)
        @test abs(pf.normalization - 1) < 1e-6
        g = greek_expansion(pf, 6; polarized = true)
        @test maximum(maximum(abs.(getfield(g, f) .- getfield(resize_moments(rayleigh, 6), f)))
                      for f in fieldnames(GreekCoefficients)) < 1e-6
        @test g ≈ greek_expansion(RayleighWaterPhase{Float64}(depolarization = 0.0), 6;
                                     polarized = true) atol = 1e-6
        p32 = TabulatedMuellerPhase(Float32.(θ); (; (k => Float32.(v) for (k, v) in pairs(rayleigh_F(cosd.(θ))))...)...)
        @test greek_expansion(p32, 6; polarized = true).α[3] ≈ 3 rtol = 1e-5
    end

    @testset "TabulatedMuellerPhase: Siewert (2000) aerosol, vSmartMOM oracle" begin
        # vSmartMOM's reconstruct_phase evaluated the Siewert PROBLEM IIA
        # coefficients (vSmartMOM convention) on a 1° grid; its f12 and f34
        # have the opposite sign to the standard convention (f12 > 0 for
        # Rayleigh). Expanding the standard-convention table must return the
        # published coefficients with vSmartMOM's flips γ = −row3, ϵ = −row5.
        t  = OceanOptics._read_csv_columns(joinpath(@__DIR__, "data", "siewert_IIA_scattering_matrix.csv"))
        pf = TabulatedMuellerPhase(t["theta_deg"]; F11 = t["f11"], F12 = .-t["f12"],
                                   F22 = t["f22"], F33 = t["f33"], F34 = .-t["f34"], F44 = t["f44"])
        g  = greek_expansion(pf, 11; polarized = true)
        siewert = (
            α = [0.0, 0.0, 3.726079, 2.202868, 1.190694, 0.391203, 0.105556, 0.020484, 0.003097, 0.000366, 3.5e-5, 3.0e-6],
            β = [1.0, 2.104031, 2.095158, 1.414939, 0.703593, 0.235001, 0.064039, 0.012837, 0.00201, 0.000246, 2.4e-5, 2.0e-6],
            γ = .-[0.0, 0.0, -0.116688, -0.20937, -0.227137, -0.144524, -0.05264, -0.0124, -0.002093, -0.000267, -2.7e-5, -2.0e-6],
            δ = [0.915207, 2.095727, 2.008624, 1.436545, 0.706244, 0.238475, 0.056448, 0.009703, 0.001267, 0.00013, 1.1e-5, 1.0e-6],
            ϵ = .-[0.0, 0.0, 0.065456, 0.221658, 0.097752, 0.052458, 0.009239, 0.001411, 0.000133, 1.1e-5, 1.0e-6, 0.0],
            ζ = [0.0, 0.0, 3.615946, 2.240516, 1.139473, 0.365605, 0.082779, 0.013649, 0.001721, 0.000172, 1.4e-5, 1.0e-6])
        # Published to 6 decimals; the 1° table and its interpolation add < 1e-6.
        for f in keys(siewert)
            @test maximum(abs.(getfield(g, f) .- getfield(siewert, f))) < 2e-6
        end
    end

    @testset "TabulatedMuellerPhase: input checks" begin
        θ = collect(0.0:1.0:180.0); F = rayleigh_F(cosd.(θ))
        @test_throws ArgumentError TabulatedMuellerPhase(θ[2:end]; (; (k => v[2:end] for (k, v) in pairs(F))...)...)
        @test_throws ArgumentError TabulatedMuellerPhase(θ; merge(F, (F12 = -2 .* F.F11,))...)  # |F12| > F11
        # A table of F/(4π) (normalization 1/(4π)) is caught unless explicitly rescaled.
        Fs = (; (k => v ./ 4π for (k, v) in pairs(F))...)
        @test_throws ArgumentError TabulatedMuellerPhase(θ; Fs...)
        @test TabulatedMuellerPhase(θ; Fs..., normalize = true).normalization ≈ 1 / 4π rtol = 1e-6
        @test_throws ArgumentError greek_expansion(TabulatedGreekPhase([1.0, 0.0, 0.5]), 4; polarized = true)
    end

    @testset "WavelengthDependentPhase" begin
        pfs = [HenyeyGreensteinPhase{Float64}(g = g) for g in (0.3, 0.8)]
        pf  = WavelengthDependentPhase([450.0, 650.0], pfs)
        @test is_wavelength_dependent(pf) && !is_polarizable(pf)
        @test phase_function_moments(pf, 2, 650.0) ≈ [1.0, 2.4, 3.2]
        @test backscatter_fraction(pf, 450.0) == backscatter_fraction(pfs[1])
        @test_throws ArgumentError phase_function_moments(pf, 2, 550.0)
        @test_throws ArgumentError phase_function_moments(pf, 2)
        @test_throws ArgumentError WavelengthDependentPhase([650.0, 450.0], pfs)
        @test_throws ArgumentError WavelengthDependentPhase([450.0, 650.0], [pf, pf])
        rw = WavelengthDependentPhase([450.0, 650.0],
                                      [RayleighWaterPhase{Float64}(depolarization = ρ) for ρ in (0.0, 0.09)])
        @test is_polarizable(rw)
        @test greek_expansion(rw, 2, 650.0; polarized = true).γ[3] <
              greek_expansion(rw, 2, 450.0; polarized = true).γ[3]
    end
end
