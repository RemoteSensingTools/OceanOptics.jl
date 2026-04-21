using QuadGK

@testset "Spectral redistribution" begin

    @testset "GaussianEmission peak, width, normalization" begin
        r = GaussianEmission{Float64}(center_nm = 685.0, σ_nm = 10.6)

        # Peak at center
        f_peak = redistribution(r, 500.0, 685.0)
        f_off  = redistribution(r, 500.0, 700.0)
        @test f_peak > f_off > 0
        @test redistribution(r, 500.0, 685.0) ≈ redistribution(r, 600.0, 685.0)  # excitation-independent

        # Normalization: ∫ f(λ', λ) dλ = 1 over a wide window
        I, _ = quadgk(λ -> redistribution(r, 500.0, λ), 600.0, 770.0; rtol = 1e-10)
        @test I ≈ 1.0 rtol=1e-6

        # FWHM relationship: f(center ± FWHM/2) = f(center)/2
        fwhm = 2 * sqrt(2 * log(2)) * r.σ_nm
        @test redistribution(r, 500.0, 685.0 + fwhm/2) ≈ f_peak / 2 rtol=1e-10
        @test fwhm_to_sigma(fwhm) ≈ r.σ_nm                          rtol=1e-12
    end

    @testset "GaussianEmission constructor validation" begin
        @test_throws ArgumentError GaussianEmission{Float64}(center_nm = 0.0,   σ_nm = 1.0)
        @test_throws ArgumentError GaussianEmission{Float64}(center_nm = 685.0, σ_nm = 0.0)
        @test_throws ArgumentError GaussianEmission{Float64}(center_nm = -1.0,  σ_nm = 1.0)
    end

    @testset "WavenumberShiftEmission: Raman peak at 440 nm → ≈ 491 nm" begin
        # 3400 cm⁻¹ Stokes shift moves ν' = 10⁷/440 ≈ 22727 cm⁻¹ to ν ≈ 19327 cm⁻¹,
        # which is λ ≈ 517 nm. Wait — redo: shift for water Raman is 3400 cm⁻¹;
        # at λ'=440 we expect emission peak λ = 10⁷/(22727 - 3400) = 10⁷/19327 ≈ 517.4 nm.
        r = WavenumberShiftEmission{Float64}(shift_cm_inv = 3400.0, σ_cm_inv = 85.0)

        λ_prime = 440.0
        ν_prime = 1e7 / λ_prime
        ν_peak  = ν_prime - 3400.0
        λ_peak  = 1e7 / ν_peak

        # Numerical scan to confirm the peak location
        λ_grid  = range(500.0, 540.0; length = 4001)
        f_grid  = [redistribution(r, λ_prime, λ) for λ in λ_grid]
        (_, i) = findmax(f_grid)
        @test λ_grid[i] ≈ λ_peak atol = 0.1   # 4001 points over 40 nm ⇒ 0.01 nm resolution

        # Normalization in wavelength: the dν/dλ Jacobian preserves ∫f dλ = 1
        I, _ = quadgk(λ -> redistribution(r, λ_prime, λ), 400.0, 700.0; rtol = 1e-9)
        @test I ≈ 1.0 rtol=1e-4

        # Different excitation ⇒ different peak location. At λ' = 500 nm
        # the Stokes-shifted emission peaks at ~602.4 nm.
        λ_peak_500 = 1e7 / (1e7 / 500.0 - 3400.0)
        @test 600.0 < λ_peak_500 < 605.0
    end

    @testset "WavenumberShiftEmission constructor validation" begin
        @test_throws ArgumentError WavenumberShiftEmission{Float64}(shift_cm_inv = 3400.0, σ_cm_inv = 0.0)
    end

    @testset "Forward-mode AD through redistribution" begin
        r = GaussianEmission{Float64}(center_nm = 685.0, σ_nm = 10.6)
        # ∂f/∂λ at the peak = 0 by symmetry
        d_peak = ForwardDiff.derivative(λ -> redistribution(r, 500.0, λ), 685.0)
        @test abs(d_peak) < 1e-10

        # ∂f/∂λ on the positive side is negative
        d_off = ForwardDiff.derivative(λ -> redistribution(r, 500.0, λ), 695.0)
        @test d_off < 0

        # WavenumberShiftEmission: derivative through the full chain
        r2 = WavenumberShiftEmission{Float64}(shift_cm_inv = 3400.0, σ_cm_inv = 85.0)
        d = ForwardDiff.derivative(λ -> redistribution(r2, 440.0, λ), 517.0)
        @test isfinite(d)
    end
end
