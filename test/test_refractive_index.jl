@testset "Refractive index of seawater" begin
    # Independent reference: pure water at 20 °C and the sodium D line
    # (589.3 nm), n = 1.33299 relative to air (CRC Handbook).
    qf = QuanFry1995(temperature = 20, salinity = 0)
    @test refractive_index(qf, 589.3) ≈ 1.33299 atol = 5e-5
    # Salinity raises the index by roughly 2e-4 per PSU; dispersion is normal.
    qs = QuanFry1995(temperature = 20, salinity = 35)
    @test 1.5e-4 < (refractive_index(qs, 589.3) - refractive_index(qf, 589.3)) / 35 < 2.2e-4
    @test refractive_index(qs, 400.0) > refractive_index(qs, 700.0)
    @test refractive_index(qs, [450.0, 550.0]) == [refractive_index(qs, 450.0), refractive_index(qs, 550.0)]

    # Quan & Fry (1995), Eq. 4, evaluated term by term with the published
    # coefficients at a point where every term contributes (10 °C, 30 PSU, 450 nm).
    T, S, λ = 10.0, 30.0, 450.0
    n_ref = 1.31405 + (1.779e-4 - 1.05e-6 * T + 1.6e-8 * T^2) * S - 2.02e-6 * T^2 +
            (15.868 + 0.01155 * S - 0.00423 * T) / λ - 4382 / λ^2 + 1.1455e6 / λ^3
    @test refractive_index(QuanFry1995(temperature = T, salinity = S), λ) ≈ n_ref rtol = 1e-14
    @test n_ref ≈ 1.3458 atol = 1e-4

    # Fitted domain only, unless extrapolation is chosen explicitly.
    @test_throws DomainError refractive_index(qs, 350.0)
    @test_throws DomainError refractive_index(qs, 865.0)
    @test_throws DomainError QuanFry1995(temperature = 35, salinity = 35)
    @test isfinite(refractive_index(QuanFry1995(temperature = 20, salinity = 35, extrapolate = true), 865.0))

    @test refractive_index(FixedRefractiveIndex(1.34), 865.0) == 1.34
    @test_throws ArgumentError FixedRefractiveIndex(-1.0)
end
