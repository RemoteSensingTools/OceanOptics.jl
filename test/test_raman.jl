using QuadGK

@testset "Water Raman" begin

    @testset "Default Haltrin-Kattawar parameterization" begin
        ram = WaterRaman{Float64}()

        # Defaults: Δν = 3400 cm⁻¹, σ_ν = 85, cross_section 2.7e-4 at 488 nm, n=5.5
        @test ram.shift_cm_inv        == 3400.0
        @test ram.σ_cm_inv            == 85.0
        @test ram.cross_section_ref   == 2.7e-4
        @test ram.cross_section_λ_ref == 488.0
        @test ram.cross_section_n     == 5.5

        # Peak wavelength: λ' = 440 nm → ≈ 517.4 nm (Stokes-shifted O-H stretch)
        @test raman_peak_wavelength(ram, 440.0) ≈ 517.4 rtol=1e-3

        # λ' = 532 nm (common Nd:YAG-doubled lidar band) → ≈ 649.5 nm
        @test raman_peak_wavelength(ram, 532.0) ≈ 649.5 rtol=1e-3
    end

    @testset "excitation_absorption: λ^{-5.5} cross-section" begin
        ram = WaterRaman{Float64}()

        # At the reference wavelength the cross-section equals `cross_section_ref`.
        @test excitation_absorption(ram, 488.0) ≈ ram.cross_section_ref rtol=1e-12

        # Spectral slope: ratio at 400 vs 500 nm equals (500/400)^5.5
        ratio = excitation_absorption(ram, 400.0) / excitation_absorption(ram, 500.0)
        @test ratio ≈ (ram.cross_section_λ_ref/400.0)^5.5 /
                      (ram.cross_section_λ_ref/500.0)^5.5  rtol=1e-12

        # UV is stronger, NIR weaker
        @test excitation_absorption(ram, 400.0) > excitation_absorption(ram, 700.0)
    end

    @testset "Emission spectrum properties" begin
        ram = WaterRaman{Float64}()

        # Normalization in wavelength space (Jacobian compensation)
        λ_prime = 440.0
        I, _ = quadgk(λ -> emission(ram, λ_prime, λ), 450.0, 650.0; rtol = 1e-9)
        @test I ≈ 1.0 rtol=1e-4

        # Emission peak is at the Raman-shifted wavelength
        λ_peak = raman_peak_wavelength(ram, λ_prime)
        λ_grid = range(λ_peak - 10.0, λ_peak + 10.0; length = 2001)
        f_grid = [emission(ram, λ_prime, λ) for λ in λ_grid]
        (_, i) = findmax(f_grid)
        @test λ_grid[i] ≈ λ_peak atol = 0.05
    end

    @testset "is_isotropic: false for Cabannes-depolarized Raman" begin
        ram = WaterRaman{Float64}()
        @test !is_isotropic(ram)
    end

    @testset "excitation_range + constructor validation" begin
        ram = WaterRaman{Float64}(excitation_range = (400.0, 650.0))
        @test excitation_range(ram) == (400.0, 650.0)

        @test_throws ArgumentError WaterRaman{Float64}(σ_cm_inv = 0.0)
        @test_throws ArgumentError WaterRaman{Float64}(cross_section_ref = -1e-4)
        @test_throws ArgumentError WaterRaman{Float64}(depolarization = 1.0)
        @test_throws ArgumentError WaterRaman{Float64}(excitation_range = (700.0, 400.0))
    end

    @testset "Forward-mode AD through WaterRaman" begin
        ram = WaterRaman{Float64}()
        # ∂a^R/∂λ' = -n · a^R / λ' (analytic from λ^{-n} power law)
        for λ′ in (400.0, 500.0, 600.0)
            d_ad  = ForwardDiff.derivative(x -> excitation_absorption(ram, x), λ′)
            d_ref = -ram.cross_section_n * excitation_absorption(ram, λ′) / λ′
            @test d_ad ≈ d_ref rtol=1e-10
        end
    end
end
