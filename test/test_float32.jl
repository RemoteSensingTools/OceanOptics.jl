@testset "Float32" begin
    FT = Float32; λ = 550f0
    @test absorption(CDOM{FT}(), λ) isa FT
    for nap in (NonAlgalParticles{FT}(), NonAlgalParticles{FT, FellMineral}()),
        f in (absorption, scattering, backscattering)
        @test f(nap, λ) isa FT
    end
    @test refractive_index(QuanFry1995(temperature = 20f0, salinity = 35f0), λ) isa FT
    for pf in scalar_phases(FT)
        @test greek_expansion(pf, 16) isa GreekCoefficients{FT}
        @test backscatter_fraction(pf) isa FT
    end
    @test phase_function_value(PetzoldPhase{FT}(), 0.3f0) isa FT

    layer = inelastic_layer(FT); λg = collect(FT, 350:10:750)
    o = layer_optics(layer, λg; ℓ_max = 8)             # element type of the layer
    @test o isa OceanLayerOptics{FT}
    @test eltype(o.λ) === FT && all(a -> eltype(a.τ) === FT, o.absorbers) &&
          all(s -> eltype(s.τ) === FT, o.scatterers)
    @test layer_optics(Float64, layer, λg; ℓ_max = 8) isa OceanLayerOptics{Float64}
    # The element type is the first positional argument; there is no keyword for it.
    @test_throws MethodError layer_optics(layer, λg; ℓ_max = 8, FT = Float64)
    column = uniform_column(layer.constituents; n_layers = 2, depth = FT(10))
    col = column_optics(column, λg; refractive_index = FixedRefractiveIndex(FT(1.34)),
                        bottom = LambertianBottom(0.1), ℓ_max = 8)
    @test typeof(col) === OceanColumnOptics{FT, LambertianBottom{FT, Vector{FT}}, Vector{FT}, Vector{Int}}
    @test col isa OceanColumnOptics{FT, <:LambertianBottom{FT}}
    E₀ = ones(FT, length(λg))
    @test absorbed_power(o, E₀) isa Matrix{FT}
    @test excitation_optical_depth(o, 2) isa Vector{FT}
    @test inelastic_emission(o, 2, E₀, FT[500, 520]) isa Vector{FT}

    @test absorbed_par(o, E₀) isa Vector{FT}
    @test absorbed_photons(o, E₀) isa Matrix{FT}
end

# Functions that read a data table or evaluate a formula with Float64
# constants compute on the host in Float64 and round the result once. With
# parameters and arguments that are exact in both precisions the Float32
# value is therefore the correctly rounded Float64 value (`===`); where a
# Float32 parameter differs from its Float64 default the two models differ by
# that rounding and the comparison carries a tolerance in units of
# eps(Float32).
@testset "Float32 against Float64" begin
    FT = Float32; ϵ = eps(FT)
    λs = collect(FT, 352:7:748)
    rounded(f, c, c64, xs) = all(x -> f(c, x) === FT(f(c64, Float64(x))), xs)

    # The depolarization 0.039 is not exact in Float32 and enters the Zhang-Hu
    # magnitude, so the Float64 twin takes the Float32 model's stored value.
    water   = PureWater{FT}()
    water64 = PureWater{Float64}(depolarization = Float64(water.depolarization))
    @test (water.temperature, water.salinity) == (water64.temperature, water64.salinity)
    @test all(x -> isapprox(scattering(water, x), scattering(PureWater{Float64}(), Float64(x)); rtol = 2ϵ), λs)
    @test absorption(water, 550f0) isa FT && scattering(water, 550f0) isa FT
    @test rounded(absorption, water, water64, λs)
    @test rounded(scattering, water, water64, λs)
    @test rounded(backscattering, water, water64, λs)

    phyto, phyto64 = Phytoplankton{FT}(Chl = 0.5f0), Phytoplankton{Float64}(Chl = 0.5)
    @test absorption(phyto, 550f0) isa FT && scattering(phyto, 550f0) isa FT
    @test rounded(absorption, phyto, phyto64, λs)
    @test rounded(scattering, phyto, phyto64, λs)

    # The same (n, μ_j), exact in both precisions.
    ff   = FournierForandPhase{FT}(backscatter_fraction = 0.018f0)
    ff64 = FournierForandPhase{Float64}(Float64(ff.n), Float64(ff.μ_j))
    μs   = collect(FT, -0.95:0.1:0.95)
    @test phase_function_value(ff, 0.3f0) isa FT
    @test rounded(phase_function_value, ff, ff64, μs)

    # Default parameters: 2.7e-4 and 0.17 are not exact in Float32, which
    # changes the Float32 model by up to eps/2 in each.
    raman, raman64 = WaterRaman{FT}(), WaterRaman{Float64}()
    λem = FT[505, 515, 520, 525, 535]
    @test emission(raman, 440f0, 520f0) isa FT
    @test all(x -> isapprox(emission(raman, 440f0, x), emission(raman64, 440.0, Float64(x)); rtol = 2ϵ), λem)
    @test raman_peak_wavelength(raman, 440f0) === FT(raman_peak_wavelength(raman64, 440.0))

    # A wider argument widens the result; the model's precision is the floor.
    @test absorption(water, 550.0) === absorption(water64, 550.0)
    @test scattering(phyto, 550.0) isa Float64
    @test phase_function_value(ff, 0.3) isa Float64
    @test phase_function_value(RayleighWaterPhase{FT}(), 0.3) isa Float64
    @test phase_function_value(TabulatedGreekPhase(FT[1, 0.3]), 0.3) isa Float64
    @test phase_function_value(PetzoldPhase{FT}(), 0.3) isa Float64
    @test phase_function_value(RayleighWaterPhase{FT}(), FT(0.3)) isa FT
    @test phase_function_value(TabulatedGreekPhase(FT[1, 0.3]), FT(0.3)) isa FT
    @test phase_function_value(PetzoldPhase{FT}(), FT(0.3)) isa FT
    @test emission(raman, 440.0, 520.0) isa Float64

    # A dual-number argument passes through a Float32 model; the derivative
    # is compared with a central difference in Float64 (step 1e-3 nm or
    # 1e-5 in μ, truncation and round-off below 1e-7 relative).
    index = QuanFry1995(temperature = 20f0, salinity = 35f0)
    dn    = ForwardDiff.derivative(x -> refractive_index(index, x), 550.0)
    @test dn ≈ (refractive_index(index, 550.001) - refractive_index(index, 549.999)) / 0.002 rtol = 1e-6
    rayleigh = RayleighWaterPhase{FT}()
    dp       = ForwardDiff.derivative(x -> phase_function_value(rayleigh, x), 0.3)
    @test dp ≈ (phase_function_value(rayleigh, 0.30001) - phase_function_value(rayleigh, 0.29999)) / 0.00002 rtol = 1e-6

    # Layer results: optical depths are rounded per component from Float64
    # coefficients; the sums over components and the trapezoidal sum of
    # absorbed_par run in Float32 (41 wavelengths, positive terms).
    layer, layer64 = inelastic_layer(FT), inelastic_layer(Float64)
    λg = collect(FT, 350:10:750)
    o, o64 = layer_optics(layer, λg; ℓ_max = 8), layer_optics(layer64, Float64.(λg); ℓ_max = 8)
    for f in (absorption_optical_depth, scattering_optical_depth, optical_depth, single_scattering_albedo)
        @test f(o) isa Vector{FT}
        @test f(o) ≈ f(o64) rtol = 4ϵ
    end
    E₀ = ones(FT, length(λg))
    par, par64 = absorbed_par(o, E₀), absorbed_par(o64, ones(length(λg)))
    @test par isa Vector{FT}
    @test count(>(0), par64) ≥ 3                       # the other parts absorb nothing in 400–700 nm
    @test all(((x, y),) -> y == 0 ? x == 0 : isapprox(x, y; rtol = 8ϵ), zip(par, par64))
end
