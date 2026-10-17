@testset "HydroLight IOP reference (Loisel et al. 2023)" begin
    # Pure-seawater IOPs and three IOP scenarios from the HydroLight input
    # database of Loisel et al. (2023), Dryad doi:10.6076/D1630T (CC0); see
    # the provenance header of the CSV. The database stores single-precision
    # values with totals rounded to about five significant figures, so sums
    # of components close only to ~5e-5 relative; tolerances reflect that.
    ref = OceanOptics._read_csv_columns(
        OceanOptics.data_path("reference/loisel2023_hydrolight_iops.csv"))
    λ   = ref["lambda_nm"]
    @test λ == collect(350.0:5.0:750.0)

    @testset "Reference data consistency" begin
        @test ref["bbw"] ≈ ref["bw"] ./ 2 rtol = 2e-4
        for case in ("clear", "median", "turbid")
            c(v) = ref["$(case)_$v"]
            @test c("a") ≈ ref["aw"] .+ c("aph") .+ c("ag") .+ c("ad") rtol = 1e-4
            @test c("b") ≈ ref["bw"] .+ c("bph") .+ c("bd") rtol = 1e-4
            # Fournier-Forand particles with b̃b = 0.01 (phytoplankton) and
            # 0.018 (non-algal particles), as stated in the paper.
            @test c("bb") ≈ ref["bbw"] .+ 0.01 .* c("bph") .+ 0.018 .* c("bd") rtol = 1e-4
        end
    end

    @testset "Pure-seawater scattering: Zhang et al. (2009), 18 °C, S = 35" begin
        # Agreement is at the single-precision level of the reference
        # (observed ≤ 5e-5 relative over 350–750 nm).
        water = PureWater{Float64}(scattering_model = ZhangHu2009(),
                                   temperature = 18.0, salinity = 35.0)
        @test scattering(water, λ) ≈ ref["bw"] rtol = 2e-4
        @test maximum(abs, scattering(water, λ) ./ ref["bw"] .- 1) < 1e-4
    end

    @testset "Pure-water absorption" begin
        # HydroLight's tables carry no temperature correction; at the tables'
        # 22 °C reference OceanOptics applies none either.
        # Pope & Fry (1997) is the IOCCG composite over 420–725 nm.
        pf  = PureWater{Float64}(absorption_model = PopeFry1997(), temperature = 22.0)
        sel = 420.0 .≤ λ .≤ 725.0
        @test absorption(pf, λ[sel]) ≈ ref["aw"][sel] rtol = 1e-6
        # The IOCCG2018 model is this composite (bundled from the same database;
        # an identity check of the model's table switch, not an independent test).
        io = PureWater{Float64}(absorption_model = IOCCG2018(), temperature = 22.0)
        @test absorption(io, λ) ≈ ref["aw"] rtol = 1e-5
        # Kou et al. (1993), the composite's source above 730 nm, matches it at
        # 730, 740 and 750 nm within 0.5 % (within 1 % over 730–750 nm).
        kou = PureWater{Float64}(absorption_model = Kou1993(), temperature = 22.0)
        @test absorption(kou, [730.0, 740.0, 750.0]) ≈ ref["aw"][end-4:2:end] rtol = 5e-3
    end

    @testset "Tabulated route reproduces the HydroLight bulk IOPs" begin
        # Bulk a as a total absorber; water, phytoplankton and non-algal
        # scattering as tabulated scatterers with their phase functions.
        for case in ("clear", "median", "turbid")
            c(v) = ref["$(case)_$v"]
            layer = OceanLayer(0.0, 1.0, AbstractOceanConstituent{Float64}[
                TabulatedAbsorber(λ, c("a"); total = true),
                TabulatedScatterer(λ, ref["bw"], RayleighWaterPhase{Float64}()),
                TabulatedScatterer(λ, c("bph"),
                                   FournierForandPhase{Float64}(backscatter_fraction = 0.01)),
                TabulatedScatterer(λ, c("bd"),
                                   FournierForandPhase{Float64}(backscatter_fraction = 0.018))])
            opt = layer_optics(layer, λ; ℓ_max = 16)
            @test single_scattering_albedo(opt) ≈ c("b") ./ (c("a") .+ c("b")) rtol = 1e-4
            @test backscatter_fraction(opt) ≈ c("bb") ./ c("b") rtol = 2e-4
            @test optical_depth(opt) ≈ c("a") .+ c("b") rtol = 1e-4
        end
    end
end
