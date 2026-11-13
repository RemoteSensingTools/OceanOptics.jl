# Device tests: result structs moved with `Adapt.adapt`, their validation and
# derived quantities on device arrays, and the array methods of the
# closed-form models.
#
# The array type under test is `DEVICE_ARRAY`. An including script may set
# it (test/gpu/runtests.jl sets `CuArray`); otherwise the reference backend
# JLArrays is used, with scalar indexing disallowed. JLArrays ≥ 0.5 compiles
# every broadcast through the GPU compiler; earlier versions run kernels as
# plain host code and accept what a real device rejects.

using Adapt: adapt

if !@isdefined(DEVICE_ARRAY)
    using JLArrays: JLArrays, JLArray
    JLArrays.allowscalar(false)
    const DEVICE_ARRAY = JLArray
    @testset "reference backend" begin
        @test pkgversion(JLArrays) ≥ v"0.5"
    end
end

# Two layers of water, CDOM, phytoplankton and non-algal particles with
# chlorophyll fluorescence and water Raman, polarized, Lambertian bottom.
function device_test_column(::Type{FT}) where {FT}
    phyto = Phytoplankton{FT}(Chl = FT(0.5))
    fluo  = chlorophyll_fluorescence(phyto; quantum_yield = FT(0.02))
    cs    = AbstractOceanConstituent{FT}[PureWater{FT}(), CDOM{FT}(), phyto, NonAlgalParticles{FT}()]
    layers = [OceanLayer(FT(0), FT(4), cs; inelastic = AbstractOceanInelasticProcess{FT}[fluo, WaterRaman{FT}()]),
              OceanLayer(FT(4), FT(15), cs)]
    return column_optics(FT, OceanColumn(layers), collect(FT, 400:20:700);
                         refractive_index = QuanFry1995(temperature = FT(20), salinity = FT(35)),
                         bottom = LambertianBottom(FT(0.1)), ℓ_max = 16, polarized = true,
                         scalar_phase_approximation = :depolarizing)
end

# A layer whose scattering component has two expansions selected by
# `phase_index`, and an absorber that vanishes at one wavelength.
function device_test_indexed_layer(::Type{FT}) where {FT}
    g1 = GreekCoefficients(FT[1, 0.3, 0.1])
    g2 = GreekCoefficients(FT[1, 1.5, 0.6])
    s  = ScatteringComponent("two phases", FT[0.2, 0.0, 0.4, 0.1], [g1, g2], FT[0.4, 0.1], :scalar;
                             phase_index = [1, 2, 2, 1])
    r  = ScatteringComponent("Rayleigh", FT[0.05, 0.0, 0.02, 0.01], [GreekCoefficients(FT[1, 0, 0.5])],
                             FT[0.5], :scalar)
    a  = AbsorptionComponent("absorber", :total, FT[0.3, 0.0, 0.1, 0.2])
    return OceanLayerOptics{FT}(FT[400, 500, 600, 700], 0, 2, [a], [s, r], InelasticComponent{FT}[], 3, false)
end

@testset "Device arrays ($DEVICE_ARRAY, $FT)" for FT in (Float32, Float64)
    Dev = DEVICE_ARRAY
    col = device_test_column(FT)
    dev = adapt(Dev, col)

    @testset "scalar indexing is disallowed" begin
        @test_throws "calar indexing" dev.refractive_index[1]
    end

    @testset "adapt: arrays move, the element type and host metadata stay" begin
        @test dev isa OceanColumnOptics{FT, <:LambertianBottom{FT}}
        @test dev.refractive_index isa Dev{FT, 1}
        @test dev.bottom.albedo isa Dev{FT, 1}
        for (l, lh) in zip(dev.layers, col.layers)
            @test l isa OceanLayerOptics{FT}
            @test l.λ isa Dev{FT, 1}
            @test (l.depth_top, l.depth_bottom, l.n_moments, l.polarized) ===
                  (lh.depth_top, lh.depth_bottom, lh.n_moments, lh.polarized)
            @test all(c -> c.τ isa Dev{FT, 1}, l.absorbers)
            @test [(c.label, c.part) for c in l.absorbers] == [(c.label, c.part) for c in lh.absorbers]
            for (c, ch) in zip(l.scatterers, lh.scatterers)
                @test c.τ isa Dev{FT, 1} && c.backscatter_fraction isa Dev{FT, 1}
                @test c.phase_index isa Dev{Int, 1}
                @test all(g -> g isa GreekCoefficients{FT, <:Dev} && g.β isa Dev{FT, 1} &&
                               g.α isa Dev{FT, 1}, c.greek)
                @test (c.label, c.polarization) == (ch.label, ch.polarization)
                @test c.phase_id == ch.phase_id          # content keys of the host copies
            end
            # Inelastic components are host data: the same objects.
            @test l.inelastic === lh.inelastic
            @test all(c -> c.efficiency isa Vector{FT} && c.excitation isa Vector{Int}, l.inelastic)
        end
        @test length(dev.layers[1].inelastic) == 2 && isempty(dev.layers[2].inelastic)
        # The excitation link of the adapted layer is read on the host, from host data.
        k = findfirst(c -> !isempty(c.excitation), col.layers[1].inelastic)
        @test dev.layers[1].inelastic[k].excitation == col.layers[1].inelastic[k].excitation
    end

    @testset "adapt: shared expansions stay shared, keys are kept" begin
        g_host = col.layers[1].scatterers[1].greek[1]
        @test g_host === col.layers[2].scatterers[1].greek[1]       # shared by column_optics
        @test dev.layers[1].scatterers[1].greek[1] === dev.layers[2].scatterers[1].greek[1]
        @test dev.layers[1].scatterers[1].greek[1].β === dev.layers[2].scatterers[1].greek[1].β
        distinct(c) = length(unique(objectid, [g for l in c.layers for s in l.scatterers for g in s.greek]))
        @test distinct(dev) == distinct(col)
        @test Array(dev.layers[2].scatterers[1].greek[1].β) == g_host.β
    end

    @testset "adapt: an adaptor that changes the element type is refused" begin
        Other = FT === Float32 ? Float64 : Float32
        for x in (col, col.layers[1], col.layers[1].scatterers[1], col.layers[1].absorbers[1],
                  col.layers[1].scatterers[1].greek[1], col.bottom,
                  SpectralTable(FT[400, 500], FT[1, 2]))
            @test_throws "changes the element type from $FT to $Other" adapt(Dev{Other}, x)
            @test_throws "changes the element type" adapt(Array{Other}, x)
        end
    end

    # A struct captured by a broadcast is rebuilt by the backend's kernel adaptor,
    # which must not run reductions. The structs without host metadata qualify;
    # components carry a label and stay containers.
    @testset "result structs pass into a device broadcast" begin
        xs = Dev(FT[1, 2, 3])
        g  = dev.layers[1].scatterers[1].greek[1]
        @test Array(broadcast((gg, x) -> length(gg.β) * x, Ref(g), xs)) == FT[17, 34, 51]
        b  = dev.bottom
        @test Array(broadcast((bb, x) -> length(bb.albedo) * x, Ref(b), xs)) == FT[1, 2, 3]
    end

    @testset "round trip and host no-op" begin
        back = adapt(Array, dev)
        @test typeof(back) === typeof(col)
        @test back.refractive_index == col.refractive_index
        @test back.bottom.albedo == col.bottom.albedo
        for (l, lh) in zip(back.layers, col.layers)
            @test l.λ == lh.λ
            @test all(((c, ch),) -> c.τ == ch.τ, zip(l.absorbers, lh.absorbers))
            @test all(((c, ch),) -> c.τ == ch.τ && c.greek == ch.greek && c.phase_id == ch.phase_id &&
                                    c.backscatter_fraction == ch.backscatter_fraction,
                      zip(l.scatterers, lh.scatterers))
            @test all(((c, ch),) -> c.efficiency == ch.efficiency, zip(l.inelastic, lh.inelastic))
        end
        same = adapt(Array, col)                         # host to host: the same arrays
        @test same.layers[1].λ === col.layers[1].λ
        @test same.layers[2].scatterers[1].greek[1].β === col.layers[2].scatterers[1].greek[1].β
        @test same.refractive_index === col.refractive_index
    end

    # The device runs the host's broadcasts; a device may contract a·b + c
    # into one fused operation, which changes a sum by at most one rounding.
    @testset "derived quantities on the device equal the host's" begin
        rtol = 4 * eps(FT)
        for (l, lh) in zip(dev.layers, col.layers),
            f in (absorption_optical_depth, scattering_optical_depth, optical_depth,
                  single_scattering_albedo, backscatter_fraction)
            y = f(l)
            @test y isa Dev{FT, 1}
            @test Array(y) ≈ f(lh) rtol = rtol
        end
        @test bottom_albedo(dev) isa Dev{FT, 1}
        @test Array(bottom_albedo(dev)) == bottom_albedo(col)
        black = adapt(Dev, column_optics(FT, uniform_column(AbstractOceanConstituent{FT}[PureWater{FT}()];
                                                            n_layers = 1, depth = FT(5)), FT[450, 550];
                                         refractive_index = FixedRefractiveIndex(FT(1.34)),
                                         bottom = BlackBottom(), ℓ_max = 4))
        @test black isa OceanColumnOptics{FT, BlackBottom}
        @test bottom_albedo(black) isa Dev{FT, 1} && Array(bottom_albedo(black)) == zeros(FT, 2)

        # phase_index gather, zero optical depth (ϖ = 0, B = 0 there)
        lh = device_test_indexed_layer(FT)
        l  = adapt(Dev, lh)
        @test l.scatterers[1].phase_index isa Dev{Int, 1}
        @test Array(single_scattering_albedo(l)) ≈ FT[0.25 / 0.55, 0, 0.42 / 0.52, 0.11 / 0.31] rtol = rtol
        @test Array(backscatter_fraction(l)) ≈
              FT[(0.2 * 0.4 + 0.05 * 0.5) / 0.25, 0, (0.4 * 0.1 + 0.02 * 0.5) / 0.42,
                 (0.1 * 0.4 + 0.01 * 0.5) / 0.11] rtol = rtol
        @test Array(backscatter_fraction(l)) ≈ backscatter_fraction(lh) rtol = rtol
        @test Array(single_scattering_albedo(l))[2] === zero(FT)
    end

    @testset "construction and validation on device arrays" begin
        @test AbsorptionComponent("w", :total, Dev(FT[0.1, 0.2])).τ isa Dev{FT, 1}
        @test LambertianBottom(Dev(FT[0.1, 0.3])).albedo isa Dev{FT, 1}
        g = GreekCoefficients(Dev(FT[1, 0.3, 0.1]))
        @test g isa GreekCoefficients{FT, <:Dev}
        @test depolarized(g).α isa Dev{FT, 1} && Array(depolarized(g).α) == zeros(FT, 3)
        @test Array(resize_moments(g, 4).β) == FT[1, 0.3, 0.1, 0, 0]
        @test phase_id(g) == phase_id(GreekCoefficients(FT[1, 0.3, 0.1]))
        @test_throws "|β_1| = 4.0 > 2ℓ+1 = 3" GreekCoefficients(Dev(FT[1, 4]))
        @test_throws "β₀ = 0.5" GreekCoefficients(Dev(FT[0.5, 0.1]))
        @test_throws "γ_1 = 0.2 must be zero" GreekCoefficients(Dev(FT[1, 0.3, 0.1]); α = Dev(FT[0, 0, 0.1]),
            γ = Dev(FT[0, 0.2, 0.1]), δ = Dev(FT[0, 0, 0.1]), ϵ = Dev(FT[0, 0, 0.1]), ζ = Dev(FT[0, 0, 0.1]))
        @test_throws "finite and non-negative" AbsorptionComponent("w", :total, Dev(FT[0.1, -0.2]))
        @test_throws "albedo must lie in [0, 1]" LambertianBottom(Dev(FT[0.1, 1.3]))
        lh = device_test_indexed_layer(FT)
        l  = adapt(Dev, lh)
        @test_throws "strictly increasing or strictly decreasing" OceanLayerOptics{FT}(
            Dev(FT[400, 600, 500, 700]), 0, 2, l.absorbers, l.scatterers, l.inelastic, 3, false)
        @test_throws "phase_index refers to missing expansions" ScatteringComponent(
            "two phases", l.scatterers[1].τ, l.scatterers[1].greek, l.scatterers[1].backscatter_fraction,
            :scalar; phase_index = Dev([1, 3, 2, 1]))
    end

    # A table moved to the device is interpolated there in one broadcast; the
    # range check uses the limits stored as numbers and throws on the host.
    @testset "SpectralTable on the device" begin
        host = SpectralTable(FT[400, 410, 430, 500], FT[1, 3, 2, 9])
        tbl  = adapt(Dev, host)
        @test tbl isa SpectralTable{FT, <:Dev{FT, 1}}
        @test (tbl.lo, tbl.hi) === (host.lo, host.hi)
        @test adapt(Array, host).λ === host.λ
        λh = FT[400, 405, 410, 425, 430, 465, 500]
        y  = linterp(tbl, Dev(λh))
        @test y isa Dev{FT, 1}
        @test Array(y) ≈ FT[1, 2, 3, 2.25, 2, 5.5, 9] rtol = 4 * eps(FT)
        @test Array(y) ≈ linterp(host, λh) rtol = 4 * eps(FT)
        @test_throws "outside the tabulated range 400.0–500.0 nm" linterp(tbl, Dev(FT[405, 501]))
        @test_throws DomainError linterp(tbl, Dev(FT[399, 405]))
    end

    # Models without tables are plain bits and their array methods are one
    # `map`. The device's `exp` and `^` differ from the host's by a few ulp.
    @testset "array methods of closed-form models" begin
        λh = collect(FT, 400:12.5:700)
        λd = Dev(λh)
        index = QuanFry1995(temperature = FT(20), salinity = FT(35))
        for (f, c) in ((absorption, CDOM{FT}()), (absorption, NonAlgalParticles{FT}()),
                       (scattering, NonAlgalParticles{FT}()), (backscattering, NonAlgalParticles{FT}()),
                       (absorption, NonAlgalParticles{FT, FellMineral}()),
                       (scattering, NonAlgalParticles{FT, FellMineral}()),
                       (refractive_index, index), (refractive_index, FixedRefractiveIndex(FT(1.34))))
            @test isbits(c)
            y = f(c, λd)
            @test y isa Dev{FT, 1}
            @test Array(y) ≈ f(c, λh) rtol = 32 * eps(FT)
            @test f(c, λh) == [f(c, x) for x in λh]
        end
    end
end
