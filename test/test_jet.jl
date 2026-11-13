using JET

@testset "JET" begin
    tm = (OceanOptics,)
    @testset "package" begin
        JET.test_package(OceanOptics; target_modules = tm)
    end
    @testset "scalar IOPs: no runtime dispatch, no possible error ($FT)" for FT in (Float64, Float32)
        for c in iop_constituents(FT), f in (absorption, scattering, backscattering, absorption_by_part)
            JET.test_opt(f, (typeof(c), FT); target_modules = tm)
            JET.test_call(f, (typeof(c), FT); target_modules = tm)
        end
        JET.test_opt(refractive_index, (QuanFry1995{FT}, FT); target_modules = tm)
    end
    @testset "expansions ($FT)" for FT in (Float64, Float32)
        for pf in scalar_phases(FT)
            JET.test_opt(greek_expansion, (typeof(pf), Int); target_modules = tm)
            JET.test_call(greek_expansion, (typeof(pf), Int); target_modules = tm)
        end
    end
    # Layer assembly iterates an abstractly typed constituent vector (dynamic dispatch by
    # design), so only possible errors are tested, not optimization.
    @testset "layer assembly: no possible error ($FT)" for FT in (Float64, Float32)
        L = typeof(inelastic_layer(FT))
        for approximation in (Nothing, Symbol)
            JET.test_call(OceanOptics._layer_optics,
                          (L, Vector{FT}, Int, Bool, approximation, Type{FT},
                           typeof(OceanOptics._expansion_cache(FT)));
                          target_modules = tm)
        end
    end
end
