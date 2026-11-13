using Aqua

@testset "Aqua" begin
    Aqua.test_all(OceanOptics)
    # Not part of `test_all`; `Base.Docs.undocumented_names` exists from Julia 1.11.
    VERSION ≥ v"1.11" && Aqua.test_undocumented_names(OceanOptics)
end
