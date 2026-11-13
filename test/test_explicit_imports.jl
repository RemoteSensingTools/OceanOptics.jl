using ExplicitImports

@testset "ExplicitImports" begin
    # `Base.Fix1`, `Base.require_one_based_indexing` and Adapt's extension
    # point `Adapt.adapt_structure` are documented API, declared `public`
    # from Julia 1.11 on; Julia 1.10 has no such declaration to check.
    ExplicitImports.test_explicit_imports(OceanOptics;
        all_qualified_accesses_are_public = (; ignore = (:Fix1, :require_one_based_indexing,
                                                         :adapt_structure)))
end
