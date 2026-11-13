# `SpectralTable` and `linterp` against values worked out by hand on a table
# with unequal spacing: λ = (400, 410, 430, 500), y = (1, 3, 2, 9).
@testset "SpectralTable" begin
    tbl = SpectralTable([400.0, 410.0, 430.0, 500.0], [1.0, 3.0, 2.0, 9.0])
    @test tbl isa SpectralTable{Float64, Vector{Float64}}
    @test (tbl.lo, tbl.hi, length(tbl)) === (400.0, 500.0, 4)

    @testset "values" begin
        @test [linterp(tbl, x) for x in tbl.λ] == tbl.y            # knots are exact
        @test linterp(tbl, 405.0) == 2.0                           # 1 + 5 · (3 − 1)/10
        @test linterp(tbl, 425.0) ≈ 2.25 rtol = 4eps()             # 3 − 15/20
        @test linterp(tbl, 465.0) ≈ 5.5 rtol = 4eps()              # 2 + 35 · 7/70
        @test linterp(tbl, [405.0, 425.0, 500.0]) ≈ [2.0, 2.25, 9.0] rtol = 4eps()
        @test linterp(tbl, Float64[]) == Float64[]
        @test @inferred(linterp(tbl, 405.0)) isa Float64
        @test @inferred(linterp(tbl, 405.0f0)) isa Float64         # the table's precision is the floor
    end

    @testset "never beyond the data" begin
        @test_throws DomainError linterp(tbl, 399.999)
        @test_throws DomainError linterp(tbl, 500.001)
        @test_throws DomainError linterp(tbl, NaN)
        @test_throws "outside the tabulated range 400.0–500.0 nm" linterp(tbl, [405.0, 501.0])
        @test_throws "outside the tabulated range 400.0–500.0 nm" linterp(tbl, [399.0, 405.0])
    end

    # Slope of the segment. At an interior knot the interpolant has a kink:
    # the derivative is the slope of one of the two segments that meet there.
    @testset "derivative" begin
        d(x) = ForwardDiff.derivative(z -> linterp(tbl, z), x)
        @test d(405.0) ≈ 0.2 rtol = 4eps()
        @test d(420.0) ≈ -0.05 rtol = 4eps()
        @test d(499.0) ≈ 0.1 rtol = 4eps()
        @test d(400.0) ≈ 0.2 rtol = 4eps()
        @test any(s -> isapprox(d(410.0), s; rtol = 4eps()), (0.2, -0.05))
        @test any(s -> isapprox(d(430.0), s; rtol = 4eps()), (-0.05, 0.1))
    end

    @testset "construction" begin
        t32 = SpectralTable(Float32[400, 500], Float32[1, 2])
        @test t32 isa SpectralTable{Float32, Vector{Float32}}
        @test linterp(t32, 450.0f0) === 1.5f0
        @test SpectralTable(400:50:500, [1, 2, 4]) isa SpectralTable{Float64, Vector{Float64}}
        @test linterp(SpectralTable(400:50:500, [1, 2, 4]), 475) == 3.0
        @test_throws DimensionMismatch SpectralTable([400.0, 500.0], [1.0])
        @test_throws "at least two points" SpectralTable([400.0], [1.0])
        @test_throws "strictly increasing" SpectralTable([400.0, 500.0, 450.0], [1.0, 2.0, 3.0])
        @test_throws "strictly increasing" SpectralTable([400.0, 450.0, 450.0], [1.0, 2.0, 3.0])
        @test_throws "finite" SpectralTable([400.0, Inf], [1.0, 2.0])
    end

    # The shipped tables are what the models interpolate: a mid-point of two
    # neighbouring rows of the Pope & Fry (1997) file is their mean.
    @testset "shipped table" begin
        pf = OceanOptics._POPE_FRY_1997_TABLE
        i  = findfirst(==(550.0), pf.λ)
        @test i !== nothing
        @test linterp(pf, (pf.λ[i] + pf.λ[i+1]) / 2) ≈ (pf.y[i] + pf.y[i+1]) / 2 rtol = 4eps()
    end

    @testset "values must be finite" begin
        @test_throws "values must be finite" SpectralTable([400.0, 500.0], [NaN, 1.0])
        @test_throws "values must be finite" SpectralTable([400.0, 500.0], [1.0, Inf])
        @test_throws "wavelengths must be finite" SpectralTable([400.0, Inf], [1.0, 2.0])
    end
end
