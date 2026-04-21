@testset "OceanColumn" begin

    @testset "Direct construction + contiguity validation" begin
        cdom = CDOM{Float64}(a_ref = 0.03)
        l1   = OceanLayer(0.0,  1.0, AbstractOceanConstituent{Float64}[cdom])
        l2   = OceanLayer(1.0,  2.5, AbstractOceanConstituent{Float64}[cdom])
        l3   = OceanLayer(2.5,  5.0, AbstractOceanConstituent{Float64}[cdom])
        col  = OceanColumn([l1, l2, l3])

        @test col isa OceanColumn{Float64}
        @test n_layers(col)    == 3
        @test thickness(col)   == 5.0
        @test bottom_depth(col) == 5.0
        @test midpoint_depths(col) ≈ [0.5, 1.75, 3.75]
        @test layer_thicknesses(col) ≈ [1.0, 1.5, 2.5]

        # Collection interface
        @test length(col) == 3
        @test col[1] === l1
        @test col[end] === l3
        @test collect(col) == [l1, l2, l3]

        # Surface-anchor violation
        bad_surface = OceanLayer(0.5, 1.0, AbstractOceanConstituent{Float64}[cdom])
        @test_throws ArgumentError OceanColumn([bad_surface])

        # Gap between layers
        g1 = OceanLayer(0.0, 1.0, AbstractOceanConstituent{Float64}[cdom])
        g2 = OceanLayer(1.5, 2.0, AbstractOceanConstituent{Float64}[cdom])
        @test_throws ArgumentError OceanColumn([g1, g2])

        # Empty column
        @test_throws ArgumentError OceanColumn(OceanLayer{Float64}[])
    end

    @testset "uniform_column" begin
        water = PureWater{Float64}()
        col   = uniform_column(AbstractOceanConstituent{Float64}[water];
                               n_layers = 10, depth = 50.0)

        @test n_layers(col) == 10
        @test thickness(col) == 50.0
        @test all(layer_thicknesses(col) .≈ 5.0)
        @test col[1].depth_top  == 0.0
        @test col[end].depth_bottom == 50.0

        # Contiguity is exact (helpers guarantee no FP drift)
        for i in 1:n_layers(col)-1
            @test col[i].depth_bottom == col[i+1].depth_top
        end

        # Sanity
        @test_throws ArgumentError uniform_column(AbstractOceanConstituent{Float64}[water];
                                                  n_layers = 0)
        @test_throws ArgumentError uniform_column(AbstractOceanConstituent{Float64}[water];
                                                  depth = -1.0)
    end

    @testset "fell_column matches DESIGN §4.5 grid" begin
        cdom = CDOM{Float64}(a_ref = 0.03)
        col  = fell_column(AbstractOceanConstituent{Float64}[cdom]; depth = 60.0)

        # 1 m (0–10) + 2 m (10–20) + 5 m (20–60) = 10 + 5 + 8 = 23 layers
        @test n_layers(col) == 23

        # Surface layer: 0–1 m, one-metre thick
        @test col[1].depth_top == 0.0 && col[1].depth_bottom == 1.0

        # Penultimate stage transitions
        @test col[10].depth_bottom == 10.0        # end of 1 m stage
        @test col[11].depth_top    == 10.0        # start of 2 m stage
        @test col[11].depth_bottom == 12.0
        @test col[15].depth_bottom == 20.0        # end of 2 m stage
        @test col[16].depth_top    == 20.0        # start of 5 m stage
        @test col[16].depth_bottom == 25.0

        # Deepest layer ends exactly at `depth`
        @test col[end].depth_bottom == 60.0
    end

    @testset "fell_column truncation / partial stages" begin
        cdom = CDOM{Float64}(a_ref = 0.03)

        # Entirely within the 1 m stage
        col5 = fell_column(AbstractOceanConstituent{Float64}[cdom]; depth = 5.0)
        @test n_layers(col5)    == 5
        @test thickness(col5)   == 5.0
        @test all(layer_thicknesses(col5) .≈ 1.0)

        # Stage-2 overhang: 10 + 2 × 2 + fractional 1 m → 13 layers
        col15 = fell_column(AbstractOceanConstituent{Float64}[cdom]; depth = 15.0)
        @test n_layers(col15) == 13
        @test col15[end].depth_bottom == 15.0

        # Stage-3 mid-range: depth = 25 → 10 + 5 + 1 = 16 layers
        col25 = fell_column(AbstractOceanConstituent{Float64}[cdom]; depth = 25.0)
        @test n_layers(col25) == 16
        @test col25[end].depth_bottom == 25.0
    end

    @testset "Column carries fluorophores through layers" begin
        phyto = Phytoplankton{Float64, Bricaud1998}(Chl = 0.5)
        sif   = chlorophyll_fluorescence(phyto)
        col   = uniform_column(AbstractOceanConstituent{Float64}[phyto];
                               n_layers = 5, depth = 25.0,
                               fluorophores = AbstractOceanInelasticProcess{Float64}[sif])
        @test all(length(layer.fluorophores) == 1 for layer in col)
        @test all(layer.fluorophores[1] === sif for layer in col)
    end
end
