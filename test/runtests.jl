using divand
using Base.Test


@testset "divand" begin
    include("test_covaris.jl");


    include("test_sparse_diff.jl");
    include("test_localize_separable_grid.jl");
    include("test_statevector.jl");

    include("test_2dvar_check.jl");
    include("test_3dvar.jl");


    lon,lat = ndgrid([0:10;],[0:5;])
    pm,pn = divand_metric(lon,lat)
    @test 111e3 < mean(1./pm) < 112e3
    @test 111e3 < mean(1./pn) < 112e3
end

