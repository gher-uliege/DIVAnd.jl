using divand
using Base.Test


@testset "divand" begin
    include("test_covaris.jl");


    include("test_sparse_diff.jl");
    include("test_localize_separable_grid.jl");
    include("test_statevector.jl");

    include("test_diagnostic_tools.jl");

    include("test_1dvar.jl");
    include("test_1D_seminormed.jl");


    include("test_2dvar_adv.jl");
	include("test_2dvar_jog.jl");
    # currently broken in v 0.6
    if VERSION == v"0.5"
        include("test_2dvar_iter.jl");
    end
    include("test_2dvar_error.jl");

    # cross-validation
    include("test_2dvar_cv.jl");
    include("test_2dvar_matfun.jl");
    include("test_2dvar_qc.jl");
    include("test_variableLandpmn.jl")

    include("test_3dvar.jl");

    include("test_4dvar.jl");
    include("test_divandgo.jl");


    # test kernel
    include("test_kernel.jl");

    # test kernel
    include("test_kernel.jl");
	
	
	include("test_averaged_bg.jl");

	
	# test divand_filter3
	A=zeros(5,5,5,5,5);A[3,3,3,3,3]=1
	z=divand_filter3(A,9999,1)
	
	@test maximum(z)≈0.00411522633744856

    # test divand_metric
    lon,lat = ndgrid([0:10;],[0:5;])
    pm,pn = divand_metric(lon,lat)
    @test 111e3 < mean(1./pm) < 112e3
    @test 111e3 < mean(1./pn) < 112e3

    # test divand_kernel
    mu,K = divand_kernel(2,[1,2,1])
    @test mu ≈ 4π
    @test K(0) ≈ 1
    @test K(1) ≈ SpecialFunctions.besselk(1,1)
end
