using divand
using Base.Test


@testset "divand" begin
    include("test_covaris.jl");
    

    include("test_sparse_diff.jl");
    include("test_localize_separable_grid.jl");
    include("test_statevector.jl");

    include("test_2dvar_check.jl");
    include("test_3dvar.jl");
    
end

