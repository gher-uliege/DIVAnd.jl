@testset "examples" begin
    files = [
             # "DIVAnd_climatology_BS0.02.jl",
             # "DIVAnd_climatology_BS0.6.jl",
             # "DIVAnd_climatology_BS.jl",
             # "DIVAnd_optimizepmn1Db.jl",
             # "DIVAnd_realistic_example.jl",
             "DIVAnd_simple_example_1D.jl",
             "DIVAnd_simple_example_1D_seminormed.jl",
             "DIVAnd_simple_example1p.jl",
             #"DIVAnd_simple_example1poptimizepmn1D.jl", # long,long
             #"DIVAnd_simple_example1poptimizepmn.jl", # long
             # "DIVAnd_simple_example_4Dbn.jl",
             # "DIVAnd_simple_example_4D.jl",
             # "DIVAnd_simple_example_bc2D.jl",
             # "DIVAnd_simple_example_big3D.jl",
             # "DIVAnd_simple_example_bug2.jl",
             # "DIVAnd_simple_example_bug.jl",
             # "DIVAnd_simple_example_climatology_4Dgo.jl",
             # "DIVAnd_simple_example_climatology_4D.jl",
             "DIVAnd_simple_example_cvepsandL.jl",
             "DIVAnd_simple_example_cveps.jl",
             # "DIVAnd_simple_example_gobig3D.jl",
             # "DIVAnd_simple_example_go.jl",
             # "DIVAnd_simple_example_go_old.jl",
             "DIVAnd_simple_example.jl",
             "DIVAnd_simple_example_jog.jl",
             "DIVAnd_simple_example_qc.jl",
             "DIVAnd_simple_example_variableLandpmn.jl",
             "DIVAnd_simple_example_waexerr.jl",
             "DIVAnd_simple_example_wcpme.jl",
             "DIVAnd_testkernelpc121.jl",
             "DIVAnd_testkernelpc2x2.jl",
             "DIVAnd_testkernelpc3Dis2Dp1D.jl",
             "DIVAnd_testkernelpc3Dis2t3t1Dbis.jl",
             "DIVAnd_testkernelpc3Dis2t3t1D.jl",
             "DIVAnd_testkernelpc4Dis2t2D.jl",
             "DIVAnd_testkernelpcB3B2.jl",
             "DIVAnd_testkernelpc.jl",
             "DIVAnd_testkernelpcjog.jl",
             "DIVAnd_testkernepc.jl"
             ]

    for file in files
        include("../examples/$file")
        close("all")
    end

end
