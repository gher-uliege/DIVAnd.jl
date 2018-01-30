@testset "examples" begin
    files = [
             # "divand_climatology_BS0.02.jl",
             # "divand_climatology_BS0.6.jl",
             # "divand_climatology_BS.jl",
             # "divand_optimizepmn1Db.jl",
             # "divand_realistic_example.jl",
             "divand_simple_example_1D.jl",
             "divand_simple_example_1D_seminormed.jl",
             "divand_simple_example1p.jl",
             #"divand_simple_example1poptimizepmn1D.jl", # long,long
             #"divand_simple_example1poptimizepmn.jl", # long
             # "divand_simple_example_4Dbn.jl",
             # "divand_simple_example_4D.jl",
             # "divand_simple_example_bc2D.jl",
             # "divand_simple_example_big3D.jl",
             # "divand_simple_example_bug2.jl",
             # "divand_simple_example_bug.jl",
             # "divand_simple_example_climatology_4Dgo.jl",
             # "divand_simple_example_climatology_4D.jl",
             "divand_simple_example_cvepsandL.jl",
             "divand_simple_example_cveps.jl",
             # "divand_simple_example_gobig3D.jl",
             # "divand_simple_example_go.jl",
             # "divand_simple_example_go_old.jl",
             "divand_simple_example.jl",
             "divand_simple_example_jog.jl",
             "divand_simple_example_qc.jl",
             "divand_simple_example_variableLandpmn.jl",
             "divand_simple_example_waexerr.jl",
             "divand_simple_example_wcpme.jl",
             "divand_testkernelpc121.jl",
             "divand_testkernelpc2x2.jl",
             "divand_testkernelpc3Dis2Dp1D.jl",
             "divand_testkernelpc3Dis2t3t1Dbis.jl",
             "divand_testkernelpc3Dis2t3t1D.jl",
             "divand_testkernelpc4Dis2t2D.jl",
             "divand_testkernelpcB3B2.jl",
             "divand_testkernelpc.jl",
             "divand_testkernelpcjog.jl",
             "divand_testkernepc.jl"
             ]

    for file in files
        include("../examples/$file")
        close("all")
    end

end
