@testset "examples" begin
    files = [
        "DIVAnd_simple_example_1D.jl",
        "DIVAnd_simple_example_1D_seminormed.jl",
        "DIVAnd_simple_example1p.jl",
        "DIVAnd_simple_example_cvepsandL.jl",
        "DIVAnd_simple_example_cveps.jl",
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
        "DIVAnd_testkernepc.jl",
    ]

    for file in files
        include("../examples/$file")
        close("all")
    end

end
