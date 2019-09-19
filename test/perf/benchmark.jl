using JSON


ndim = parse(Int, ARGS[1])
method = ARGS[2]



#sizes = 10:10:20

sizes = if ndim == 2
    100:100:1000
    #10:10:20
elseif ndim == 3
    10:10:100
elseif ndim == 4
    4:4:20
end

@show sizes

memory = zeros(length(sizes))
runtime = zeros(length(sizes))
juliaexec = "julia"




# run
for i = 1:length(sizes)
    fileout = "res_$(ndim)_$(method)_$(sizes[i]).out"
    fileerr = "res_$(ndim)_$(method)_$(sizes[i]).err"
    # must be one line otherwise the result of time -v is too difficult to parse

    if method == "sparse"
        cmd = "using DIVAnd; using JSON; include(\"test_2dvar_benchmark.jl\"); print(JSON.json(benchmark_nd_repeat($(ndim),$(sizes[i]),10; operatortype=Val{:sparse})))"
    elseif method == "varanalysis"
        cmd = "using DIVAnd; using JSON; include(\"test_2dvar_benchmark.jl\"); print(JSON.json(benchmark_nd_repeat($(ndim),$(sizes[i]),5; tol = 1e-3)))"
    else
        cmd = "using DIVAnd; using JSON; include(\"test_2dvar_benchmark.jl\"); print(JSON.json(benchmark_nd_repeat($(ndim),$(sizes[i]),5; inversion=:pcg, operatortype=Val{:MatFun}, tol=1e-3, maxit = 100000)))"
    end
    run(pipeline(
        `/usr/bin/time -v $(juliaexec) --eval $(cmd)`,
        stdout = fileout,
        stderr = fileerr,
    ))
end
