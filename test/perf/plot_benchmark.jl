using PyPlot
using JSON


method = "matfun"
method = "varanalysis"

for ndim = 2:4

    sizes = [parse(Int,split(filename,r"[_.]")[4]) for filename in filter(Regex("res_$(ndim)_$(method)_.*\.out"),readdir())]
    sort!(sizes)

    memory = zeros(length(sizes))
    runtime = zeros(length(sizes))

    # load results

    for i = 1:length(sizes)
        fileout = "res_$(ndim)_$(method)_$(sizes[i]).out"
        fileerr = "res_$(ndim)_$(method)_$(sizes[i]).err"
        @show fileerr


        try
            open(fileerr) do f
                lines = readlines(f);
                results = Dict([(split(strip(line),':')...) for line in lines])
                memory[i] = parse(results["Maximum resident set size (kbytes)"]) * 1000
            end
        catch err
            @show err
            memory[i] = 0.
        end

        try
            runtime[i] = JSON.parsefile(fileout)["median"]
        catch
            runtime[i] = 0.
        end
    end

    figure()
    subplot(2,1,1)

    sel = (runtime .> 0) .& (memory .> 0)
    plot(sizes[sel],runtime[sel]/60,"-o")
    #xlabel("size of the $(ndim)D domain")
    ylabel("time [minutes]")


    subplot(2,1,2)
    plot(sizes[sel],memory[sel]/1e9,"-o")
    xlabel("size of the $(ndim)D domain")
    ylabel("memory [GB]")

    savefig("benchmark_$(ndim)d_$(method).png")

end
