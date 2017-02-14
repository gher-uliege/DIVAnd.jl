using PyPlot
using JSON


ndim = 3
#sizes = 10:10:20

sizes = if ndim == 2
    100:100:1000
elseif ndim == 3
    10:10:100
elseif ndim == 4
    5:5:20
end

@show sizes

memory = zeros(length(sizes))
runtime = zeros(length(sizes))
juliaexec = "julia"

# run
for i = 1:length(sizes)
    fileout = "res_$(ndim)_$(sizes[i]).out"
    fileerr = "res_$(ndim)_$(sizes[i]).err"
    # must be one line otherwise the result of time -v is too difficult to parse
    cmd = "using divand; using JSON; include(\"test_2dvar_benchmark.jl\"); print(JSON.json(benchmark_nd_repeat($(ndim),$(sizes[i]),10)))"

    run(pipeline(`/usr/bin/time -v $(juliaexec) --eval $(cmd)`, stdout=fileout, stderr=fileerr))
end

# load results

for i = 1:length(sizes)
    fileout = "res_$(ndim)_$(sizes[i]).out"
    fileerr = "res_$(ndim)_$(sizes[i]).err"

    open(fileerr) do f
        lines = readlines(f);
        results = Dict([(split(strip(line),':')...) for line in lines])
        memory[i] = parse(results["Maximum resident set size (kbytes)"])
    end

    runtime[i] = JSON.parsefile(fileout)["median"]
end

subplot(2,1,1)
plot(sizes,runtime/60,"-o")
xlabel("size of the $(ndim)D domain")
ylabel("time [minutes]")


subplot(2,1,2)
plot(sizes,memory/1e6,"-o")
xlabel("size of the $(ndim)D domain)")
ylabel("memory [MB]")

savefig("benchmark_$(ndim)d.png")
