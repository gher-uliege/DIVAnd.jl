# Perform a DIVAnd benchmark
# A 2d-variational analysis is performed over a square domain
# as described in http://www.geosci-model-dev.net/7/225/2014/gmd-7-225-2014.pdf
# with different number of grid points (100 x 100, 200 x 200, ... 600 x 600).
# For every domain size, the benchmark is run 10 times and the median time is
# computed.
#
# Input (optional):
#   name: if name is present then the results will be saved in the file called
#     name
#
# Output:
#   median_time: the median of the run-time in seconds
#   ng: number of grid points along one dimension
#   time: 2-d array with all results
#
# To run this benchmark you need to install DIVAnd, for example
#
# pkg install -forge DIVAnd
# pkg load DIVAnd

# Alexander Barth
# GPLv2 or later

function test_2dvar_benchmark(name)

    @printf("Running DIVAnd benchmark in 2 dimensions\n")

    # domain sizes
    ng = 100:100:500

    time = zeros(10, length(ng))
    RMS = zeros(10, length(ng))

    for j = 1:10
        for i = 1:length(ng)
            time[i, j], RMS[i, j] = benchmark2d(ng[i])
            if (RMS[i, j] > 0.2)
                error("unexpected large RMS. Results might be wrong")
            end

            @printf("size %5d time %10.4f \n", ng[i], time[i, j])
        end
    end

    @printf("\nMedian results\n")

    median_time = median(time, 2)

    for i = 1:length(ng)
        @printf("size %5d time %10.4f \n", ng[i], median_time[i])
    end

    fname = "test_2dvar_benchmark_$(name).mat"
    @printf("save result in file %s\n", fname)
    matwrite(fname, Dict("time" => time, "RMS" => RMS, "ng" => ng))

    return median_time, ng, time
end

function benchmark_nd_repeat(n, ng, ntimes; kwargs...)
    times = zeros(ntimes)
    RMS = zeros(ntimes)

    times[1], RMS[1] = benchmark_nd(n, ng; kwargs...)

    for i = 1:ntimes
        times[i], RMS[i] = benchmark_nd(n, ng; kwargs...)
    end

    mad(x) = median(abs.(x - median(x)))

    stat = Dict{String,Any}([(string(f), f(times)) for f in [
        mean,
        std,
        median,
        mad,
        minimum,
        maximum,
    ]])
    stat["samples"] = length(times)
    stat["times"] = times
    stat["RMS"] = RMS
    return stat
end


function benchmark_nd(n, ng; kwargs...)
    mg = max(ceil(Int, ng / 5), 2)
    len = 10 / ng

    epsilon2 = 0.05

    f(xy...) = .*([cos.(2 * pi * ng * x / 20) for x in xy]...)

    # grid of background field
    mask, pmn, xyi = DIVAnd_squaredom(n, linspace(0, 1, ng))
    vi = f(xyi...)

    # grid of observations
    xy = ndgrid([linspace(1e-6, 1 - 1e-6, mg) for i = 1:n]...)
    v = f([x[:] for x in xy]...)

    t1 = time_ns()
    va, s = varanalysis(mask, pmn, xyi, xy, v, len, epsilon2; kwargs...)
    #va,s = DIVAndrun(mask,pmn,xyi,xy,v,len,epsilon2; kwargs...);
    t2 = time_ns()
    time = (t2 - t1) / 1e9
    RMS = rms(va, vi)

    return time, RMS
end


function rms(x, y)

    d = x - y

    m = .!isnan.(d)
    r = mean(d[m].^2)

    r = sqrt.(r)
    return r
end
