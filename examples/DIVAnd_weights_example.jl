
using DIVAnd
using Compat: @info, range
using PyPlot
if VERSION >= v"0.7"
    using Random
end

include("./prep_dirs.jl")

# fix seed of random number generator
if VERSION >= v"0.7.0-beta.0"
   Random.seed!(12345)
else
   srand(12345)
end

# observations
# uniformly distributed data with a cluster at (0.2,0.3)

x = [rand(75); rand(75)/10 .+ 0.2]
y = [rand(75); rand(75)/10 .+ 0.3]

# length-scale to consider clustered data
len = (0.01,0.01)

# compute weigths
weight = DIVAnd.weight_RtimesOne((x,y),len)

# Plot the results
scatter(x,y,10,weight)
colorbar()
xlabel("x")
ylabel("y")
title("weight of observations\n based on their redundancy (method 'RtimesOne')")

figname = joinpath(figdir,basename(replace(@__FILE__,r".jl$" => ".png")));
savefig(figname)
@info "Saved figure as " * figname
