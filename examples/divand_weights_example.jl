
using divand
using PyPlot

# fix seed of random number generator
srand(12345)

# observations
# uniformly distributed data with a cluster at (0.2,0.3)

x = [rand(75); rand(75)/10 + 0.2]
y = [rand(75); rand(75)/10 + 0.3]

# length-scale to consider clustered data
len = (0.01,0.01)

# compute weigths
weight = divand.weight_RtimesOne((x,y),len)

# Plot the results
scatter(x,y,10,weight)
colorbar()
xlabel("x")
ylabel("y")
title("weight of observations\n based on their redundancy (method 'RtimesOne')")
