using Base.Test
using divand

coord = [0.0853508 0.939756; 0.784134 0.080227; 0.999551 0.784304; 0.636594 0.7699; 0.357327 0.891722; 0.101827 0.856188; 0.862349 0.0555934; 0.992086 0.97036; 0.702955 0.591252; 0.685006 0.23132]'

#x = (coord[:,1],coord[:,2])

LS = (0.1,0.1)

dist2(x,y,len) = sum(((x-y)./len).^2)

x = ones(ndata)
Rx1 = zeros(ndata)
Rx = zeros(ndata)

"""
naive version, just for comparision
"""
function Rtimesx1!(coord,LS,x,Rx)
    ndata = size(coord,2)
    len = [LS...]
    cov = zeros(ndata,ndata)

    for j = 1:ndata
        for i = 1:ndata
            d2 = dist2(coord[:,i],coord[:,j],len)
            cov[i,j] = exp(-d2)
        end
    end

    Rx[:] = cov*x
end



# compare naive and optimized method
Rtimesx1!(coord,LS,x,Rx1)
divand.Rtimesx!(coord,LS,x,Rx)

@test Rx1 ≈ Rx


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

# the weight of the data points inside the cluster
# should be smaller than outside
@test mean(weight[1:75]) > mean(weight[76:end])

# "large" benchmark

# # large
# ndata = 20000
# ndim = 2

# coord = randn(ndata,ndim)'
# x = zeros(ndata)
# Rx = zeros(ndata)
# LS = ntuple(i -> 0.1,ndim)
# @time Rtimesx!(coord,LS,x,Rx)

# large 2D
#=
  memory estimate:  18.13 MiB
  allocs estimate:  418500
  --------------
  minimum time:     1.170 s (0.48% GC)
  median time:      1.201 s (0.46% GC)
  mean time:        1.204 s (0.37% GC)
  maximum time:     1.245 s (0.44% GC)
  --------------
  samples:          5
  evals/sample:     1

=#
