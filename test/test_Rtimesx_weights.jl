using Test
using Random
using Statistics
using DIVAnd
using Base.Threads

#using LibSpatialIndex
#const SI = LibSpatialIndex

coord = copy([
    0.0853508 0.939756;
    0.784134 0.080227;
    0.999551 0.784304;
    0.636594 0.7699;
    0.357327 0.891722;
    0.101827 0.856188;
    0.862349 0.0555934;
    0.992086 0.97036;
    0.702955 0.591252;
    0.685006 0.23132
]')

#x = (coord[:,1],coord[:,2])

LS = (0.1, 0.1)

dist2(x, y, len) = sum(((x - y) ./ len).^2)

ndata = size(coord, 2)
x = ones(ndata)
Rx1 = zeros(ndata)
Rx = zeros(ndata)

"""
naive version, just for comparision
"""
function Rtimesx1!(coord, LS, x, Rx)
    ndata = size(coord, 2)
    len = [LS...]
    cov = zeros(ndata, ndata)

    for j = 1:ndata
        for i = 1:ndata
            d2 = dist2(coord[:, i], coord[:, j], len)
            #cov[i, j] = exp(-d2)
            cov[i, j] = DIVAnd.approximate_gaussian(-d2)
        end
    end

    Rx[:] = cov * x
end


function Rtimesx2!(coord, len::NTuple{ndim,T}, w, Rx) where {T} where {ndim}
    factor = 3

    n = size(coord, 1)
    Nobs = size(coord, 2)

    maxcap = 10
    qt = DIVAnd.Quadtrees.QT(coord, collect(1:Nobs))::DIVAnd.Quadtrees.QT{T,Int,ndim}
    DIVAnd.Quadtrees.rsplit!(qt, maxcap, len)
    @debug "Quadtree depth: $(DIVAnd.Quadtrees.maxdepth(qt))"

    ilen = 1 ./ len

    #index_buffer = zeros(Int, Nobs)
    index_buffer_all = zeros(Int, Nobs, Threads.nthreads())

    @show Threads.nthreads()

    @inbounds Threads.@threads for i = 1:Nobs
        index_buffer = @view index_buffer_all[:,Threads.threadid()]

        xmin = ntuple(j -> coord[j, i] - factor * len[j],Val(ndim))
        xmax = ntuple(j -> coord[j, i] + factor * len[j],Val(ndim))

        nindex = DIVAnd.Quadtrees.within_buffer!(qt, xmin, xmax, index_buffer)

        Rx[i] = 0.

        #for ii in @view index_buffer[1:nindex]
        @inbounds for j = 1:nindex
            ii = index_buffer[j]

            dist2 = 0.
            for j = 1:ndim
                dist2 += ((coord[j, i] - coord[j, ii]) * ilen[j])^2
            end

            #cov = @fastmath exp(-dist2)
            cov = DIVAnd.approximate_gaussian(-dist2)
            Rx[i] += cov * w[ii]
        end
    end
end




function Rtimesx3!(coord, len::NTuple{ndim,T}, w, Rx) where {T} where {ndim}
    factor = 3

    n = size(coord, 1)
    Nobs = size(coord, 2)

    rtree = SI.RTree(n)
    for i = 1:size(coord, 2)
        SI.insert!(rtree, i, coord[:, i], coord[:, i])
    end



    ilen = 1 ./ len

    index_buffer = zeros(Int, Nobs)

    @fastmath @inbounds for i = 1:Nobs
        for j = 1:ndim
            xmin[j] = coord[j, i] - factor * len[j]
            xmax[j] = coord[j, i] + factor * len[j]
        end

        index_buffer = SI.intersects(rtree, xmin, xmax)
        nindex = length(index_buffer)

        Rx[i] = 0.

        #for ii in @view index_buffer[1:nindex]
        @inbounds for j = 1:nindex
            ii = index_buffer[j]

            dist2 = 0.
            for j = 1:ndim
                dist2 += ((coord[j, i] - coord[j, ii]) * ilen[j])^2
            end

            #cov = exp(-dist2)
            cov = DIVAnd.approximate_gaussian(-dist2)
            Rx[i] += cov * w[ii]
        end
    end
end



# compare naive and optimized method
Rtimesx1!(coord, LS, x, Rx1)
DIVAnd.Rtimesx!(coord, LS, x, Rx)

@test Rx1 ≈ Rx

Rtimesx2!(coord, LS, x, Rx)
@test Rx1 ≈ Rx rtol = 1e-4

#Rtimesx3!(coord, LS, x, Rx)
#@test Rx1 ≈ Rx rtol = 1e-4

# fix seed of random number generator
Random.seed!(12345)

# observations
# uniformly distributed data with a cluster at (0.2,0.3)

x = [rand(75); rand(75) / 10 .+ 0.2]
y = [rand(75); rand(75) / 10 .+ 0.3]

# length-scale to consider clustered data
len = (0.01, 0.01)

# compute weigths
weight = DIVAnd.weight_RtimesOne((x, y), len)

# the weight of the data points inside the cluster
# should be smaller than outside
@test mean(weight[1:75]) > mean(weight[76:end])



# "large" benchmark

# large
ndata = 70000
ndim = 2

ndata = 70000*2*2*2
ndim = 2

ndata = 30000*2*2*2*2*2
ndim = 4

coord = randn(ndim, ndata)
x = ones(ndata)
Rx = zeros(ndata)
LS = ntuple(i -> 0.1, ndim)
@time DIVAnd.Rtimesx!(coord, LS, x, Rx)

Rx2 = zeros(ndata)
@time Rtimesx2!(coord, LS, x, Rx2)

len = LS
n = size(coord, 1)
Nobs = size(coord, 2)

maxcap = 10
T = Float64
qt = DIVAnd.Quadtrees.QT(coord, collect(1:Nobs))::DIVAnd.Quadtrees.QT{T,Int,ndim}
DIVAnd.Quadtrees.rsplit!(qt, maxcap)
@show DIVAnd.Quadtrees.maxdepth(qt)

nothing

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
