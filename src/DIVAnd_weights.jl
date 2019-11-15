using Base.Threads


@inline function _grid_index(coord, i, coordmin, ilenmax, sz::NTuple{ndim,Int}) where ndim
    return CartesianIndex(
        ntuple( j -> min(max(round(Int, (coord[j, i] - coordmin[j]) * ilenmax[j]) + 1,1),sz[j]), Val(ndim))
    )
end

"""
    Rtimesx!(coord,LS,x,Rx)

Gaussian type `R` matirx in `ndim` dimensions applied to vector `x` of length
`ndata`. The Gaussian scale differs in each direction `k` : `LS[k]`
Coordinates of point i are `coord[i,1],coord[i,2],...,coord[i,ndim]`
To avoid an ndata² complexity a grid is set up first so as to allow only the calculation
of covariances when distances are smaller than `3*LS`

Adapted from DIVA3D/src/Fortran/Util/Rtimesx_weighting.f90
"""
function Rtimesx!(coord, LS::NTuple{ndim,T}, x, Rx) where {T} where {ndim}
    ndata = size(coord, 2)
    len = [LS...]
    coordmin = minimum(coord, dims = 2)[:, 1]
    coordmax = maximum(coord, dims = 2)[:, 1]

    ilen = 1 ./ len
    ilenmax = 1 ./ (3 * len)

    # Slightly enlarge bounding box to be sure all points remain
    # in box even when rounding occurs

    range = coordmax - coordmin
    coordmin -= range * eps(eltype(coord))
    coordmax += range * eps(eltype(coord))

    # Now number of grid points in each direction
    sz =  let coordmax=coordmax, coordmin=coordmin, ilenmax=ilenmax
        ntuple(j -> (round(Int, (coordmax[j] - coordmin[j]) * ilenmax[j]) + 1),Val(ndim))
    end

    # now allocate the arrays
    NP = zeros(Int, sz)
    NG = zeros(Int, ndim)
    gridindex = Vector{CartesianIndex{ndim}}(undef,ndata)

    # First dummy loop, identify the number of points which fall into
    # any bin of a regular grid

    for i = 1:ndata
        NGind = _grid_index(coord, i, coordmin, ilenmax, sz)
        NP[NGind] += 1
    end

    # Now we can allocate the array which indexes points that fall into the grid

    IP = Array{Vector{Int},ndim}(undef, sz)

    for i in eachindex(IP)
        IP[i] = Vector{Int}(undef, NP[i])
    end

    # For each grid point collect index all points which fall into bin

    NP[:] .= 0

    for i = 1:ndata
        NGind = _grid_index(coord, i, coordmin, ilenmax, sz)

        NP[NGind] += 1
        NPP = NP[NGind]

        IP[NGind][NPP] = i

        # For all points get index of grid bin where if falls
        gridindex[i] = NGind
    end

    # Ok, now finally calculate covariances and application

    Threads.@threads for i = 1:ndata
    #@inbounds for i = 1:ndata
        # Find grid indexes
        NGind = gridindex[i]

        # Now all boxes around this one
        Rx[i] = 0


        @inbounds for ind in CartesianIndices(ntuple(j -> max(1,NGind[j]-1):min(sz[j],NGind[j]+1), ndim)::NTuple{
            ndim,
            UnitRange{Int},
        })

            @inbounds for ii in IP[ind]
                dist = 0.
                @inbounds for j = 1:ndim
                    dist += ((coord[j, i] - coord[j, ii]) * ilen[j])^2
                end

                cov = @fastmath exp(-dist)
                Rx[i] += cov * x[ii]
            end
        end
    end
end



"""
     weights = weight_RtimesOne(x,len)

Compute the weight of the observations at location `x` to reduce the influence
of locally clustered data. `x` is a tuple with n elements: every element
represents a coordinate of the observations. `len` is a tuple of arrays
representing the correlation length. `len[i]` is the correlation length in the
i-th dimension.
"""
function weight_RtimesOne(x::NTuple{ndim,Vector{T}}, len) where {T} where {ndim}
    coord = copy(hcat(x...)')

    # geometric mean
    geomean(v) = prod(v)^(1 / length(v))
    LS = ([geomean(l) for l in len]...,)

    vec = ones(T, length(x[1]))
    Rvec = ones(T, length(x[1]))

    Rtimesx!(coord, LS, vec, Rvec)

    return 1 ./ Rvec
end
