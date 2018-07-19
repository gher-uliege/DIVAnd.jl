"""
    Rtimesx!(coord,LS,x,Rx)

Gaussian type `R` matix in `ndim` dimensions applied to vector `x` of length
`ndata`. The Gaussian scale differs in each direction `k` : `LS[k]`
Coordinates of point i are `coord[i,1],coord[i,2],...,coord[i,ndim]`
To avoid an ndata² complexity a grid is set up first so as to allow only to calculate
covarances when distances are smaller than `3*LS`

Adapted from DIVA3D/src/Fortran/Util/Rtimesx_weighting.f90
"""
function Rtimesx!(coord,LS::NTuple{ndim,T},x,Rx) where T where ndim
    ndata = size(coord,2)
    len = [LS...]
    coordmin = Compat.minimum(coord, dims = 2)[:,1]
    coordmax = Compat.maximum(coord, dims = 2)[:,1]

    ilen = 1 ./ len
    ilenmax = 1 ./ (3*len)

    # Slightly enlarge bounding box to be sure all points remain
    # in box even when rounding occurs

    range = coordmax - coordmin
    coordmin -= range * eps(eltype(coord))
    coordmax += range * eps(eltype(coord))


    # Now number of grid points in each direction
    sz = (round.(Int, (coordmax - coordmin) .* ilenmax) .+ 1 ...,) :: NTuple{ndim,Int}

    # now allocate the arrays
    NP = zeros(Int,sz)
    NG = zeros(Int,ndim)
    gridindex = zeros(Int,ndata,ndim)

    # First dummy loop, identify the number of points which fall into
    # any bin of a regular grid

    for i = 1:ndata
        for j = 1:ndim
            NG[j] = round(Int,(coord[j,i] - coordmin[j]) * ilenmax[j]) + 1
        end
        NP[NG...] += 1
    end

    # Now we can allocate the array which indexes points that fall into the grid

    IP = Array{Vector{Int},ndim}(undef,sz)

    for i in eachindex(IP)
        IP[i] = Vector{Int}(undef,NP[i])
    end

    # For each grid point collect index all points which fall into bin

    NP[:] .= 0

    for i = 1:ndata
        for j = 1:ndim
            NG[j] = round(Int,(coord[j,i] - coordmin[j]) * ilenmax[j]) + 1
        end

        NGind = CartesianIndex{ndim}(NG...)

        NP[NGind] += 1
        NPP = NP[NGind]

        IP[NGind][NPP] = i

        # For all points get index of grid bin where if falls
        gridindex[i,:] = NG
    end

    # Ok, now finally calculate covariances and application
    Rx[:] .= 0

    for i = 1:ndata
        # Find grid indexes
        NG = gridindex[i,:]

        # Now all boxes around this one
        Rx[i] = 0

        istart = CartesianIndex( (max.( 1,NG.-1)...,) :: NTuple{ndim,Int} )
        iend   = CartesianIndex( (min.(sz,NG.+1)...,) :: NTuple{ndim,Int} )

        for ind in CartesianIndices(ntuple(i -> istart[i]:iend[i], ndim) :: NTuple{ndim,UnitRange{Int}})

            for ipoint = 1:NP[ind]
                ii = IP[ind][ipoint]

                dist = 0.
                for j = 1:ndim
                    dist += ((coord[j,i]-coord[j,ii]) * ilen[j])^2
                end

                cov = exp(-dist)
                Rx[i] += cov * x[ii]
            end
        end
    end
end



"""
     weights = weight_RtimesOne(x,len)

Compute the weight of the observations at location `x` to reduce the influence
of locally clustered data.  `x` is a tuple with n elements. Every element
represents a coordinate of the observations. `len` is a tuple of arrays
representing the correlation length. `len[i]` is the correlation length in the
i-th dimension.
"""
function weight_RtimesOne(x::NTuple{ndim,Vector{T}},len) where T where ndim
    coord = hcat(x...)'

    # geometric mean
    geomean(v) = prod(v)^(1/length(v))
    LS = ([geomean(l) for l in len]...,)

    vec = ones(T,length(x[1]))
    Rvec = ones(T,length(x[1]))

    #@code_warntype Rtimesx!(coord,LS,vec,Rvec)
    Rtimesx!(coord,LS,vec,Rvec)

    return 1 ./ Rvec
end
