using Base.Test


coord = [0.0853508 0.939756; 0.784134 0.080227; 0.999551 0.784304; 0.636594 0.7699; 0.357327 0.891722; 0.101827 0.856188; 0.862349 0.0555934; 0.992086 0.97036; 0.702955 0.591252; 0.685006 0.23132]'

#x = (coord[:,1],coord[:,2])



LS = (0.1,0.1)

dist2(x,y,len) = sum(((x-y)./len).^2)
ndata = size(coord,2)
x = ones(ndata)
Rx1 = zeros(ndata)
Rx = zeros(ndata)

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


function Rtimesx!(coord,LS::NTuple{ndim,T},x,Rx) where T where ndim
    ndata = size(coord,2)
    len = [LS...]
    coordmin = minimum(coord,2)
    coordmax = maximum(coord,2)

    ilen = 1./len
    ilenmax = 1./(3*len)
    
    # Slightly enlarge bounding box to be sure all points remain
    # in box even when rounding occurs

    range = coordmax - coordmin
    coordmin -= range * eps(eltype(coord))
    coordmax += range * eps(eltype(coord))


    # Now number of grid points in each direction
    nx = round.(Int, (coordmax - coordmin) .* ilenmax)+1
    nxt = (nx...) :: NTuple{ndim,Int}

    # now allocate array
    #NP = zeros(Int,nx[1],nx[2])
    NP = zeros(Int,nxt)
    #NP2 = zeros(Int,nxt)

    NG = zeros(Int,ndim)
    gridindex = zeros(Int,ndata,ndim)
    #gridindex2 = Vector{NTuple{ndim,Int}}(ndata)

    # First dummy loop, identify the maximum number of points which fall into any bin of a regular grid

    for i=1:ndata
        for j = 1:ndim
            NG[j] = round(Int,(coord[j,i]-coordmin[j]) * ilenmax[j]) + 1
        end
        NP[NG...] += 1
    end

    NPMAX = maximum(NP)
    #@show maximum(NP),maximum(NP2)
    # Now we can allocate the array which indexes points that fall into the grid

    #ip_size = (nx...,NPMAX) :: NTuple{ndim+1,Int}
    IP = zeros(Int, nx[1] , nx[2], NPMAX )
    # For each grid point collect index all points which fall into bin

    NP[:] = 0
    for i = 1:ndata
        for j = 1:ndim
            NG[j] = round(Int,(coord[j,i]-coordmin[j]) * ilenmax[j]) + 1
        end
        
        #NG = ((round.(Int,(coord[:,i]-coordmin)./(3*len))+1)...) :: NTuple{ndim,Int}
        #@show NG2
        #NP[NG[1],NG[2]] += 1
        NP[NG...] += 1
        
        NPP = NP[NG...]
        IP[NG...,NPP] = i

        # For all points get index of grid bin where if falls
        gridindex[i,:] = NG
        #gridindex2[i] = NG2
    end

    # Ok , now finally calculate covariances and application

    Rx[:] = 0

    for i=1:ndata
        # Find grid indexes
        NG = gridindex[i,:]

        # Now all boxes around this one

        Rx[i] = 0
        # for i1 = max(1,NG[1]-1):min(nx[1],NG[1]+1)
        #     for i2 = max(1,NG[2]-1):min(nx[2],NG[2]+1)
        #     for i7=1:NP[i1,i2]
        #         ii=IP[i1,i2,i7]

        # loop over all indices between ng-1 and ng+1 (but still within bounds)
        istart = CartesianIndex(max.(1,NG-1)...) :: CartesianIndex{ndim}
        iend = CartesianIndex(min.(nx,NG+1)...)  :: CartesianIndex{ndim}

        for ind in CartesianRange(istart,iend)
            #Now for each point in the box calculate contribution
            for i7=1:NP[ind]
               ii=IP[ind,i7]

                COV = 1.
                dist = 0.
                for j = 1:ndim
                    dis = (coord[j,i]-coord[j,ii]) * ilen[j]
                    dist = dist + dis^2
                end
                COV = exp(-dist)
                Rx[i] += COV * x[ii]
            end
        #end
            end
    end

    return IP

end


@time Rtimesx1!(coord,LS,x,Rx1)
@time Rtimesx!(coord,LS,x,Rx)


#@test Rx1 ≈ Rx


# large
ndata = 20000
ndim = 2
coord = randn(ndata,ndim)'
x = zeros(ndata)
Rx = zeros(ndata)
LS = ntuple(i -> 0.1,ndim)
@time Rtimesx!(coord,LS,x,Rx)

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
