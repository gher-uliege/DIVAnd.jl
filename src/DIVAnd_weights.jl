using Base.Threads


# https://web.archive.org/web/20191115112938/https://codingforspeed.com/using-faster-exponential-approximation/
@inline function approximate_exp(x::T) where T
    x = T(1) + x / T(256)
    x *= x; x *= x; x *= x; x *= x;
    x *= x; x *= x; x *= x; x *= x;
    return max(x,T(0))
end

"""
    x2 is x squared
"""
function approximate_gaussian(x2)
    return @fastmath 1 / (1 + x2 + x2*x2)
end

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

                #cov = @fastmath exp(-dist)
                #cov = @fastmath approximate_exp(-dist)
                cov = approximate_gaussian(dist)

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


function weight_RtimesOne_binning(x, len)
    n = length(x)

    # grid finer than a factor of 10
    dx = ntuple(i -> len[i]/10,Val(n))

    gridx = ntuple(i -> minimum(x[i]):dx[i]:maximum(x[i])+dx[i], Val(n))
    sz = length.(gridx)

    v = ones(size(x[1]))
    m2, count2, vb2, nout2 = binning(gridx, x, v)

    mask = trues(sz)
    pmn = ntuple(i -> ones(sz)/dx[1], Val(n))

    c0 = Float64.(count2);
    c = zeros(size(c0))

    # adusted length
    coef = sqrt(2)
    len_adjusted = ntuple(i -> coef * len[i], Val(n))

    # "diffusion" coefficient
    nu = ntuple(i -> fill(len_adjusted[i].^2,sz),Val(n))

    # compute inverses of cell volumne and staggered scaled coefficients
    ivol, nus = DIVAnd.DIVAnd_laplacian_prepare(mask,pmn,nu)

    # maximum allowed time step
    α0 = 1 / (2 * sum(ntuple(i -> maximum(pmn[i].^2 .* nu[i]),Val(n))))

    # 10% safety margin
    α = α0 / 1.1

    # number of iterations 1/(2*α) (rounded)
    nmax = round(Int, 1 / (2 * α))

    # 4* L² α*nmax ≈ 2 L² = L'²
    @debug "α0: $α0, α: $α, nmax: $nmax"

    # ∂c/∂t =  ∇ ⋅ (D ∇ c)
    # G(x,x',t) = det(D)^(-½) (4π t)^(-n/2)  exp( - (x -x')ᵀ D⁻¹ (x -x')ᵀ / (4t))

    # G(x,x',t) = det(D)^(-½) (4π t)^(-n/2)  exp( - (x -x')ᵀ D⁻¹ (x -x')ᵀ / (4t))

    @debug "sum(c0): $(sum(c0))"

    DIVAnd.diffusion!(ivol, nus, α, nmax, c0, c)

    @debug "sum(c): $(sum(c))"

    detD = prod(len_adjusted.^2)
    t = α * nmax
    c = c * sqrt((4π * t)^n / detD)

    @debug "range of c: $(extrema(c))"

    itp = LinearInterpolation(gridx,c,extrapolation_bc = NaN);
    ci = itp.(x...)

    weighti = 1 ./ ci
    clamp!(weighti, 0, 1)

    return weighti
end

