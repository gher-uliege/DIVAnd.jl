"""
Computes superobservations to reduce the number of observation point.
IMPORTANT: assumes all points are located in the domain of interest.
Also, if topology (not taken into account in the present version) is complicated, a superobs could fall on land.
In that case possible solution: take as coordinates those of the closest real observation (in water).
Also consider that when using superobs for heatmap calculations, the automatic calculation of bandwidth can become too optimistic (and yield to small bandwidth).

newx,newval,sumw,varp,idx=DIVAnd_superobs(x,val,nmax;weights=[],intensive=true)

# Input:
* `x`   : tuple of coordinates of the original observations
* `val` : array which contains the values of the observation
* `nmax`: approximate number of superobservations to calculate (real number will depend on distribution in space)
* `weights` : array of weight used to calculate weighted averages (including when calculating positions of the superobs)
*               Default is unit weights.
* `intensive` : boolean which when true (default) calculates the average of values, otherwise the sum

# Output:
* `newx` : tuple of coordinates of the new super observations
* `newval`: array of values of the new super observations
* `sumw`: array containing for each superobservation value the sum of weights that where involved. Can be used for defining weight of superobservation
* `varp`: array of variance (variance around the mean defining the new super observation). Can be used for defining weight of superobservation
* `idx` : array containing the inverse of grid spacing (one in each direction) that were used for the binning

Normal use with default values creates classical superobservations

Alternative use to prepare heatmaps:

    newx,newinflation,sumw,varp=DIVAnd_superobs(x,ones(size(inflation)),nmax;weights=inflation,intensive=false)

creates a new inflation array by summing up the inflation from the points involved in the superobservation calculation.
For the moment in this case sumw, varp are not meaningful.


"""
function DIVAnd_superobs(x, val, nmax; weights = [], intensive = true)
#
    # Array containing coordinates
    coord = hcat(x...)'
    # Number of data
    ndata = size(coord, 2)
    # Dimension
    ndim = size(coord, 1)
    # how many cells in each direction for binning
    NPD = Int(ceil(nmax^(1.0 / ndim)))
    # Calculate coordinate range of observations

    coordmin = Compat.minimum(coord, dims = 2)[:, 1]
    coordmax = Compat.maximum(coord, dims = 2)[:, 1]
    range = coordmax - coordmin
    coordmin -= range * eps(eltype(coord))
    coordmax += range * eps(eltype(coord))

    # Calculate cell size for averaging
    ilenmax = NPD ./ range


    # Now number of grid points in each direction
    sz = (round.(Int, (coordmax - coordmin) .* ilenmax) .+ 1...,)::NTuple{ndim,Int}
    szb = (ndim, sz...)


    # now allocate the arrays

    if weights == []
        weights = ones(Float64, ndata)
    end

    # Working arrays
    VP = zeros(Float64, sz)
    VARP = zeros(Float64, sz)
    XP = zeros(Float64, szb)
    WP = zeros(Float64, sz)
    NG = zeros(Int, ndim)
    gridindex = zeros(Int, ndata, ndim)


    for i = 1:ndata
        for j = 1:ndim
            # array index of the point
            NG[j] = round(Int, (coord[j, i] - coordmin[j]) * ilenmax[j]) + 1

        end
        for j = 1:ndim
            # summing coordinates
            XP[j, NG...] += weights[i] * coord[j, i]

        end
        # summing values
        VP[NG...] += weights[i] * val[i]
        # prepares variance
        VARP[NG...] += weights[i] * val[i] * val[i]
        # bookkeeping of weights
        WP[NG...] += weights[i]
    end

    # Now just scale as needed and put out grid points with non-zero WP

    # Already here take out the grid and keep only points where WP.>0
    #
    sel = WP .> 0
    sza = sum(sel)


    # For arrays easy to extract relevant points
    VPa = VP[sel]
    VARPa = VARP[sel]
    WPa = WP[sel]

    #For coordinates a little bit more tricky
    XPa = zeros(Float64, (ndim, sza))
    for j = 1:ndim
        work = XP[j, [(:) for jj = 1:ndim]...]
        XPa[j, :] = work[sel]
    end

    # Now simple average and variance calculation if "intensive variable"
    if intensive
        VPa = VPa ./ WPa
        VARPa = VARPa ./ WPa .- VPa.^2
    end

    # Coordinates are always average values
    for j = 1:ndim
        XPa[j, :] = XPa[j, :] ./ WPa
    end
    # Finished, return coordinate tuple and arrays
    return ([XPa[j, :] for j = 1:ndim]...,), VPa, WPa, VARPa, ilenmax
end