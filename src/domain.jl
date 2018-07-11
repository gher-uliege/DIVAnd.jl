
"""
    Δcoord = localresolution(coord)

Estimate the local resolution of the coordinate `coord` by finite differences.
"""
function localresolution(coord)
    Δcoord = similar(coord)

    if length(coord) == 1
        error("coord has just one element $(coord)")
    end

    for i = 2:length(coord)-1
        Δcoord[i] = (coord[i+1] - coord[i-1])/2
    end
    Δcoord[1] = coord[2] - coord[1]
    Δcoord[end] = coord[end] - coord[end-1]

    return Δcoord
end


"""
    mask,pmn,xyi = DIVAnd_squaredom(n,coord)

Create a "square" domain in `n` dimensions with the coordinates `coord`
assuming a Catersian metric. This functions returns
the mask `mask`, the coordinates `(xi,yi,...)` and the metric `(pm,pn...)`.

# Example

mask,(pm,pn),(xi,yi) = DIVAnd_squaredom(2,range(0,stop=1,length=50))
"""
function DIVAnd_squaredom(n,coord)
    coords = ntuple(i-> coord,n)
    return DIVAnd_rectdom(coords...)
end


"""
    mask,pmn,xyi = DIVAnd_rectdom(coord1,coord2,...)

Create a "rectangular" domain in `n` dimensions with the coordinates `coord1`
`coord2`... assuming a Catersian metric. This functions returns
the mask `mask`, the coordinates `(xi,yi,...)` and the metric `(pm,pn...)`.

For example:

```julia-repl
julia> mask,(pm,pn),(xi,yi) = DIVAnd_rectdom(range(0,stop=1,length=50),linspace(0,stop=1,length=50))
```
"""
function DIVAnd_rectdom(coords...)
    # grid of background field
    xyi = ndgrid(coords...)

    # mask (all points are valid)
    mask = trues(size(xyi[1]))

    # metric (inverse of the resolution)
    pmn = ndgrid([1 ./ localresolution(coords[i]) for i = 1:length(coords)]...)

    return mask,pmn,xyi
end



"""
    mask,(pm,pn),(xi,yi) = domain(bathname,bathisglobal,lonr,latr)

Generate a 2D geospatial domain based on the topography from the NetCDF file
`bathname`.
"""
function domain(bathname,bathisglobal,lonr,latr)
    mask,(pm,pn),(xi,yi) = DIVAnd.DIVAnd_rectdom(lonr,latr)

    mxi,myi,mask = DIVAnd.load_mask(bathname,bathisglobal,lonr,latr,0.)

    pm,pn = DIVAnd.DIVAnd_metric(xi,yi)

    return mask,(pm,pn),(xi,yi)
end

"""
    mask,(pm,pn,po),(xi,yi,zi) = domain(bathname,bathisglobal,lonr,latr,depthr)

Generate a 3D geospatial domain based on the topography from the NetCDF file
`bathname`. If `zlevel` is `:surface`, then `depthr` is zero for the sea surface and 
positive in water (positive is down). If `zlevel` is `:floor`, then `depthr` is 
zero for the sea floor and positive in water (positive is up)
"""
function domain(bathname,bathisglobal,lonr,latr,depthr; zlevel = :surface)

    mask,(pm,pn,po),(xi,yi,zi) = DIVAnd.DIVAnd_rectdom(lonr,latr,depthr)

    pm[:,:,1],pn[:,:,1] = DIVAnd.DIVAnd_metric(xi[:,:,1,1],yi[:,:,1,1])
    for k = 1:size(pm,3)
        pm[:,:,k] = pm[:,:,1]
        pn[:,:,k] = pn[:,:,1]
    end

    dx = lonr[2] - lonr[1]
    dy = latr[2] - latr[1]

    if zlevel == :surface
        mxi,myi,mask[:,:,:] = DIVAnd.load_mask(bathname,bathisglobal,lonr,latr,depthr)
    else        
        bxi,byi,bi = load_bath(bathname,bathisglobal,lonr,latr)
        
        for k = 1:size(mask,3)
            mask[:,:,k] = depthr[k] .< bi
        end
    end

    return mask,(pm,pn,po),(xi,yi,zi)
end


"""
    mask,(pm,pn,po,pp),(xi,yi,zi,ti) = domain(bathname,bathisglobal,lonr,latr,depthr,timer)

Generate a geospatial domain based on the topography from the NetCDF file
`bathname`.
"""
function domain(bathname,bathisglobal,lonr,latr,depthr,timer)

    mask,(pm,pn,po,pp),(xi,yi,zi,ti) = DIVAnd.DIVAnd_rectdom(lonr,latr,depthr,timer)

    pm[:,:,1,1],pn[:,:,1,1] = DIVAnd.DIVAnd_metric(xi[:,:,1,1],yi[:,:,1,1])
    for n = 1:size(pn,4)
        for k = 1:size(pm,3)
            pm[:,:,k,n] = pm[:,:,1,1]
            pn[:,:,k,n] = pn[:,:,1,1]
        end
    end
    #pm[:] = 4.15379e-5
    #pn[:] = 4.15379e-5

    dx = lonr[2] - lonr[1]
    dy = latr[2] - latr[1]

    mxi,myi,mask2 = DIVAnd.load_mask(bathname,bathisglobal,lonr,latr,depthr)

    mask3 = repeat(mask2,inner = (1,1,1,length(timer)))

    return mask3,(pm,pn,po,pp),(xi,yi,zi,ti)
end



