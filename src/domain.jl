"""
mask,xyi,pmn = divand_squaredom(n,coord)

Create a "square" domain in `n` dimensions with the coordinates `coord`
assuming a Catersian metric. This functions returns
the mask `mask`, the coordinates `(xi,yi,...)` and the metric `(pm,pn...)`.

# Example

mask,(pm,pn),(xi,yi) = divand_squaredom(2,linspace(0,1,50))
"""
function divand_squaredom(n,coord)
    coords = ([coord for i = 1:n]...)
    return divand_rectdom(coords...)
end


"""
mask,xyi,pmn = divand_squaredom(n,coord)

Create a "square" domain in `n` dimensions with the coordinates `coord`
assuming a Catersian metric. This functions returns
the mask `mask`, the coordinates `(xi,yi,...)` and the metric `(pm,pn...)`.

# Example

mask,(pm,pn),(xi,yi) = divand_rectdom(linspace(0,1,50),linspace(0,1,50))
"""
function divand_rectdom(coords...)
    # grid of background field
    xyi = ndgrid(coords...)

    # mask (all points are valid)
    mask = trues(xyi[1])

    # metric (inverse of the resolution)
    pmn = ([ones(size(mask)) / (coords[i][2]-coords[i][1]) for i = 1:length(coords)]...)

    return mask,pmn,xyi
end



"""
    mask,(pm,pn),(xi,yi) = domain(bathname,bathisglobal,lonr,latr)

Generate a 2D geospatial domain based on the topography from the NetCDF file
`bathname`.
"""

function domain(bathname,bathisglobal,lonr,latr)
    mask,(pm,pn),(xi,yi) = divand.divand_rectdom(lonr,latr)

    mxi,myi,mask = divand.load_mask(bathname,bathisglobal,lonr,latr,0.)

    pm,pn = divand.divand_metric(xi,yi)

    return mask,(pm,pn),(xi,yi)
end

"""
    mask,(pm,pn,po),(xi,yi,zi) = domain(bathname,bathisglobal,lonr,latr,depthr)

Generate a 3D geospatial domain based on the topography from the NetCDF file
`bathname`.
if zlevel = :surface, then depthr is zero for the sea surface and positive in water (positive is down)
if zlevel = :floor, then depthr is zero for the sea floor and positive in water (positive is up)

"""

function domain(bathname,bathisglobal,lonr,latr,depthr; zlevel = :surface)

    mask,(pm,pn,po),(xi,yi,zi) = divand.divand_rectdom(lonr,latr,depthr)

    pm[:,:,1],pn[:,:,1] = divand.divand_metric(xi[:,:,1,1],yi[:,:,1,1])
    for k = 1:size(pm,3)
        pm[:,:,k] = pm[:,:,1]
        pn[:,:,k] = pn[:,:,1]
    end

    dx = lonr[2] - lonr[1]
    dy = latr[2] - latr[1]

    if zlevel == :surface
        mxi,myi,mask[:,:,:] = divand.load_mask(bathname,bathisglobal,lonr,latr,depthr)
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

    mask,(pm,pn,po,pp),(xi,yi,zi,ti) = divand.divand_rectdom(lonr,latr,depthr,timer)

    pm[:,:,1,1],pn[:,:,1,1] = divand.divand_metric(xi[:,:,1,1],yi[:,:,1,1])
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

    mxi,myi,mask2 = divand.load_mask(bathname,bathisglobal,lonr,latr,depthr)

    mask3 = repeat(mask2,inner = (1,1,1,length(timer)))

    return mask3,(pm,pn,po,pp),(xi,yi,zi,ti)
end



