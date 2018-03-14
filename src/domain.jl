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
            mask[:,:,k] = bi + depthr[k] .< 0
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



