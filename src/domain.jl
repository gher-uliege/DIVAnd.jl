"""
    mask,(pm,pn),(xi,yi) = domain(bathname,isglobal,lonr,latr)

Generate a 2D geospatial domain based on the topography from the NetCDF file
`bathname`.
"""

function domain(bathname,isglobal,lonr,latr)
    mask,(pm,pn),(xi,yi) = divand.divand_rectdom(lonr,latr)

    mxi,myi,mask = divand.load_mask(bathname,isglobal,lonr,latr,1.)

    pm,pn = divand.divand_metric(xi,yi)

    return mask,(pm,pn),(xi,yi)
end


"""
    mask,(pm,pn,po,pp),(xi,yi,zi,ti) = domain(bathname,isglobal,lonr,latr,depthr,timer)

Generate a geospatial domain based on the topography from the NetCDF file
`bathname`.
"""

function domain(bathname,isglobal,lonr,latr,depthr,timer)

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


    @show size(mask)
    dx = lonr[2] - lonr[1]
    dy = latr[2] - latr[1]

    mxi,myi,mask2 = divand.load_mask(bathname,isglobal,minimum(lonr),maximum(lonr),dx,minimum(latr),maximum(latr),dy,depthr)

    mask3 = repeat(mask2,inner = (1,1,1,length(timer)))

    return mask3,(pm,pn,po,pp),(xi,yi,zi,ti)
end



